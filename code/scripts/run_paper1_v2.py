"""Paper 1 (v2) rerun: frozen QEC-SB confirmation + paired CP-SA and MPS baselines.

Everything runs on the *same* 12 confirmation batches (d=5,7,9; p=.10,.12; two
replicates each; 12,000 shots) produced by ``build_confirmation_plan``. Sample
and solver seeds are those of the original confirmation, so the frozen decoder
columns must reproduce the published table (checked in ``analyze``).

Subcommands
-----------
  check      code-mapping self test (ours <-> qecsim), d=3..9; seconds
  mps-d3     MPS vs exact maximum-likelihood decoding at d=3 (validates MPS as MLD reference)
  ours       rerun frozen P4.2a pipeline (MWPM, CP-SB MAP, local SOE, importance, Bennett)
  baselines  PyMatching, CP-SA MAP / CP-SA local SOE, MPS (chi list) on the same shots
  analyze    merge, reproduction check, LER table, paired tests, MPS gap closure
  all        check -> mps-d3 -> ours + baselines (one shared process pool) -> analyze

Example (Linux server):
  python scripts/run_paper1_v2.py all --out results/paper1_v2 --workers 32
"""

from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per worker process: essential on many-core servers

import argparse
import hashlib
import json
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from qec_sb.adaptive_soe import fixed_soe_prediction  # noqa: E402
from qec_sb.scalable import collect_cp_class_pools_no_oracle  # noqa: E402
from qec_sb.cp_sa import CPSAParams, run_constraint_preserving_sa  # noqa: E402
from qec_sb.cp_sb import build_stabilizer_hamiltonian  # noqa: E402
from qec_sb.logical_classes import logical_class_seed  # noqa: E402
from qec_sb.p42a import _syn_key, run_p42a  # noqa: E402
from qec_sb.p5.baselines import PyMatchingCSSDecoder  # noqa: E402
from qec_sb.p5.config import build_frozen_p42a_config, frozen_config_sha256  # noqa: E402
from qec_sb.p5.confirmation import DEFAULT_BASE_SEED, build_confirmation_plan  # noqa: E402
from qec_sb.p5.statistics import adjust_p_values, exact_paired_statistics, wilson_interval  # noqa: E402
from qec_sb.rotated_code import build_rotated_surface_code, sample_depolarizing  # noqa: E402
from qec_sb.scalable import CSSMWPMDecoder, correction_class, cp_map_prediction  # noqa: E402

# Published confirmation (paper/RESULTS.md): (d, p) -> (final QEC-SB LER, MWPM LER), pooled shots.
PUBLISHED = {
    (5, .10): (.0633, .0930), (5, .12): (.1203, .1543),
    (7, .10): (.0355, .0775), (7, .12): (.0805, .1190),
    (9, .10): (.0270, .0640), (9, .12): (.0520, .1300),
}
OURS_METHODS = ("MWPM", "CP_SB_MAP", "CP_local_SOE", "P4_1_importance", "P4_2_Bennett")
SA_PARAMS = CPSAParams(chains=8, sweeps=100, temperature_start=2.0, temperature_end=0.05,
                       restart_modes=("warm", "random"))  # pilot main budget: 16 chains/class
SA_SEED_OFFSET = 5_000_011  # own namespace, disjoint from SB (+0), importance (+7000003), Bennett (+9000007)
D3_SEED = 2027201300        # separate namespace for the d=3 MPS-vs-exact check


# --------------------------------------------------------------------------- utils
def gf2_rank(m):
    a = (np.asarray(m, dtype=np.uint8) & 1).copy()
    rank, rows, cols = 0, a.shape[0], a.shape[1]
    for c in range(cols):
        pivot = np.flatnonzero(a[rank:, c])
        if not len(pivot):
            continue
        r = rank + pivot[0]
        a[[rank, r]] = a[[r, rank]]
        hit = np.flatnonzero(a[:, c])
        hit = hit[hit != rank]
        a[hit] ^= a[rank]
        rank += 1
        if rank == rows:
            break
    return rank


def same_rowspace(a, b):
    ra, rb = gf2_rank(a), gf2_rank(b)
    return ra == rb == gf2_rank(np.vstack([a, b]))


def sample_sha256(samples):
    h = hashlib.sha256()
    for key in ("pauli", "sx", "sz", "logical"):
        h.update(np.ascontiguousarray(samples[key]).tobytes())
    return h.hexdigest()


def env_manifest():
    import numpy, scipy, pymatching, stim, qecsim  # noqa: E401
    return {"python": sys.version.split()[0], "platform": platform.platform(),
            "numpy": numpy.__version__, "scipy": scipy.__version__, "pandas": pd.__version__,
            "pymatching": pymatching.__version__, "stim": stim.__version__, "qecsim": qecsim.__version__,
            "cpu_count": os.cpu_count(), "frozen_config_sha256": frozen_config_sha256()}


# --------------------------------------------------------- ours <-> qecsim mapping
class QecsimBridge:
    """Bijection between our rotated code and qecsim's RotatedPlanarCode.

    Found by search over the 8 lattice symmetries and an optional X<->Z swap and
    accepted only if the two stabilizer groups coincide exactly (GF(2) row
    space). Logical classes are always evaluated in *our* convention after
    mapping a recovery back, so only stabilizer-group equality is required.
    """

    def __init__(self, code):
        from qecsim.models.rotatedplanar import RotatedPlanarCode
        d = code.distance
        self.code, self.n = code, code.n
        self.qcode = RotatedPlanarCode(d, d)
        qidx = {}
        for x in range(d):
            for y in range(d):
                qidx[(x, y)] = int(np.flatnonzero(self.qcode.new_pauli().site("X", (x, y)).to_bsf())[0])
        ours_xy = [((cx - 1) // 2, (cy - 1) // 2) for cx, cy in code.data_coords]
        n = code.n
        ours_stab = np.vstack([
            np.hstack([code.Hx, np.zeros_like(code.Hx)]),   # X-type checks
            np.hstack([np.zeros_like(code.Hz), code.Hz]),   # Z-type checks
        ])
        transforms = [
            lambda x, y: (x, y), lambda x, y: (d - 1 - x, y), lambda x, y: (x, d - 1 - y),
            lambda x, y: (d - 1 - x, d - 1 - y), lambda x, y: (y, x), lambda x, y: (d - 1 - y, x),
            lambda x, y: (y, d - 1 - x), lambda x, y: (d - 1 - y, d - 1 - x),
        ]
        self.perm = self.swap = None
        for t in transforms:
            perm = np.asarray([qidx[t(*xy)] for xy in ours_xy])
            for swap in (False, True):
                mapped = self._to_q(ours_stab[:, :n], ours_stab[:, n:], perm, swap)
                if same_rowspace(mapped, self.qcode.stabilizers):
                    self.perm, self.swap = perm, swap
                    break
            if self.perm is not None:
                break
        if self.perm is None:
            raise RuntimeError(f"no stabilizer-preserving map to qecsim found for d={d}")
        self.inv = np.argsort(self.perm)

    def _to_q(self, ex, ez, perm=None, swap=None):
        perm = self.perm if perm is None else perm
        swap = self.swap if swap is None else swap
        ex, ez = np.atleast_2d(ex), np.atleast_2d(ez)
        qx = np.zeros_like(ex); qz = np.zeros_like(ez)
        qx[:, perm] = ez if swap else ex
        qz[:, perm] = ex if swap else ez
        return np.hstack([qx, qz]).astype(np.uint8)

    def to_qecsim(self, ex, ez):
        return self._to_q(ex, ez)[0]

    def from_qecsim(self, bsf):
        bsf = np.asarray(bsf, dtype=np.uint8)
        qx, qz = bsf[: self.n][self.perm], bsf[self.n:][self.perm]
        return (qz, qx) if self.swap else (qx, qz)

    def qsyndrome(self, ex, ez):
        from qecsim import paulitools as pt
        return pt.bsp(self.to_qecsim(ex, ez), self.qcode.stabilizers.T)


class MPSDecoder:
    def __init__(self, code, p, chi):
        from qecsim.models.generic import DepolarizingErrorModel
        from qecsim.models.rotatedplanar import RotatedPlanarMPSDecoder
        self.bridge, self.code, self.p = QecsimBridge(code), code, float(p)
        self.dec, self.em = RotatedPlanarMPSDecoder(chi=int(chi)), DepolarizingErrorModel()

    def decode_class(self, ex_feasible, ez_feasible, sx, sz):
        """Decode from the syndrome only (a feasible correction just carries it across)."""
        syn = self.bridge.qsyndrome(ex_feasible, ez_feasible)
        rec = self.dec.decode(self.bridge.qcode, syn, error_model=self.em, error_probability=self.p)
        rx, rz = self.bridge.from_qecsim(rec)
        if not (np.array_equal((rx @ self.code.Hz.T) & 1, sx) and np.array_equal((rz @ self.code.Hx.T) & 1, sz)):
            raise RuntimeError("MPS recovery does not reproduce the syndrome after mapping")
        return int(correction_class(self.code, rx, rz))


# ------------------------------------------------------------------- CP-SA solver
def cp_sa_predictions(code, ex0, ez0, p, seed):
    hams, pools, proposals = {}, {}, 0
    for c in range(4):
        cex, cez = logical_class_seed(code, ex0, ez0, c)
        ham = build_stabilizer_hamiltonian(code, cex, cez)
        res = run_constraint_preserving_sa(ham, SA_PARAMS, int(seed + c * 100003))
        u = np.unique(np.vstack([np.zeros((1, ham.n_spins), np.uint8), res.terminal_bits]), axis=0)
        wt = ham.energy_bits(u)
        pools[c], hams[c] = u[np.argsort(wt, kind="stable")], ham
        proposals += res.proposals
    map_cls, _ = cp_map_prediction(pools, hams)
    soe_cls, *_ = fixed_soe_prediction(pools, hams, p, 1, 1)
    return int(map_cls), int(soe_cls), proposals


# -------------------------------------------------------------------------- tasks
def task_ours(out, batch):
    bdir = Path(out) / "ours" / batch["batch_id"]
    if (bdir / "metadata.json").exists():
        return f"ours {batch['batch_id']} cached"
    tmp = bdir.with_name(bdir.name + ".tmp")
    t0 = time.time()
    run_p42a(tmp, build_frozen_p42a_config(batch["distance"], batch["p"], batch["expected_shots"],
                                           batch["seed"], batch["solver_seed"]))
    if bdir.exists():
        raise RuntimeError(f"{bdir} appeared during run")
    tmp.rename(bdir)
    return f"ours {batch['batch_id']} {time.time() - t0:.0f}s"


def import_bposd():
    """ldpc's package __init__ imports sinter (only needed for its sinter wrappers)."""
    try:
        from ldpc import BpOsdDecoder
    except ModuleNotFoundError as exc:
        if exc.name != "sinter":
            raise
        import types
        stub = types.ModuleType("sinter"); stub.Decoder = object; stub.CompiledDecoder = object
        sys.modules["sinter"] = stub
        from ldpc import BpOsdDecoder
    return BpOsdDecoder


class BPOSDDecoder:
    """BP+OSD-CS (order 7, min-sum, max_iter=n) on the X and Z components, q=2p/3 each."""

    def __init__(self, code, p):
        import scipy.sparse as sp
        BpOsd = import_bposd()
        q = 2.0 * float(p) / 3.0
        kw = dict(error_rate=q, max_iter=code.n, bp_method="minimum_sum", osd_method="osd_cs", osd_order=7)
        self.code = code
        self.xdec = BpOsd(sp.csr_matrix(code.Hz), **kw)
        self.zdec = BpOsd(sp.csr_matrix(code.Hx), **kw)

    def decode(self, sx, sz):
        ex = np.asarray(self.xdec.decode(np.asarray(sx, np.uint8)), np.uint8)
        ez = np.asarray(self.zdec.decode(np.asarray(sz, np.uint8)), np.uint8)
        if not (np.array_equal((ex @ self.code.Hz.T) & 1, sx) and np.array_equal((ez @ self.code.Hx.T) & 1, sz)):
            raise RuntimeError("BP+OSD returned a syndrome-incompatible correction")
        return ex, ez


class MethodRunner:
    """One baseline/ablation method; decode(sx, sz, key, solver_seed) -> dict of columns.

    pym    PyMatching (independent X/Z)          bposd  BP+OSD-CS (independent X/Z)
    sa     CP-SA MAP and CP-SA local SOE         mpsN   MPS decoder, chi=N
    r2     CP-SB pools rebuilt exactly as in the frozen pipeline; local SOE radius 1 (check) and 2
    mld    exact maximum likelihood (d=3 only)
    """

    def __init__(self, code, p, method):
        self.code, self.p, self.method = code, float(p), method
        self.mwpm = CSSMWPMDecoder(code)
        if method == "pym":
            self.dec = PyMatchingCSSDecoder(code, p)
        elif method == "bposd":
            self.dec = BPOSDDecoder(code, p)
        elif method.startswith("mps"):
            self.dec = MPSDecoder(code, p, int(method[3:]))
        elif method == "r2":
            self.cp_cfg = build_frozen_p42a_config(code.distance, p, 1, 0, 0).cp
        elif method == "mld":
            from qec_sb.oracle import build_d3_oracle
            self.oracle = build_d3_oracle(code)
        elif method != "sa":
            raise ValueError(f"unknown method {method}")

    def decode(self, sx, sz, key, solver_seed):
        code, p, m = self.code, self.p, self.method
        if m in ("pym", "bposd"):
            ex, ez = self.dec.decode(sx, sz)
            return {"pred_PyMatching" if m == "pym" else "pred_BP_OSD": int(correction_class(code, ex, ez))}
        if m == "mld":
            _, masses = self.oracle.exact_mld_class(key, p)
            return {"pred_MLD_exact": int(np.argmax(masses))}
        ex0, ez0 = self.mwpm.decode(sx, sz)
        if m == "sa":
            a, b, props = cp_sa_predictions(code, ex0, ez0, p, solver_seed + SA_SEED_OFFSET + key * 1009)
            return {"pred_CP_SA_MAP": a, "pred_CP_SA_local_SOE": b, "sa_proposals": props}
        if m == "r2":
            pools, hams, _ = collect_cp_class_pools_no_oracle(code, ex0, ez0, solver_seed + key * 1009, self.cp_cfg)
            cmap, _ = cp_map_prediction(pools, hams)
            r1, *_ = fixed_soe_prediction(pools, hams, p, 1, 1)
            r2, *_ = fixed_soe_prediction(pools, hams, p, 1, 2)
            return {"pred_CP_SB_MAP_rebuilt": int(cmap), "pred_CP_local_SOE_rebuilt": int(r1),
                    "pred_CP_local_SOE_r2": int(r2)}
        return {f"pred_MPS_chi{m[3:]}": self.dec.decode_class(ex0, ez0, sx, sz)}


def decode_samples(samples, lo, hi, runner, solver_seed, extra):
    """Decode shots [lo, hi) with per-syndrome caching; returns one row per shot."""
    cache, rows = {}, []
    for i in range(lo, hi):
        sx, sz = samples["sx"][i], samples["sz"][i]
        key = _syn_key(sx, sz)
        if key not in cache:
            ts = time.perf_counter()
            pred = runner.decode(sx, sz, key, solver_seed)
            pred[f"time_{runner.method}_s"] = time.perf_counter() - ts
            cache[key] = pred
        truth = int(samples["logical"][i, 0]) * 2 + int(samples["logical"][i, 1])
        rows.append({**extra, "shot": i, "syndrome_key": str(key), "truth_class": truth,
                     "error_weight": int(np.count_nonzero(samples["pauli"][i])), **cache[key]})
    return rows


def write_atomic(df, path):
    tmp = Path(str(path) + ".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(path)


def task_baseline(out, batch, method, lo, hi, chis):
    path = Path(out) / "baselines" / f"{batch['batch_id']}__{method}__{lo:05d}.csv"
    if path.exists():
        return f"{path.name} cached"
    code = build_rotated_surface_code(batch["distance"])
    samples = sample_depolarizing(code, batch["p"], batch["expected_shots"], batch["seed"])
    t0 = time.time()
    runner = MethodRunner(code, batch["p"], method)
    rows = decode_samples(samples, lo, hi, runner, batch["solver_seed"],
                          {"batch_id": batch["batch_id"], "sample_sha256": sample_sha256(samples)})
    write_atomic(pd.DataFrame(rows), path)
    return f"{path.name} {time.time() - t0:.0f}s"


SMOKE = 0  # >0: truncate every batch to this many shots (pipeline test only)
BASE_SEED = DEFAULT_BASE_SEED  # --base-seed: seed family of the confirmation batches


def get_plan():
    plan = [asdict(b) for b in build_confirmation_plan(BASE_SEED)]
    if SMOKE:
        for b in plan:
            b["expected_shots"] = min(SMOKE, b["expected_shots"])
    return plan


def build_tasks(out, chis, chunk, which):
    plan = get_plan()
    tasks = []
    if "ours" in which:
        tasks += [(task_ours, (out, b)) for b in plan]
    if "baselines" in which:
        methods = ["pym", "bposd", "sa", "r2"] + [f"mps{c}" for c in chis]
        for b in plan:
            for m in methods:
                for lo in range(0, b["expected_shots"], chunk):
                    tasks.append((task_baseline, (out, b, m, lo, min(lo + chunk, b["expected_shots"]), chis)))
    # longest first: d=9 SA chunks and the ours batches
    tasks.sort(key=lambda t: -(t[1][1]["distance"] ** 3) * (3 if t[0] is task_ours else 1))
    return tasks


def run_pool(tasks, workers):
    failures = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, *args): (fn.__name__, args[1]["batch_id"]) for fn, args in tasks}
        for k, f in enumerate(as_completed(futs), 1):
            try:
                print(f"[{k}/{len(futs)}] {f.result()}", flush=True)
            except Exception as e:  # keep going, report at the end
                failures.append((futs[f], repr(e)))
                print(f"[{k}/{len(futs)}] FAILED {futs[f]}: {e!r}", flush=True)
    if failures:
        raise SystemExit(f"{len(failures)} task(s) failed: {failures[:5]}")


# ------------------------------------------------------------------- self checks
def cmd_check(out):
    for d in (3, 5, 7, 9):
        code = build_rotated_surface_code(d)
        br = QecsimBridge(code)
        s = sample_depolarizing(code, .12, 200, 99)
        mw = CSSMWPMDecoder(code)
        for i in range(200):
            ex0, ez0 = mw.decode(s["sx"][i], s["sz"][i])
            assert np.array_equal(br.qsyndrome(s["ex"][i], s["ez"][i]), br.qsyndrome(ex0, ez0))
            rx, rz = br.from_qecsim(br.to_qecsim(s["ex"][i], s["ez"][i]))
            assert np.array_equal(rx, s["ex"][i]) and np.array_equal(rz, s["ez"][i])
        print(f"check d={d}: stabilizer groups match (swap={br.swap}); 200-shot syndrome/round-trip OK")


def cmd_mps_d3(out, shots, chis):
    from qec_sb.oracle import build_d3_oracle
    path = Path(out) / "mps_d3_check.csv"
    code = build_rotated_surface_code(3)
    oracle = build_d3_oracle(code)
    mw = CSSMWPMDecoder(code)
    rows = []
    for j, p in enumerate((.10, .12)):
        s = sample_depolarizing(code, p, shots, D3_SEED + j)
        decs = {c: MPSDecoder(code, p, c) for c in chis}
        cache = {}
        for i in range(shots):
            sx, sz = s["sx"][i], s["sz"][i]
            key = _syn_key(sx, sz)
            if key not in cache:
                _, masses = oracle.exact_mld_class(key, p)
                top = np.sort(masses)[::-1]
                ex0, ez0 = mw.decode(sx, sz)
                cache[key] = {"pred_MLD": int(np.argmax(masses)),
                              "mld_tie": int(np.isclose(top[0], top[1], rtol=1e-12, atol=0)),
                              **{f"pred_MPS_chi{c}": decs[c].decode_class(ex0, ez0, sx, sz) for c in chis}}
            truth = int(s["logical"][i, 0]) * 2 + int(s["logical"][i, 1])
            rows.append({"p": p, "shot": i, "truth_class": truth, **cache[key]})
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    lines = ["# d=3: MPS vs exact maximum-likelihood decoding", "",
             "| p | shots | exact-MLD LER | " + " | ".join(f"MPS chi={c} LER | agree w/ MLD (non-tie shots)" for c in chis) + " |",
             "|---|---|---|" + "---|---|" * len(chis)]
    for p, g in df.groupby("p"):
        nt = g[g.mld_tie == 0]
        cells = [f"{(g.pred_MLD != g.truth_class).mean():.4%}"]
        for c in chis:
            cells += [f"{(g[f'pred_MPS_chi{c}'] != g.truth_class).mean():.4%}",
                      f"{(nt[f'pred_MPS_chi{c}'] == nt.pred_MLD).mean():.4%} ({len(nt)})"]
        lines.append(f"| {p:.2f} | {len(g)} | " + " | ".join(cells) + " |")
    (Path(out) / "MPS_D3_CHECK.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


# ----------------------------------------------------------------------- analysis
def load_merged(out):
    plan = get_plan()
    frames = []
    for b in plan:
        o = pd.read_csv(Path(out) / "ours" / b["batch_id"] / "shot_level.csv", dtype={"syndrome_key": str})
        o["batch_id"], o["distance"], o["p"] = b["batch_id"], b["distance"], b["p"]
        smp = sample_depolarizing(build_rotated_surface_code(b["distance"]), b["p"], b["expected_shots"], b["seed"])
        o["error_weight"] = np.count_nonzero(smp["pauli"], axis=1)
        parts = sorted((Path(out) / "baselines").glob(f"{b['batch_id']}__*.csv"))
        by_method = {}
        for f in parts:
            by_method.setdefault(f.name.split("__")[1], []).append(pd.read_csv(f))
        for m, dfs in by_method.items():
            bd = pd.concat(dfs).sort_values("shot").reset_index(drop=True)
            if bd.shot.tolist() != list(range(b["expected_shots"])):
                raise ValueError(f"{b['batch_id']} {m}: incomplete shot coverage")
            if bd.sample_sha256.nunique() != 1:
                raise ValueError(f"{b['batch_id']} {m}: sample hash mismatch across chunks")
            if not np.array_equal(bd.truth_class.to_numpy(), o.truth_class.to_numpy()):
                raise ValueError(f"{b['batch_id']} {m}: truth differs from frozen-decoder rerun")
            keep = [c for c in bd.columns if c.startswith(("pred_", "time_", "sa_"))]
            o = pd.concat([o, bd[keep]], axis=1)
        frames.append(o)
    df = pd.concat(frames, ignore_index=True)
    for c in [c for c in df.columns if c.startswith("pred_")]:
        df["correct_" + c[5:]] = (df[c] == df.truth_class).astype(np.uint8)
    return df


def gap_closure(mwpm_ok, x_ok, ref_ok, reps=4000, seed=1):
    """(LER_MWPM - LER_X) / (LER_MWPM - LER_ref) with paired shot bootstrap."""
    def g(a, b, c):
        den = c.mean() - a.mean()
        return np.nan if den <= 0 else (b.mean() - a.mean()) / den
    rng = np.random.default_rng(seed)
    n = len(mwpm_ok)
    boots = []
    for _ in range(reps):
        idx = rng.integers(0, n, n)
        boots.append(g(mwpm_ok[idx], x_ok[idx], ref_ok[idx]))
    boots = np.asarray(boots)
    boots = boots[np.isfinite(boots)]
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) else (np.nan, np.nan))
    return g(mwpm_ok, x_ok, ref_ok), lo, hi


MAIN_CANDIDATES = ("CP_local_SOE", "CP_local_SOE_r2")


def select_main_method(df, alpha=0.05):
    """Pre-declared rule (fixed before any grid data exist).

    The paper's main decoder is CP-SB + local SOE. Radius 2 replaces radius 1 only if,
    pooled over all 12,000 confirmation shots, it has strictly fewer failures and the
    exact paired McNemar test gives p < alpha. Importance / Bennett are ablations and
    are never eligible as the main method.
    """
    r1, r2 = (f"correct_{m}" for m in MAIN_CANDIDATES)
    if r2 not in df:
        return {"main": MAIN_CANDIDATES[0], "reason": "r2 not available"}
    st = exact_paired_statistics(df[r2].to_numpy(), df[r1].to_numpy())
    better = st["discordant_wins"] > st["discordant_losses"] and st["exact_p_value"] < alpha
    return {"main": MAIN_CANDIDATES[1] if better else MAIN_CANDIDATES[0],
            "rule": select_main_method.__doc__.strip(), "r2_wins": st["discordant_wins"],
            "r2_losses": st["discordant_losses"], "exact_p": st["exact_p_value"], "shots": st["shots"]}


def effective_distance_table(df, methods, dcol):
    """Smallest physical error weight that caused a logical failure, per cell and method."""
    rows = []
    for (d, p), g in df.groupby([dcol, "p"]):
        half = (int(d) + 1) // 2
        for m in methods:
            c = g[f"correct_{m}"]
            fail = g.loc[c.notna() & (c == 0), "error_weight"]
            rows.append({"d": d, "p": p, "method": m, "failures": len(fail),
                         "min_failing_weight": int(fail.min()) if len(fail) else None,
                         "failures_below_ceil_d_half": int((fail < half).sum()), "ceil_d_half": half})
    return pd.DataFrame(rows)


def per_unique_syndrome_time(df, keys):
    """Mean wall time per unique syndrome (each unique syndrome counted once per task chunk)."""
    timecols = [c for c in df.columns if c.startswith("time_")]
    if not timecols or "syndrome_key" not in df:
        return pd.DataFrame()
    u = df.drop_duplicates(keys + ["syndrome_key"])
    return u.groupby(keys[:1])[timecols].mean().reset_index()


def cmd_analyze(out, chis):
    out = Path(out)
    df = load_merged(out)
    df.to_csv(out / "shot_level_merged.csv", index=False)
    ref = f"MPS_chi{max(chis)}"
    final = "P4_2_Bennett"           # frozen pipeline output, used only for the reproduction check
    sel = select_main_method(df)
    main = sel["main"]
    (out / "main_method.json").write_text(json.dumps(sel, indent=2), encoding="utf-8")
    methods = [m for m in ("MWPM", "PyMatching", "BP_OSD", "CP_SA_MAP", "CP_SB_MAP", "CP_SA_local_SOE", "CP_local_SOE",
                           "CP_local_SOE_r2", "P4_1_importance", final, *[f"MPS_chi{c}" for c in chis])
               if f"correct_{m}" in df]
    ler_rows, pair_rows, repro, gaps, chi_rows = [], [], [], [], []
    # comparison families; Holm is applied within each (a, b) family across the 6 cells
    pairs = [(main, "MWPM"), (main, "PyMatching"), (main, "BP_OSD"), (main, ref),
             ("CP_SA_local_SOE", "CP_local_SOE"), ("CP_SA_MAP", "CP_SB_MAP"),
             ("CP_local_SOE_r2", "CP_local_SOE"), ("CP_local_SOE", "CP_SB_MAP"),
             ("P4_1_importance", "CP_local_SOE"), (final, "P4_1_importance"),
             (final, "MWPM"), ("PyMatching", "MWPM"), ("BP_OSD", ref), ("CP_SA_local_SOE", ref)]
    pairs = list(dict.fromkeys(pairs))
    for (d, p), g in df.groupby(["distance", "p"]):
        n = len(g)
        for m in methods:
            k = int(n - g[f"correct_{m}"].sum())
            lo, hi = wilson_interval(k, n)
            ler_rows.append({"d": d, "p": p, "method": m, "shots": n, "failures": k,
                             "LER": k / n, "LER_ci95_lo": lo, "LER_ci95_hi": hi})
        for a, b in pairs:
            if f"correct_{a}" in g and f"correct_{b}" in g and a != b:
                st = exact_paired_statistics(g[f"correct_{a}"].to_numpy(), g[f"correct_{b}"].to_numpy())
                pair_rows.append({"d": d, "p": p, "a": a, "b": b,
                                  "LER_a_minus_LER_b": -st["paired_difference"],
                                  "ci95_lo": -st["paired_ci95_hi"], "ci95_hi": -st["paired_ci95_lo"],
                                  "a_right_b_wrong": st["discordant_wins"], "a_wrong_b_right": st["discordant_losses"],
                                  "exact_p": st["exact_p_value"]})
        pub = PUBLISHED[(d, round(p, 2))]
        got = (1 - g[f"correct_{final}"].mean(), 1 - g["correct_MWPM"].mean())
        repro.append({"d": d, "p": p, "published_final": pub[0], "rerun_final": got[0],
                      "published_MWPM": pub[1], "rerun_MWPM": got[1],
                      "match": bool(abs(pub[0] - got[0]) < 6e-5 and abs(pub[1] - got[1]) < 6e-5),
                      "r2_pools_rebuilt_identically": (bool((g.pred_CP_SB_MAP_rebuilt == g.pred_CP_SB_MAP).all()
                                                            and (g.pred_CP_local_SOE_rebuilt == g.pred_CP_local_SOE).all())
                                                       if "pred_CP_SB_MAP_rebuilt" in g else None)})
        if f"correct_{ref}" in g:
            for m in dict.fromkeys((main, final, "CP_local_SOE", "CP_local_SOE_r2", "CP_SA_local_SOE",
                                    "CP_SB_MAP", "CP_SA_MAP", "BP_OSD")):
                if f"correct_{m}" in g:
                    val, lo, hi = gap_closure(g.correct_MWPM.to_numpy(), g[f"correct_{m}"].to_numpy(),
                                              g[f"correct_{ref}"].to_numpy())
                    gaps.append({"d": d, "p": p, "method": m, "reference": ref,
                                 "gap_closed": val, "ci95_lo": lo, "ci95_hi": hi})
            if len(chis) > 1:
                lo_chi = f"pred_MPS_chi{min(chis)}"
                chi_rows.append({"d": d, "p": p, "chi_low": min(chis), "chi_high": max(chis),
                                 "class_disagreements": int((g[lo_chi] != g[f"pred_{ref}"]).sum()), "shots": n})
    pairs_df = pd.DataFrame(pair_rows)
    for (a, b), idx in pairs_df.groupby(["a", "b"]).groups.items():
        pairs_df.loc[idx, "holm_p"] = adjust_p_values(pairs_df.loc[idx, "exact_p"].tolist(), "holm")
    tables = {"ler_table": pd.DataFrame(ler_rows), "paired_comparisons": pairs_df,
              "reproduction_check": pd.DataFrame(repro), "gap_closure": pd.DataFrame(gaps),
              "mps_chi_convergence": pd.DataFrame(chi_rows),
              "effective_distance": effective_distance_table(df, methods, "distance"),
              "time_per_unique_syndrome_CONTENDED": per_unique_syndrome_time(df, ["distance", "batch_id"])}
    for name, t in tables.items():
        t.to_csv(out / f"{name}.csv", index=False)
    md = ["# Paper 1 v2 results", "",
          f"Main method (pre-declared rule): **{main}** — r2 vs r1 wins/losses {sel.get('r2_wins')}/{sel.get('r2_losses')}, "
          f"exact p={sel.get('exact_p')}.", f"MPS reference: {ref}. Frozen pipeline output: {final} (reproduction only).",
          "Timing below was measured under full server load; use `run_paper1_grid.py timing` for paper numbers.", ""]
    for name, t in tables.items():
        md += [f"## {name}", "", t.to_markdown(index=False, floatfmt=".4g") if len(t) else "(empty)", ""]
    (out / "RESULTS_V2.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md[:4]))
    print(f"main method -> {main}  (written to {out / 'main_method.json'})")


# ---------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "mps-d3", "ours", "baselines", "decode", "analyze", "all"])
    ap.add_argument("--shard", type=int, default=0, help="this job's index (Slurm array: $SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--nshards", type=int, default=1, help="number of array jobs sharing the decode tasks")
    ap.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED,
                    help="base seed of the confirmation family (the published post-correction family used its own)")
    ap.add_argument("--out", default="results/paper1_v2")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--chi", type=int, nargs="+", default=[16, 32])
    ap.add_argument("--chunk", type=int, default=100, help="shots per baseline task")
    ap.add_argument("--d3-shots", type=int, default=5000)
    ap.add_argument("--smoke", type=int, default=0,
                    help="pipeline test: only the first N shots of every batch (reproduction check will not match)")
    a = ap.parse_args()
    global SMOKE, BASE_SEED
    SMOKE, BASE_SEED = a.smoke, a.base_seed
    if SMOKE:
        a.d3_shots = min(a.d3_shots, 200)
    out = Path(a.out)
    for sub in ("", "ours", "baselines"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    manifest = out / "run_manifest.json"
    if not manifest.exists():
        manifest.write_text(json.dumps({"env": env_manifest(), "base_seed": BASE_SEED,
                                        "sa_params": asdict(SA_PARAMS), "sa_seed_offset": SA_SEED_OFFSET,
                                        "chi": a.chi, "started": time.strftime("%Y-%m-%dT%H:%M:%S")},
                                       indent=2), encoding="utf-8")
    recorded = json.loads(manifest.read_text(encoding="utf-8")).get("base_seed")
    if recorded != BASE_SEED:
        raise SystemExit(f"{out} was produced with base seed {recorded}, not {BASE_SEED}: use a different --out")
    if a.command in ("check", "all"):
        cmd_check(out)
    if a.command in ("mps-d3", "all"):
        cmd_mps_d3(out, a.d3_shots, a.chi)
    which = {"ours": {"ours"}, "baselines": {"baselines"}, "decode": {"ours", "baselines"},
             "all": {"ours", "baselines"}}.get(a.command)
    if which:
        tasks = build_tasks(str(out), a.chi, a.chunk, which)[a.shard::a.nshards]  # round-robin over cost-sorted list
        run_pool(tasks, a.workers)
    if a.command in ("analyze", "all"):
        cmd_analyze(out, a.chi)


if __name__ == "__main__":
    main()
