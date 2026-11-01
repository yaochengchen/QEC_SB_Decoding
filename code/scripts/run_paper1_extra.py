"""Paper 1 extra experiments on the existing grid shots (paired with results/paper1_grid).

  A) bposd   correlated BP+OSD-CS on ALL grid cells: one decoder on the stacked matrix
             [[Hz, 0, Hz], [0, Hx, Hx]] over (X, Z, Y) error columns, each with prior p/3,
             so Y correlations are modelled (the grid's BP_OSD decodes X and Z separately).
  B) budget  search-budget sweep at d = 9, 11, 13 and p = .08 .. .16, first BUDGET_SHOTS shots
             of each grid cell:
               CP-SB  (replicas, steps) in SB_CONFIGS   -> MAP, local SOE r1, local SOE r2
               CP-SA  (chains/mode, sweeps) in SA_CONFIGS -> MAP, local SOE r1, local SOE r2
             The base configs (SB 8x100, SA 8x100) must reproduce the grid predictions shot by
             shot (checked in `analyze`).
     v4:     the sweep is extended to p = .18 for the five configs in EXT_CONFIGS (threshold
             recovery with a scaled budget), and `budget_paired_matched.csv` compares SB and SA
             at the SAME SOE radius.
  C) scaled  (v4) budget-scaled rows for d >= 11: CP-SB 8x400 (pre-declared rule: steps = 400 for
             d >= 11, replicas unchanged) on the grid cells themselves -> MAP, SOE r1, SOE r2,
             paired with the frozen main method, CP-SA and the reference. All grid shots for
             p >= FULL_P_MIN; the low-p extension cells (p <= .06) are capped at FULL_LOWP_CAP
             shots (first shots of the cell, so the estimate is unbiased, only the CI is wider).
             The sweep's 8x400 files (first 10k shots, p=.08-.16) are reused as-is.

Subcommands
  run       all tasks (resumable; --shard/--nshards for Slurm arrays)
  hard64    MPS chi=64 re-check of the grid's hard shots (d>=11)
  analyze   join with the grid, LER / gap closure / paired tests / consistency checks / figures
  timing    single-process latency of every budget config + MPS chi=16/32 + PyMatching
            (run on an idle, exclusive node; one p per job, see slurm/submit_followup.sh)

  python scripts/run_paper1_extra.py run --grid results/paper1_grid --out results/paper1_extra --workers 512
  python scripts/run_paper1_extra.py analyze --grid results/paper1_grid --out results/paper1_extra
  python scripts/run_paper1_extra.py timing --out results/paper1_extra_timing_p0.10 --timing-p 0.10
"""

from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paper1_grid as grid  # noqa: E402
from run_paper1_v2 import (SA_SEED_OFFSET, import_bposd, env_manifest, write_atomic, gap_closure,  # noqa: E402
                           exact_paired_statistics, wilson_interval, adjust_p_values, _syn_key,
                           build_rotated_surface_code, sample_depolarizing, sample_sha256)
from qec_sb.adaptive_soe import fixed_soe_prediction  # noqa: E402
from qec_sb.cp_sa import CPSAParams, run_constraint_preserving_sa  # noqa: E402
from qec_sb.cp_sb import build_stabilizer_hamiltonian  # noqa: E402
from qec_sb.logical_classes import logical_class_seed  # noqa: E402
from qec_sb.p5.config import FROZEN_P42A  # noqa: E402
from qec_sb.scalable import (CSSMWPMDecoder, ScalableCPConfig, collect_cp_class_pools_no_oracle,  # noqa: E402
                             correction_class, cp_map_prediction)

BUDGET_D = (9, 11, 13)
BUDGET_P = (.08, .10, .12, .14, .16)
BUDGET_SHOTS = 10_000
SUB = {9: 250, 11: 125, 13: 50}       # shots per budget task (longest task ~20 min)
SB_CONFIGS = [(8, 100), (16, 100), (32, 100), (64, 100), (8, 200), (8, 400), (32, 400)]
SA_CONFIGS = [(8, 25), (8, 50), (8, 100), (8, 200), (16, 100), (32, 200)]
TIMING_SEED = 2_939_000_000

# ---- v4 additions -------------------------------------------------------------------------
BUDGET_P_EXT = (.18,)                                               # sweep extension (threshold recovery)
EXT_CONFIGS = ("sbR8S100", "sbR8S400", "sbR32S400", "saC8W100", "saC32W200")
FULL_D = (11, 13)                                                   # budget-scaled grid rows
FULL_CONFIGS = ("sbR8S400",)                                        # pre-declared: steps 400 for d>=11
FULL_P_MIN = .08                                                    # all grid shots for p >= this
FULL_LOWP_CAP = int(os.environ.get("QECSB_FULL_LOWP_CAP", "100000"))  # low-p extension cells: first N shots
FULL_SUB = {11: 500, 13: 200}                                       # shots per task outside the sweep range
TIMING_EXTRA = ("mps16", "mps32", "pym")                            # grid MethodRunner methods timed alongside


def sb_name(r, s):
    return f"sbR{r}S{s}"


def sa_name(c, w):
    return f"saC{c}W{w}"


ALL_BUDGET = [sb_name(*c) for c in SB_CONFIGS] + [sa_name(*c) for c in SA_CONFIGS]
# solver-seed replicates of the two base configs (same shots, independent solver randomness)
REPS = ["sbR8S100_k1", "sbR8S100_k2", "saC8W100_k1", "saC8W100_k2"]
REP_P = (.10, .14)
REP_SEED_STRIDE = 13_000_019


def _in(p, ps):
    return any(np.isclose(p, q) for q in ps)


# ------------------------------------------------------------------------- decoders
class CorrelatedBPOSD:
    """BP+OSD-CS on [[Hz,0,Hz],[0,Hx,Hx]] with columns (X, Z, Y), prior p/3 each."""

    def __init__(self, code, p, osd_order=7):
        import scipy.sparse as sp
        Hx, Hz = code.Hx.astype(np.uint8), code.Hz.astype(np.uint8)
        zx, zz = np.zeros_like(Hz), np.zeros_like(Hx)
        H = np.block([[Hz, zx, Hz], [zz, Hx, Hx]]).astype(np.uint8)
        n = code.n
        self.code, self.n = code, n
        self.dec = import_bposd()(sp.csr_matrix(H), error_channel=[float(p) / 3.0] * (3 * n),
                                  max_iter=3 * n, bp_method="minimum_sum", osd_method="osd_cs", osd_order=osd_order)

    def decode(self, sx, sz):
        e = np.asarray(self.dec.decode(np.concatenate([sx, sz]).astype(np.uint8)), np.uint8)
        x, z, y = e[: self.n], e[self.n: 2 * self.n], e[2 * self.n:]
        ex, ez = x ^ y, z ^ y
        if not (np.array_equal((ex @ self.code.Hz.T) & 1, sx) and np.array_equal((ez @ self.code.Hx.T) & 1, sz)):
            raise RuntimeError("correlated BP+OSD returned a syndrome-incompatible correction")
        return ex, ez


def sb_cfg(replicas, steps):
    f = FROZEN_P42A.cp
    return ScalableCPConfig(replicas=replicas, steps=steps, dt=f.dt, coupling_scale=f.coupling_scale,
                            warm_strength=f.warm_strength, random_restarts=f.random_restarts,
                            warm_restarts=f.warm_restarts, coupling_multipliers=tuple(f.coupling_multipliers))


def _infer(pools, hams, p):
    t0 = time.perf_counter(); m, _ = cp_map_prediction(pools, hams)
    t1 = time.perf_counter(); r1, *_ = fixed_soe_prediction(pools, hams, p, 1, 1)
    t2 = time.perf_counter(); r2, *_ = fixed_soe_prediction(pools, hams, p, 1, 2)
    t3 = time.perf_counter()
    return int(m), int(r1), int(r2), t1 - t0, t2 - t1, t3 - t2


class BudgetRunner:
    def __init__(self, code, p, method):
        self.code, self.p, self.method = code, float(p), method
        self.mwpm = CSSMWPMDecoder(code)
        base, _, rep = method.partition("_k")
        self.rep = int(rep) if rep else 0
        method = base
        if method == "bposd3":
            self.dec = CorrelatedBPOSD(code, p)
        elif method.startswith("sbR"):
            r, s = method[3:].split("S")
            self.cfg = sb_cfg(int(r), int(s))
        elif method.startswith("saC"):
            c, w = method[3:].split("W")
            self.sa = CPSAParams(chains=int(c), sweeps=int(w), temperature_start=2.0, temperature_end=0.05,
                                 restart_modes=("warm", "random"))
        else:
            raise ValueError(method)

    def decode(self, sx, sz, key, solver_seed):
        m, code, p = self.method, self.code, self.p
        t0 = time.perf_counter()
        if m == "bposd3":
            ex, ez = self.dec.decode(sx, sz)
            return {"pred_BP_OSD_corr": int(correction_class(code, ex, ez)), "time_bposd3_s": time.perf_counter() - t0}
        ex0, ez0 = self.mwpm.decode(sx, sz)
        t_seed = time.perf_counter() - t0
        t0 = time.perf_counter()
        if m.startswith("sbR"):   # identical seed arithmetic to the frozen pipeline
            pools, hams, _ = collect_cp_class_pools_no_oracle(code, ex0, ez0, solver_seed + key * 1009
                                                              + self.rep * REP_SEED_STRIDE, self.cfg)
        else:                     # identical seed arithmetic to the grid's CP-SA
            seed = solver_seed + SA_SEED_OFFSET + key * 1009 + self.rep * REP_SEED_STRIDE
            pools, hams = {}, {}
            for c in range(4):
                cex, cez = logical_class_seed(code, ex0, ez0, c)
                ham = build_stabilizer_hamiltonian(code, cex, cez)
                res = run_constraint_preserving_sa(ham, self.sa, int(seed + c * 100003))
                u = np.unique(np.vstack([np.zeros((1, ham.n_spins), np.uint8), res.terminal_bits]), axis=0)
                pools[c], hams[c] = u[np.argsort(ham.energy_bits(u), kind="stable")], ham
        t_search = time.perf_counter() - t0
        cmap, r1, r2, tm, t1, t2 = _infer(pools, hams, p)
        return {f"pred_{m}_MAP": cmap, f"pred_{m}_r1": r1, f"pred_{m}_r2": r2,
                f"t_{m}_seed_s": t_seed, f"t_{m}_search_s": t_search, f"t_{m}_r1_s": t1, f"t_{m}_r2_s": t2}


# ---------------------------------------------------------------------------- tasks
def grid_chunks():
    grid.SCALE = 1.0
    return {(t["d"], t["p"], t["chunk"]): t for t in grid.plan() if t["method"] == "mps16"}


def plan_tasks():
    tasks, seen = [], set()

    def add(d, p, k, lo, hi, m):
        sig = (d, round(float(p), 3), k, lo, hi, m)
        if sig not in seen and hi > lo:
            seen.add(sig)
            tasks.append({"d": d, "p": p, "chunk": k, "lo": lo, "hi": hi, "method": m})

    for (d, p, k), t in grid_chunks().items():
        n = t["hi"] - t["lo"]
        add(d, p, k, 0, n, "bposd3")
        # B) sweep (first BUDGET_SHOTS shots) -- original cells with every config, extension cells with EXT_CONFIGS
        sweep_cfgs = None
        if d in BUDGET_D and _in(p, BUDGET_P):
            sweep_cfgs = ALL_BUDGET + (REPS if _in(p, REP_P) else [])
        elif d in BUDGET_D and _in(p, BUDGET_P_EXT):
            sweep_cfgs = list(EXT_CONFIGS)
        covered = 0
        if sweep_cfgs is not None and t["lo"] < BUDGET_SHOTS:
            covered = min(n, BUDGET_SHOTS - t["lo"])
            for m in sweep_cfgs:
                for lo in range(0, covered, SUB[d]):
                    add(d, p, k, lo, min(lo + SUB[d], covered), m)
        # C) budget-scaled rows on the grid cells themselves (d>=11)
        if d in FULL_D:
            cap = n if p >= FULL_P_MIN - 1e-9 else max(0, min(n, FULL_LOWP_CAP - t["lo"]))
            for m in FULL_CONFIGS:
                start = covered if m in (sweep_cfgs or []) else 0
                for lo in range(start, cap, FULL_SUB[d]):
                    add(d, p, k, lo, min(lo + FULL_SUB[d], cap), m)
    return tasks


def task_path(out, t):
    return (Path(out) / "chunks" / grid.cell_id(t["d"], t["p"])
            / f"c{t['chunk']:05d}_s{t['lo']:05d}__{t['method']}.csv")


def cost(t):
    m = t["method"].partition("_k")[0]
    if m == "bposd3":
        w = .01
    elif m.startswith("sbR"):
        r, s = map(int, m[3:].split("S")); w = r * s / 800 + 4
    else:
        c, w_ = map(int, m[3:].split("W")); w = c * w_ / 200 + 4
    return -(t["d"] ** 3) * (t["hi"] - t["lo"]) * w


def run_one(out, t):
    path = task_path(out, t)
    if path.exists():
        return f"{path.name} cached"
    g = grid_chunks()[(t["d"], t["p"], t["chunk"])]
    code = build_rotated_surface_code(t["d"])
    s = sample_depolarizing(code, t["p"], g["hi"] - g["lo"], g["sample_seed"])
    runner = BudgetRunner(code, t["p"], t["method"])
    cache, rows, t0 = {}, [], time.time()
    for i in range(t["lo"], t["hi"]):
        key = _syn_key(s["sx"][i], s["sz"][i])
        if key not in cache:
            cache[key] = runner.decode(s["sx"][i], s["sz"][i], key, g["solver_seed"])
        rows.append({"d": t["d"], "p": t["p"], "chunk": t["chunk"], "shot": i, "syndrome_key": str(key),
                     "truth_class": int(s["logical"][i, 0]) * 2 + int(s["logical"][i, 1]),
                     "sample_sha256": sample_sha256(s), **cache[key]})
    path.parent.mkdir(parents=True, exist_ok=True)
    write_atomic(pd.DataFrame(rows), path)
    return f"{grid.cell_id(t['d'], t['p'])}/{path.name} {time.time() - t0:.0f}s"


def cmd_run(out, workers, shard, nshards):
    tasks = sorted(plan_tasks(), key=cost)[shard::nshards]
    tasks = [t for t in tasks if not task_path(out, t).exists()]
    grid.run_parallel(run_one, [(out, t) for t in tasks], workers, f"extra shard {shard}/{nshards}")


def cmd_estimate(out):
    """Print the remaining work (tasks not yet on disk) by group, in shots and rough core-hours."""
    per_shot = {("sbR8S400", 11): 2.1, ("sbR8S400", 13): 5.3, ("sbR32S400", 11): 4.7, ("sbR32S400", 13): 10.1,
                ("saC8W100", 11): 3.0, ("saC8W100", 13): 6.3, ("saC32W200", 11): 16.4, ("saC32W200", 13): 25.0,
                ("sbR8S100", 11): 1.3, ("sbR8S100", 13): 4.0}   # idle end-to-end (seed+search+r2) seconds, timing_budget.csv
    rows = {}
    for t in plan_tasks():
        if task_path(out, t).exists():
            continue
        if _in(t["p"], BUDGET_P_EXT) and t["lo"] < BUDGET_SHOTS and t["d"] in BUDGET_D:
            grp = "B) sweep extension p=.18"
        elif t["d"] in FULL_D and t["method"] in FULL_CONFIGS:
            grp = "C) budget-scaled rows"
        else:
            grp = "A/B) original (should be cached)"
        k = (grp, t["d"], t["method"])
        r = rows.setdefault(k, {"tasks": 0, "shots": 0, "core_h": 0.0})
        r["tasks"] += 1; r["shots"] += t["hi"] - t["lo"]
        r["core_h"] += (t["hi"] - t["lo"]) * per_shot.get((t["method"], t["d"]), 2.0) / 3600
    if not rows:
        print("nothing to do: every planned task is on disk"); return
    df = pd.DataFrame([{"group": g, "d": d, "method": m, **v} for (g, d, m), v in sorted(rows.items())])
    print(df.to_string(index=False))
    print(f"\nremaining ~{df.core_h.sum():.0f} core-hours (rough; idle-node timings, d<=9 counted at 2 s/shot)")


# ------------------------------------------------- MPS chi=64 on the grid's hard shots
def hard64_jobs(grid_out, out):
    jobs = []
    for f in sorted((Path(grid_out) / "hard").glob("d*_p*/c*__mps32.csv")):
        dst = Path(out) / "hard64" / f.parent.name / f.name.replace("mps32", "mps64")
        if not dst.exists():
            jobs.append((str(f), str(dst)))
    return jobs


def hard64_task(src, dst):
    from run_paper1_v2 import MethodRunner
    h = pd.read_csv(src)
    cell = Path(src).parent.name                      # d13_p0.120
    d, p = int(cell[1:cell.index("_")]), float(cell.split("_p")[1])
    k = int(h.chunk.iloc[0])
    g = next(t for (dd, pp, kk), t in grid_chunks().items() if dd == d and np.isclose(pp, p) and kk == k)
    code = build_rotated_surface_code(d)
    s = sample_depolarizing(code, g["p"], g["hi"] - g["lo"], g["sample_seed"])
    r = MethodRunner(code, g["p"], "mps64")
    rows = [{"chunk": k, "shot": int(i), **r.decode(s["sx"][i], s["sz"][i], _syn_key(s["sx"][i], s["sz"][i]), 0)}
            for i in h.shot]
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    write_atomic(pd.DataFrame(rows), dst)
    return f"{cell}/{Path(dst).name} {len(rows)} shots"


def cmd_hard64(grid_out, out, workers, shard, nshards):
    grid.run_parallel(hard64_task, hard64_jobs(grid_out, out)[shard::nshards], workers, f"hard64 shard {shard}/{nshards}")


# --------------------------------------------------------------------------- analyze
def load_extra(out):
    by_method = {}
    for f in (Path(out) / "chunks").glob("d*_p*/c*__*.csv"):
        by_method.setdefault(f.name.split("__")[1][:-4], []).append(pd.read_csv(f, dtype={"syndrome_key": str}))
    keys = ["d", "p", "chunk", "shot"]
    merged = None
    for m, dfs in sorted(by_method.items()):
        df = pd.concat(dfs, ignore_index=True)
        if df.duplicated(keys).any():
            raise ValueError(f"{m}: duplicated shots")
        cols = keys + ["truth_class", "sample_sha256"] + [c for c in df.columns if c.startswith(("pred_", "t_", "time_"))]
        df = df[cols].rename(columns={"truth_class": f"truth_{m}", "sample_sha256": f"sha_{m}"})
        merged = df if merged is None else merged.merge(df, on=keys, how="outer")
    return merged


def _has(g, col):
    return col in g and g[col].notna().any()


def _ler_row(g, col, ref_col="correct_REF"):
    c = g[col].dropna()
    n = len(c); k = int(n - c.sum()); lo, hi = wilson_interval(k, n)
    return {"shots": n, "failures": k, "LER": k / n, "LER_ci95_lo": lo, "LER_ci95_hi": hi,
            "LER_over_REF": (k / n) / max(1e-12, 1 - g.loc[c.index, ref_col].mean())}


def _pair(g, a, b):
    h = g[g[f"correct_{a}"].notna() & g[f"correct_{b}"].notna()]
    st = exact_paired_statistics(h[f"correct_{a}"].to_numpy(np.uint8), h[f"correct_{b}"].to_numpy(np.uint8))
    return {"shots": len(h), "LER_a_minus_LER_b": -st["paired_difference"], "ci95_lo": -st["paired_ci95_hi"],
            "ci95_hi": -st["paired_ci95_lo"], "a_right_b_wrong": st["discordant_wins"],
            "a_wrong_b_right": st["discordant_losses"], "exact_p": st["exact_p_value"]}


def _holm(df, by):
    if len(df):
        for key, idx in df.groupby(by).groups.items():
            df.loc[idx, "holm_p"] = adjust_p_values(df.loc[idx, "exact_p"].tolist(), "holm")
    return df


def load_timing(out):
    """Concatenate every timing_budget.csv next to the extra results (one directory per timing p)."""
    out = Path(out)
    cands = sorted(out.parent.glob("paper1_extra_timing*/timing_budget.csv")) + [out / "timing_budget.csv"]
    frames = [pd.read_csv(c) for c in cands if c.exists()]
    if not frames:
        return None
    t = pd.concat(frames, ignore_index=True)
    return t.drop_duplicates(["d", "p", "config", "stage"], keep="last")


def cmd_analyze(grid_out, out):
    out = Path(out)
    G = grid.load_grid(grid_out)
    G["p"] = G.p.round(3)
    G["correct_REF"] = np.nan
    for (d, p), idx in G.groupby(["d", "p"]).groups.items():
        ref = grid.reference_for(d, G.loc[idx])
        G.loc[idx, "correct_REF"] = G.loc[idx, f"correct_{ref}"].to_numpy()
    X = load_extra(out)
    X["p"] = X.p.round(3)
    keep = ["d", "p", "chunk", "shot", "truth_class", "sample_sha256", "error_weight", "correct_REF", "correct_MWPM",
            "pred_CP_SA_MAP", "pred_CP_SA_local_SOE", "pred_CP_SB_MAP", "pred_CP_local_SOE", "pred_CP_local_SOE_r2",
            "correct_CP_local_SOE_r2", "correct_CP_local_SOE", "correct_CP_SB_MAP", "correct_CP_SA_MAP",
            "correct_CP_SA_local_SOE", "correct_BP_OSD", "correct_PyMatching"]
    M = G[keep].merge(X, on=["d", "p", "chunk", "shot"], how="inner")
    for c in [c for c in M.columns if c.startswith("truth_") and c != "truth_class"]:
        bad = M[c].notna() & (M[c] != M.truth_class)
        if bad.any():
            raise ValueError(f"{c}: truth differs from the grid (wrong shots)")
        sha = "sha_" + c[6:]
        if (M[sha].notna() & (M[sha] != M.sample_sha256)).any():
            raise ValueError(f"{sha}: sample hash differs from the grid")
    for c in [c for c in M.columns if c.startswith("pred_") and c not in keep]:
        M["correct_" + c[5:]] = np.where(M[c].isna(), np.nan, (M[c] == M.truth_class).astype(float))
    cell_shots = G.groupby(["d", "p"]).size()

    # ---- A) correlated BP+OSD on the full grid
    A_rows, A_pairs = [], []
    for (d, p), g in M[M.pred_BP_OSD_corr.notna()].groupby(["d", "p"]):
        n = len(g)
        for m in ("BP_OSD", "BP_OSD_corr", "PyMatching", "CP_local_SOE_r2", "REF"):
            k = int(n - g[f"correct_{m}"].sum()); lo, hi = wilson_interval(k, n)
            row = {"d": d, "p": p, "method": m, "shots": n, "LER": k / n, "LER_ci95_lo": lo, "LER_ci95_hi": hi}
            if m not in ("REF",):
                row["gap_closed"] = gap_closure(g.correct_MWPM.to_numpy(), g[f"correct_{m}"].to_numpy(),
                                                g.correct_REF.to_numpy(), reps=300)[0]
            A_rows.append(row)
        for a, b in (("CP_local_SOE_r2", "BP_OSD_corr"), ("BP_OSD_corr", "BP_OSD")):
            A_pairs.append({"d": d, "p": p, "a": a, "b": b, **_pair(g, a, b)})
    A_pairs = _holm(pd.DataFrame(A_pairs), ["a", "b"])

    # ---- B) budget sweep (first BUDGET_SHOTS shots; p=.18 cells only carry EXT_CONFIGS)
    B = M[M.d.isin(BUDGET_D) & M[f"pred_{ALL_BUDGET[0]}_r2"].notna()].copy()
    checks = []
    for (d, p), g in B.groupby(["d", "p"]):
        checks.append({"d": d, "p": p, "shots": len(g),
                       "sbR8S100_matches_grid": bool((g.pred_sbR8S100_MAP == g.pred_CP_SB_MAP).all()
                                                     and (g.pred_sbR8S100_r1 == g.pred_CP_local_SOE).all()
                                                     and (g.pred_sbR8S100_r2 == g.pred_CP_local_SOE_r2).all()),
                       "saC8W100_matches_grid": bool((g.pred_saC8W100_MAP == g.pred_CP_SA_MAP).all()
                                                     and (g.pred_saC8W100_r1 == g.pred_CP_SA_local_SOE).all())})
    B_rows = []
    for (d, p), g in B.groupby(["d", "p"]):
        n = len(g)
        for m in ALL_BUDGET:
            if not _has(g, f"correct_{m}_r2"):
                continue
            for inf in ("MAP", "r1", "r2"):
                c = g[f"correct_{m}_{inf}"]
                k = int(n - c.sum()); lo, hi = wilson_interval(k, n)
                B_rows.append({"d": d, "p": p, "solver": "SB" if m.startswith("sb") else "SA", "config": m,
                               "inference": inf, "shots": n, "LER": k / n, "LER_ci95_lo": lo, "LER_ci95_hi": hi,
                               "gap_closed": gap_closure(g.correct_MWPM.to_numpy(), c.to_numpy(),
                                                         g.correct_REF.to_numpy(), reps=300)[0],
                               "LER_over_REF": (k / n) / max(1e-12, 1 - g.correct_REF.mean()),
                               "contended_mean_search_s": g[f"t_{m}_search_s"].mean()})
    B_rows = pd.DataFrame(B_rows)
    pairs = []
    for (d, p), g in B.groupby(["d", "p"]):
        for m in ALL_BUDGET:
            if _has(g, f"correct_{m}_r2"):
                pairs.append({"d": d, "p": p, "a": f"{m}_r2", "b": "saC8W100_r1 (grid CP-SA)",
                              **_pair(g, f"{m}_r2", "CP_SA_local_SOE")})
    pairs = _holm(pd.DataFrame(pairs), "a")
    # (v4) SB vs SA at the SAME SOE radius -- the comparison the paper's SB-vs-SA sentence rests on
    MATCHED = [("sbR8S400", "saC8W100"), ("sbR8S400", "saC8W50"), ("sbR32S400", "saC32W200"), ("sbR32S400", "saC8W200"),
               ("sbR8S100", "saC8W25"), ("sbR8S400", "sbR32S100"), ("sbR8S400", "sbR64S100")]
    mpairs = []
    for (d, p), g in B.groupby(["d", "p"]):
        for a, b in MATCHED:
            for inf in ("MAP", "r1", "r2"):
                if _has(g, f"correct_{a}_{inf}") and _has(g, f"correct_{b}_{inf}"):
                    mpairs.append({"d": d, "p": p, "a": f"{a}_{inf}", "b": f"{b}_{inf}", "family": f"{a} vs {b}",
                                   **_pair(g, f"{a}_{inf}", f"{b}_{inf}")})
    mpairs = _holm(pd.DataFrame(mpairs), ["family", "a"])

    # ---- C) solver-seed replicates of the base configs
    R_rows = []
    for (d, p), g in M[M.d.isin(BUDGET_D)].groupby(["d", "p"]):
        for base in ("sbR8S100", "saC8W100"):
            cols = [f"pred_{base}_r1"] + [f"pred_{base}_k{j}_r1" for j in (1, 2)]
            if not all(c in g for c in cols) or g[cols[1]].isna().all():
                continue
            h = g[g[cols[1]].notna() & g[cols[2]].notna() & g[cols[0]].notna()]
            lers = [1 - (h[c] == h.truth_class).mean() for c in cols]
            R_rows.append({"d": d, "p": p, "config": base + " + SOE r1", "shots": len(h),
                           "LER_seed0": lers[0], "LER_seed1": lers[1], "LER_seed2": lers[2],
                           "LER_spread": max(lers) - min(lers),
                           "all_three_agree": float(((h[cols[0]] == h[cols[1]]) & (h[cols[0]] == h[cols[2]])).mean())})
    # ---- D) MPS chi=64 vs chi=32 on the grid's hard shots (d>=11)
    H_rows = []
    for f in sorted((out / "hard64").glob("d*_p*/c*__mps64.csv")):
        h64 = pd.read_csv(f)
        h32 = pd.read_csv(Path(grid_out) / "hard" / f.parent.name / f.name.replace("mps64", "mps32"))
        j = h32.merge(h64, on=["chunk", "shot"])
        H_rows.append({"cell": f.parent.name, "shots": len(j),
                       "chi32_vs_chi64_disagreements": int((j.pred_MPS_chi32 != j.pred_MPS_chi64).sum())})
    H = pd.DataFrame(H_rows)
    if len(H):
        H = H.groupby("cell", as_index=False).sum()

    # ---- E) (v4) budget-scaled rows on the grid cells: d>=11, CP-SB 8x400 (MAP / r1 / r2)
    E_rows, E_pairs = [], []
    for cfg in FULL_CONFIGS:
        col = f"correct_{cfg}_r2"
        if col not in M:
            continue
        for (d, p), g in M[M.d.isin(FULL_D) & M[col].notna()].groupby(["d", "p"]):
            full = int(cell_shots.get((d, p), 0))
            for name, c in ((f"{cfg}_MAP", f"correct_{cfg}_MAP"), (f"{cfg}_r1", f"correct_{cfg}_r1"),
                            (f"{cfg}_r2", col), ("CP_SB_MAP (frozen)", "correct_CP_SB_MAP"),
                            ("CP_local_SOE (frozen r1)", "correct_CP_local_SOE"),
                            ("CP_local_SOE_r2 (frozen main)", "correct_CP_local_SOE_r2"),
                            ("CP_SA_MAP (grid)", "correct_CP_SA_MAP"), ("CP_SA_local_SOE (grid, r1)", "correct_CP_SA_local_SOE"),
                            ("REF", "correct_REF")):
                row = _ler_row(g, c)
                row["gap_closed"] = (np.nan if name == "REF" else
                                     gap_closure(g.correct_MWPM.to_numpy(), g[c].to_numpy(), g.correct_REF.to_numpy(), reps=300)[0])
                E_rows.append({"d": d, "p": p, "method": name, "cell_shots": full, "coverage": len(g) / max(1, full), **row})
            for a, b in ((f"{cfg}_r2", "CP_local_SOE_r2"), (f"{cfg}_r1", "CP_local_SOE"), (f"{cfg}_MAP", "CP_SB_MAP"),
                         (f"{cfg}_r1", "CP_SA_local_SOE"), (f"{cfg}_r2", "CP_SA_local_SOE"), (f"{cfg}_MAP", "CP_SA_MAP"),
                         (f"{cfg}_r2", f"{cfg}_r1"), (f"{cfg}_r2", "REF"), (f"{cfg}_r1", "REF")):
                E_pairs.append({"d": d, "p": p, "a": a, "b": b, "family": f"{a} vs {b}", **_pair(g, a, b)})
    E_rows, E_pairs = pd.DataFrame(E_rows), _holm(pd.DataFrame(E_pairs), "family")

    # ---- F) (v4) LER vs d with the budget-scaled main method (d<=9 frozen 8x100 r2, d>=11 8x400 r2)
    F_rows = []
    gl = pd.read_csv(Path(grid_out) / "grid_ler.csv") if (Path(grid_out) / "grid_ler.csv").exists() else None
    if gl is not None:
        gl["p"] = gl.p.round(3)
        for (d, p), g in gl.groupby(["d", "p"]):
            ref = g.reference.iloc[0]
            def pick(m):
                s = g[g.method == m]
                return (float(s.LER.iloc[0]), float(s.LER_ci95_lo.iloc[0]), float(s.LER_ci95_hi.iloc[0]), int(s.shots.iloc[0])) if len(s) else (np.nan,) * 4
            fr, sa, rf = pick("CP_local_SOE_r2"), pick("CP_SA_local_SOE"), pick(ref)
            sc = fr
            if d in FULL_D and len(E_rows):
                e = E_rows[(E_rows.d == d) & (E_rows.p == p) & (E_rows.method == f"{FULL_CONFIGS[0]}_r2")]
                if len(e):
                    sc = (float(e.LER.iloc[0]), float(e.LER_ci95_lo.iloc[0]), float(e.LER_ci95_hi.iloc[0]), int(e.shots.iloc[0]))
            F_rows.append({"d": d, "p": p, "budget_scaled_LER": sc[0], "budget_scaled_lo": sc[1], "budget_scaled_hi": sc[2],
                           "budget_scaled_shots": sc[3], "frozen_main_LER": fr[0], "CP_SA_local_SOE_LER": sa[0],
                           "REF_LER": rf[0], "budget_scaled_config": f"{FULL_CONFIGS[0]}_r2" if (d in FULL_D and sc is not fr) else "CP_local_SOE_r2"})
    F_rows = pd.DataFrame(F_rows)

    timing = load_timing(out)
    tables = {"bposd_corr_ler": pd.DataFrame(A_rows), "bposd_corr_paired": A_pairs,
              "budget_consistency": pd.DataFrame(checks), "budget_ler": B_rows, "budget_paired_vs_SA": pairs,
              "budget_paired_matched": mpairs, "budget_scaled_ler": E_rows, "budget_scaled_paired": E_pairs,
              "ler_vs_d_budget_scaled": F_rows,
              "seed_replicates": pd.DataFrame(R_rows), "mps_chi64_hard_check": H}
    for k, t in tables.items():
        t.to_csv(out / f"{k}.csv", index=False)
    plot_budget(B_rows, timing, out / "fig_budget_accuracy_vs_cost.pdf", pd.DataFrame(A_rows))
    if len(F_rows):
        plot_ler_vs_d(F_rows, out / "fig_ler_vs_d_budget_scaled.pdf")
    if gl is not None:
        plot_decomposition(gl, E_rows, out / "fig_decomposition.pdf")
    md = ["# Paper 1 extra experiments", "",
          "A) correlated BP+OSD on every grid cell; B) CP-SB / CP-SA search-budget sweep (d=9,11,13; p=.18 added in v4);",
          "C) budget-scaled rows (d>=11, CP-SB 8x400) on the grid cells; D) chi=64 check; seed replicates.",
          "REF = exact MLD (d=3), MPS chi=32 (d<=9), MPS chi=16 + chi=32 on hard shots (d>=11).",
          "Search times in `budget_ler` were measured under full load; the cost axis of the figure uses the idle-node "
          "`timing_budget.csv` files when present.", ""]
    for k, t in tables.items():
        md += [f"## {k}", "", t.to_markdown(index=False, floatfmt=".4g") if len(t) else "(empty)", ""]
    (out / "RESULTS_EXTRA.md").write_text("\n".join(md), encoding="utf-8")
    print(f"wrote {out / 'RESULTS_EXTRA.md'}")


# ---------------------------------------------------------------------------- figures
_COL = {"SB": "#2a78d6", "SA": "#eb6834", "BP": "#1baf7a", "MPS": "#0b0b0b", "MAP": "#6250d6", "r1": "#eb6834", "r2": "#2a78d6"}


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _timing_lookup(timing, d, p):
    """stage-> median seconds for (d, nearest timing p); None if no timing."""
    if timing is None:
        return None
    t = timing[timing.d == d]
    if not len(t):
        return None
    pt = t.p.unique()[np.argmin(np.abs(t.p.unique() - p))]
    t = t[np.isclose(t.p, pt)]
    return {(c, s): v for c, s, v in zip(t.config, t.stage, t.median_s)}


def plot_budget(B, timing, path, A=None):
    """LER/REF vs END-TO-END time (seed + search + SOE) for SB and SA with r1 and r2; MPS/BP+OSD points when timed."""
    plt = _plt()
    ps = sorted(B.p.unique())
    fig, axes = plt.subplots(len(BUDGET_D), len(ps), figsize=(3.0 * len(ps), 2.8 * len(BUDGET_D)), sharey="row", squeeze=False)
    for i, d in enumerate(BUDGET_D):
        for j, p in enumerate(ps):
            ax = axes[i, j]
            tl = _timing_lookup(timing, d, p)
            for solver in ("SB", "SA"):
                for inf, ls, mk in (("r1", "-", "o"), ("r2", "--", "s")):
                    s = B[(B.d == d) & np.isclose(B.p, p) & (B.solver == solver) & (B.inference == inf)].copy()
                    if not len(s):
                        continue
                    if tl is not None:
                        s["x"] = [tl.get((c, "seed"), np.nan) + tl.get((c, "search"), np.nan) + tl.get((c, f"soe_{inf}"), np.nan)
                                  for c in s.config]
                    else:
                        s["x"] = s.contended_mean_search_s
                    s = s.dropna(subset=["x"]).sort_values("x")
                    ax.plot(s.x, s.LER_over_REF, ls, marker=mk, ms=4, lw=1.8, color=_COL[solver], label=f"CP-{solver} + SOE {inf}")
            if tl is not None:
                pts = [("mps16", "MPS chi=16 (reference class)", "*", _COL["MPS"], 1.0)]
                if A is not None and len(A):
                    a = A[(A.d == d) & np.isclose(A.p, p)]
                    if len(a) and (a.method == "BP_OSD_corr").any() and (a.method == "REF").any():
                        pts.append(("bposd3", "BP+OSD corr.", "D", _COL["BP"],
                                    float(a[a.method == "BP_OSD_corr"].LER.iloc[0]) / max(1e-12, float(a[a.method == "REF"].LER.iloc[0]))))
                for cfg, lab, mk, col, y in pts:
                    if (cfg, "total") in tl and np.isfinite(y):
                        ax.plot([tl[(cfg, "total")]], [y], mk, ms=10 if mk == "*" else 6, color=col, ls="none", label=lab)
            ax.axhline(1, color="#9a9994", lw=1, ls=":")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=.25, which="both")
            ax.set_title(f"d={d}, p={p:.2f}", fontsize=9)
            if i == len(BUDGET_D) - 1:
                ax.set_xlabel("end-to-end time per syndrome (s)" + ("" if timing is not None else ", contended search only"))
        axes[i, 0].set_ylabel("LER / LER(ref)")
    axes[0, -1].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=160)


def plot_ler_vs_d(F, path):
    """LER vs d at fixed p: frozen main (8x100 everywhere) vs budget-scaled main (8x400 for d>=11) vs REF."""
    plt = _plt()
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    for ax, (col, title) in zip(axes, (("frozen_main_LER", "CP-SB 8x100 + SOE r2 (frozen)"),
                                       ("budget_scaled_LER", "budget-scaled: 8x400 + SOE r2 for d>=11"),
                                       ("REF_LER", "reference (MLD / MPS)"))):
        for p, s in F.groupby("p"):
            s = s.sort_values("d")
            ax.plot(s.d, s[col].clip(lower=1e-7), marker="o", ms=3, lw=1.5, label=f"p={p:.2f}")
        ax.set_yscale("log"); ax.set_xlabel("d"); ax.set_title(title, fontsize=9); ax.grid(alpha=.25, which="both")
    axes[0].set_ylabel("logical error rate"); axes[-1].legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=160)


def plot_decomposition(gl, E, path):
    """Where the gap to ML comes from: MAP (search) vs SOE r1 vs SOE r2, frozen budget (solid) and 8x400 (dashed, d>=11)."""
    plt = _plt()
    ds = sorted(gl.d.unique())
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6), sharey=True, squeeze=False)
    for ax, d in zip(axes.flat, ds):
        g = gl[gl.d == d]
        ref = g[g.method == g.reference].set_index("p").LER
        for m, lab, key in (("CP_SB_MAP", "CP-SB MAP (search only)", "MAP"), ("CP_local_SOE", "+ local SOE r1", "r1"),
                            ("CP_local_SOE_r2", "+ local SOE r2", "r2")):
            s = g[g.method == m].set_index("p").LER.sort_index()
            ax.plot(s.index, (s / ref.reindex(s.index)).to_numpy(), "-", marker="o", ms=3, lw=1.8, color=_COL[key], label=lab + " (8x100)")
        if d in FULL_D and len(E):
            e = E[E.d == d]
            for inf, lab in (("MAP", "8x400 MAP"), ("r1", "8x400 + r1"), ("r2", "8x400 + r2")):
                s = e[e.method == f"{FULL_CONFIGS[0]}_{inf}"].sort_values("p")
                if len(s):
                    ax.plot(s.p, s.LER_over_REF, "--", marker="s", ms=3, lw=1.5, color=_COL[inf], label=lab)
        ax.axhline(1, color="#9a9994", lw=1, ls=":"); ax.set_yscale("log"); ax.grid(alpha=.25, which="both")
        ax.set_title(f"d={d}", fontsize=9); ax.set_xlabel("p")
    for ax in axes[:, 0]:
        ax.set_ylabel("LER / LER(ref)")
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=160)


# ---------------------------------------------------------------------------- timing
def cmd_timing(out, shots, ps):
    """Single process, one p per call. Budget configs report stages; MPS/PyMatching report `total` only."""
    from run_paper1_v2 import MethodRunner
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for p in ps:
        for d in BUDGET_D:
            code = build_rotated_surface_code(d)
            s = sample_depolarizing(code, p, shots, TIMING_SEED + d)
            keys = []
            seen = set()
            for i in range(shots):
                key = _syn_key(s["sx"][i], s["sz"][i])
                if key in seen or key == 0:
                    continue
                seen.add(key); keys.append(i)
            for m in ["bposd3"] + ALL_BUDGET + list(TIMING_EXTRA):
                if m in TIMING_EXTRA:
                    r = MethodRunner(code, p, m)
                    rec = []
                    for i in keys:
                        t0 = time.perf_counter()
                        r.decode(s["sx"][i], s["sz"][i], _syn_key(s["sx"][i], s["sz"][i]), TIMING_SEED + 500_000 + d)
                        rec.append({"t_total": time.perf_counter() - t0})
                    df = pd.DataFrame(rec); stages = {"total": ["t_total"]}
                else:
                    r = BudgetRunner(code, p, m)
                    df = pd.DataFrame([r.decode(s["sx"][i], s["sz"][i], _syn_key(s["sx"][i], s["sz"][i]), TIMING_SEED + 500_000 + d)
                                       for i in keys])
                    stages = {"total": [c for c in df.columns if c.startswith(("t_", "time_"))]}
                    if m != "bposd3":
                        stages.update({"search": [f"t_{m}_search_s"], "seed": [f"t_{m}_seed_s"],
                                       "soe_r1": [f"t_{m}_r1_s"], "soe_r2": [f"t_{m}_r2_s"],
                                       "end_to_end_r1": [f"t_{m}_seed_s", f"t_{m}_search_s", f"t_{m}_r1_s"],
                                       "end_to_end_r2": [f"t_{m}_seed_s", f"t_{m}_search_s", f"t_{m}_r2_s"]})
                for st, cols in stages.items():
                    x = df[cols].sum(axis=1).to_numpy()
                    rows.append({"d": d, "p": p, "config": m, "stage": st, "syndromes": len(x),
                                 "median_s": float(np.median(x)), "p90_s": float(np.percentile(x, 90))})
                tot = next(r_ for r_ in reversed(rows) if r_["stage"] == "total")
                print(f"p={p} d={d} {m:10s} total median {tot['median_s']*1e3:9.1f} ms  ({len(keys)} syndromes)", flush=True)
            pd.DataFrame(rows).to_csv(out / "timing_budget.csv", index=False)   # checkpoint after every d
    t = pd.DataFrame(rows)
    t.to_csv(out / "timing_budget.csv", index=False)
    (out / "TIMING_BUDGET.md").write_text("```json\n" + json.dumps(env_manifest(), indent=1) + "\n```\n\n"
                                          + t.to_markdown(index=False, floatfmt=".4g") + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "hard64", "analyze", "timing", "estimate"])
    ap.add_argument("--grid", default="results/paper1_grid")
    ap.add_argument("--out", default="results/paper1_extra")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--timing-shots", type=int, default=100)
    ap.add_argument("--timing-p", type=float, nargs="+", default=[.10])
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.command == "timing":
        cmd_timing(out, a.timing_shots, a.timing_p)
        return
    if a.command == "estimate":
        cmd_estimate(str(out))
        return
    m = out / "run_manifest.json"
    if not m.exists():
        m.write_text(json.dumps({"env": env_manifest(), "grid": a.grid, "budget_d": BUDGET_D, "budget_p": BUDGET_P,
                                 "budget_shots": BUDGET_SHOTS, "sb_configs": SB_CONFIGS, "sa_configs": SA_CONFIGS},
                                indent=2), encoding="utf-8")
    m4 = out / "run_manifest_v4.json"
    if not m4.exists():
        m4.write_text(json.dumps({"env": env_manifest(), "budget_p_ext": BUDGET_P_EXT, "ext_configs": EXT_CONFIGS,
                                  "full_d": FULL_D, "full_configs": FULL_CONFIGS, "full_p_min": FULL_P_MIN,
                                  "full_lowp_cap": FULL_LOWP_CAP, "full_sub": FULL_SUB,
                                  "rule": "pre-declared before any d>=11 budget-scaled shot was decoded: CP-SB replicas 8, "
                                          "steps 400 for d >= 11 (chosen on the 10k-shot sweep, p=.08-.16); "
                                          "inference r2 (main, as fixed on the confirmation set) with r1 reported alongside."},
                                 indent=2), encoding="utf-8")
    if a.command == "run":
        cmd_run(str(out), a.workers, a.shard, a.nshards)
    elif a.command == "hard64":
        cmd_hard64(a.grid, str(out), a.workers, a.shard, a.nshards)
    else:
        cmd_analyze(a.grid, out)


if __name__ == "__main__":
    main()
