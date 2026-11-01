"""Paper-1 appendix: locality of the circuit-level DEM stabilizer group (CPU only, no decoding).

  python scripts/run_circuit_locality.py run --out results/circuit_locality --workers 512 --workers-big 96
  python scripts/run_circuit_locality.py summary --out results/circuit_locality

Instances (default): rotated surface code d = 3,5,7,9,11,13 with rounds = d (Stim generated rotated_memory_z,
after_clifford_depolarization = before_round_data_depolarization = before_measure_flip = after_reset_flip = p), and the
[[72,12,6]] BB code with the Bravyi et al. depth-8 cycle, 6 rounds, same p.  p = 0.003.
Per instance (resumable, each step skipped if its output exists; a failing instance is logged and skipped):
  build   -> H.npz, L.npz, mech.npz, circuit.stim, dem.dem, info.json        (Stim DEM, hyperedges kept)
  t1      -> G.npz, Lg.npz, generator_stats.json, T1_REPORT.md              (self-decoded generating set of ker([H;L]))
  inspect -> T1_global_elements.csv, T1_INSPECT.md                           (what the heavy / global elements are)
summary (rewritten after every instance) -> <out>/summary/{locality_summary.csv, LOCALITY_SUMMARY.md, locality_table.tex}
"Term order" of mechanism i = number of generators containing it = spin order of its energy term in generator coordinates
(code capacity in the paper: <= 4).
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from qec_sb.circuit_level import circuits, codes  # noqa: E402
from qec_sb.circuit_level.dem import DEMMatrices, dem_matrices  # noqa: E402
from qec_sb.circuit_level.kernel_basis import build_generator_set, verdict  # noqa: E402

DEFAULT_CODES = ["sc3", "sc5", "sc7", "sc9", "sc11", "sc13", "bb72B"]


def log(msg, out=None):
    line = f"{time.strftime('%F %T')}  {msg}"
    print(line, flush=True)
    if out is not None:
        with open(Path(out) / "locality.log", "a") as f:
            f.write(line + "\n")


def instance_name(code, p, rounds):
    return f"{code}_p{p:.4f}_r{rounds}"


def rounds_of(code, a):
    if code.startswith("sc"):
        return int(code[2:]) if a.rounds_sc <= 0 else a.rounds_sc
    return a.rounds_bb


def build_circuit(code, p, rounds):
    if code.startswith("sc"):
        return circuits.surface_memory(int(code[2:]), p, rounds), None
    if code.endswith("B"):
        base = code[:-1]
        c = codes.bb_code(base)
        return circuits.bb_memory_z_bravyi(base, c, p, rounds), c
    raise ValueError(f"unknown code {code} (use scD or bb72B)")


# ------------------------------------------------------------------ I/O
def save_dem(path: Path, dm: DEMMatrices, circuit_text: str):
    path.mkdir(parents=True, exist_ok=True)
    sp.save_npz(path / "H.npz", dm.H.tocsr())
    sp.save_npz(path / "L.npz", dm.L.tocsr())
    np.savez(path / "mech.npz", probs=dm.probs, weights=dm.weights, tick=dm.tick)
    (path / "circuit.stim").write_text(circuit_text)
    (path / "dem.dem").write_text(str(dm.dem))


def load_dem(path: Path) -> DEMMatrices:
    import stim
    z = np.load(path / "mech.npz")
    return DEMMatrices(sp.load_npz(path / "H.npz").tocsr(), sp.load_npz(path / "L.npz").tocsr(), z["probs"],
                       z["weights"], z["tick"], stim.DetectorErrorModel((path / "dem.dem").read_text()))


def load_generators(path: Path):
    return sp.load_npz(path / "G.npz").tocsr(), sp.load_npz(path / "Lg.npz").tocsr()


def write_json(path: Path, obj):
    def _d(o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=_d))
    tmp.replace(path)


# ------------------------------------------------------------------ steps
def step_build(inst: Path, code, p, rounds):
    t0 = time.time()
    circuit, c = build_circuit(code, p, rounds)
    dm = dem_matrices(circuit)
    save_dem(inst, dm, str(circuit))
    info = {"instance": inst.name, "code": code, "p": p, "rounds": rounds, "n_qubits_circuit": circuit.num_qubits,
            "n_det": dm.n_det, "n_obs": dm.n_obs, "n_mech": dm.n_mech,
            "det_per_mech_mean": float(np.asarray(dm.H.sum(axis=0)).mean()),
            "ticks_known_frac": float((dm.tick >= 0).mean()), "n_ticks": int(dm.tick.max() + 1),
            "noise": ("stim rotated_memory_z, after_clifford_depolarization = before_round_data_depolarization = "
                      "before_measure_flip_probability = after_reset_flip_probability = p") if code.startswith("sc") else
                     ("Bravyi et al. depth-8 cycle: DEPOLARIZE2 after CNOT, DEPOLARIZE1 on idle data, flip after reset "
                      "and before measurement, all = p"),
            "seconds": time.time() - t0}
    if c is not None:
        info.update({"n": c.n, "k": c.k, "schedule": "bravyi_depth8"})
    write_json(inst / "info.json", info)
    return info


def step_t1(inst: Path, workers: int, completion: str):
    dm = load_dem(inst)
    chunk = max(10, dm.n_mech // (4 * max(workers, 1)))           # keep all workers busy
    gs = build_generator_set(dm, workers=workers, osd_order=7, jitter_trials=1, chunk=chunk, completion=completion,
                             log=lambda m: print(m, flush=True))
    ok, msg = verdict(gs.stats)
    gs.stats.update({"verdict": msg, "workers": workers, "chunk": chunk})
    sp.save_npz(inst / "G.npz", gs.G.tocsr())
    sp.save_npz(inst / "Lg.npz", gs.Lg.tocsr())
    write_json(inst / "generator_stats.json", gs.stats)
    hist = np.bincount(gs.weights.astype(int))
    s = gs.stats
    lines = [f"# Locality of the DEM stabilizer generating set: {inst.name}", "", f"**Gate (Phase-0 definition): {msg}**", "",
             f"- mechanisms {s['n_mech']}, detectors {s['n_det']}, observables {s['n_obs']}; stabilizer dimension {s['dim_stabilizer']}",
             f"- generator weight: mean {s['gen_weight_mean']:.2f}, p50 {s['gen_weight_p50']:.0f}, p99 {s['gen_weight_p99']:.0f}, "
             f"max {s['gen_weight_max']}; self-decoded max {s['selfdecoded_weight_max']}",
             f"- null-space completions (global elements) {s['n_completion']} (max weight {s['completion_weight_max']}, mode {completion})",
             f"- generators per mechanism: mean {s['col_weight_mean']:.2f}, max {s['col_weight_max']}",
             f"- logical representative weights {s['logical_weights']}; wall time {s['seconds']:.0f}s with {workers} workers",
             "", "| weight | count |", "|---|---|"] + [f"| {k} | {int(v)} |" for k, v in enumerate(hist) if v]
    (inst / "T1_REPORT.md").write_text("\n".join(lines) + "\n")
    return s


def step_inspect(inst: Path, heavy: int = 12):
    dm = load_dem(inst)
    G, _ = load_generators(inst)
    st = json.loads((inst / "generator_stats.json").read_text())
    w = np.diff(G.indptr)
    m = G.shape[0]
    n_comp = int(st.get("n_completion", 0))
    comp_rows = set(range(m - n_comp, m))
    rows_idx = sorted(set(np.nonzero(w > heavy)[0].tolist()) | comp_rows)
    coords = dm.dem.get_detector_coordinates()
    ndim = max((len(v) for v in coords.values()), default=0)
    Hc = dm.H.tocsc()
    rows = []
    for g in rows_idx:
        mech = G.indices[G.indptr[g]:G.indptr[g + 1]]
        dets = Counter()
        for i in mech:
            dets.update(Hc.indices[Hc.indptr[i]:Hc.indptr[i + 1]].tolist())
        touched = sorted(dets)
        C = np.asarray([(list(coords.get(d, [])) + [np.nan] * ndim)[:ndim] for d in touched], float) if touched \
            else np.zeros((0, ndim))
        tk = dm.tick[mech]
        tk = tk[tk >= 0]
        r = {"row": int(g), "weight": int(w[g]), "is_completion": g in comp_rows, "n_det_touched": len(touched),
             "odd_touch": int(sum(c % 2 for c in dets.values())),
             "tick_min": int(tk.min()) if len(tk) else -1, "tick_max": int(tk.max()) if len(tk) else -1}
        for k in range(ndim):
            col = C[:, k][np.isfinite(C[:, k])] if len(C) else np.zeros(0)
            r[f"c{k}_n_distinct"] = int(len(np.unique(col)))
            r[f"c{k}_min"] = float(col.min()) if len(col) else np.nan
            r[f"c{k}_max"] = float(col.max()) if len(col) else np.nan
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(inst / "T1_global_elements.csv", index=False)
    conv = "BB [check, round, type 0=Z 1=X]" if inst.name.startswith("bb") else "surface code [x, y, t]"
    lines = [f"# Heavy and global generators: {inst.name}", "",
             f"- generators {m}; global (null-space completions) {n_comp}; heavier than {heavy}: {int((w > heavy).sum())}",
             f"- detector coordinates: {conv}; odd_touch must be 0 (kernel element)", ""]
    if len(df):
        assert int(df.odd_touch.sum()) == 0, "a listed generator is not a kernel element"
        lines += [df.to_string(index=False)]
    (inst / "T1_INSPECT.md").write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ summary
def _q(x, p):
    return float(np.percentile(x, p)) if len(x) else float("nan")


def summarize_one(inst: Path) -> dict:
    dm = load_dem(inst)
    G, _ = load_generators(inst)
    st = json.loads((inst / "generator_stats.json").read_text())
    info = json.loads((inst / "info.json").read_text())
    m = G.shape[0]
    nc = int(st.get("n_completion", 0))
    ns = m - nc
    w = np.diff(G.indptr)
    ws, wg = w[:ns], w[ns:]
    colw = np.diff(G.tocsc().indptr)
    colw_s = np.diff(G[:ns].tocsc().indptr)
    glob_mech = np.unique(G[ns:].indices) if nc else np.zeros(0, int)
    t_idx = 1 if inst.name.startswith("bb") else 2
    t_of_det = np.full(dm.n_det, np.nan)
    for d, c in dm.dem.get_detector_coordinates().items():
        if len(c) > t_idx:
            t_of_det[d] = c[t_idx]
    n_rounds_det = int(len(np.unique(t_of_det[np.isfinite(t_of_det)])))
    Hc = dm.H.tocsc()
    rspan = []
    for g in range(ns, m):
        mech = G.indices[G.indptr[g]:G.indptr[g + 1]]
        dets = np.unique(np.concatenate([Hc.indices[Hc.indptr[i]:Hc.indptr[i + 1]] for i in mech]))
        tt = t_of_det[dets]
        rspan.append(int(len(np.unique(tt[np.isfinite(tt)]))))
    return {
        "instance": inst.name, "code": info["code"], "p": info["p"], "rounds": info["rounds"], "n_det": dm.n_det,
        "n_obs": dm.n_obs, "n_mech": dm.n_mech, "dim_spins": m, "completion_mode": st.get("completion_mode"),
        "t1_seconds": st.get("seconds"), "t1_workers": st.get("workers"),
        "self_w_p50": _q(ws, 50), "self_w_p99": _q(ws, 99), "self_w_max": int(ws.max()),
        "self_frac_w_le4": float((ws <= 4).mean()), "self_frac_w_le8": float((ws <= 8).mean()),
        "n_global": nc, "global_frac_of_dim": nc / m,
        "global_w_min": int(wg.min()) if nc else 0, "global_w_median": _q(wg, 50) if nc else 0,
        "global_w_max": int(wg.max()) if nc else 0,
        "global_round_span_min": int(min(rspan)) if nc else 0, "global_round_span_median": _q(rspan, 50) if nc else 0,
        "n_detector_rounds": n_rounds_det, "mech_touched_by_global_frac": len(glob_mech) / dm.n_mech,
        "all_frac_w_le8": float((w <= 8).mean()), "all_w_max": int(w.max()),
        "order_mean": float(colw.mean()), "order_p50": _q(colw, 50), "order_p99": _q(colw, 99),
        "order_max": int(colw.max()), "order_frac_le4": float((colw <= 4).mean()),
        "order_self_p99": _q(colw_s, 99), "order_self_max": int(colw_s.max()),
        "logical_rep_weight_min": int(min(st.get("logical_weights") or [0])),
    }


def write_summary(out: Path, insts):
    rows = [summarize_one(i) for i in insts if (i / "G.npz").exists() and (i / "info.json").exists()]
    if not rows:
        return
    sd = out / "summary"
    sd.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(sd / "locality_summary.csv", index=False)
    cols = ["instance", "n_mech", "dim_spins", "self_w_p50", "self_w_p99", "self_w_max", "self_frac_w_le8", "n_global",
            "global_w_min", "global_w_max", "global_round_span_min", "n_detector_rounds", "mech_touched_by_global_frac",
            "order_mean", "order_p99", "order_max", "order_self_max", "order_frac_le4", "logical_rep_weight_min",
            "completion_mode", "t1_seconds", "t1_workers"]
    (sd / "LOCALITY_SUMMARY.md").write_text(
        "# Locality of the circuit-level DEM stabilizer group (paper-1 appendix)\n\n"
        "self = self-decoded generators; global = null-space completions; order = generators per mechanism (spin order of "
        "that energy term); global_round_span_min = fewest detector rounds touched by a global element.\n\n"
        + df[cols].to_markdown(index=False, floatfmt=".4g") + "\n")
    tex = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
           r"code & mech. & spins & weight $p_{50}/p_{99}$ & $\le 8$ & global (weight) & order mean/$p_{99}$/max\\", r"\midrule"]
    for r in rows:
        lab = f"surface $d={r['code'][2:]}$" if r["code"].startswith("sc") else r"BB $[[72,12,6]]$"
        glob = f"{r['n_global']} ({r['global_w_min']}--{r['global_w_max']})" if r["n_global"] else "0"
        mech = f"{r['n_mech']:,}".replace(",", "{,}")
        spins = f"{r['dim_spins']:,}".replace(",", "{,}")
        tex.append(f"{lab} & {mech} & {spins} & {r['self_w_p50']:.0f}/{r['self_w_p99']:.0f} & "
                   f"{100 * r['self_frac_w_le8']:.2f}\\% & {glob} & {r['order_mean']:.2f}/{r['order_p99']:.0f}/{r['order_max']}\\\\")
    tex += [r"\bottomrule", r"\end{tabular}"]
    (sd / "locality_table.tex").write_text("\n".join(tex) + "\n")


# ------------------------------------------------------------------ CLI
def cmd_run(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    insts = []
    for code in a.codes:
        rounds = rounds_of(code, a)
        inst = out / instance_name(code, a.p, rounds)
        insts.append(inst)
        big = code.startswith("bb") or (code.startswith("sc") and int(code[2:]) >= a.big_from)
        workers = a.workers_big if big else a.workers
        completion = a.completion_bb if code.startswith("bb") else a.completion_sc
        try:
            if not (inst / "info.json").exists():
                info = step_build(inst, code, a.p, rounds)
                log(f"{inst.name} built: {info['n_det']} detectors, {info['n_mech']} mechanisms, {info['seconds']:.0f}s", out)
            if not (inst / "G.npz").exists():
                t0 = time.time()
                s = step_t1(inst, workers, completion)
                log(f"{inst.name} t1: dim {s['dim_stabilizer']}, completions {s['n_completion']}, self-decoded max weight "
                    f"{s['selfdecoded_weight_max']}, {time.time() - t0:.0f}s ({workers} workers, {completion})", out)
            if not (inst / "T1_INSPECT.md").exists():
                step_inspect(inst)
            write_summary(out, insts)
            log(f"{inst.name} done", out)
        except Exception as exc:                      # includes BrokenProcessPool (worker killed, e.g. out of memory)
            log(f"{inst.name} FAILED: {exc!r} -- rerun the same command (resumes); if a worker was killed, lower --workers-big", out)
            traceback.print_exc()
    write_summary(out, insts)
    log(f"summary -> {out / 'summary'}", out)


def cmd_summary(a):
    out = Path(a.out)
    insts = sorted([p for p in out.iterdir() if p.is_dir() and (p / "G.npz").exists()],
                   key=lambda p: json.loads((p / "info.json").read_text())["n_mech"])
    write_summary(out, insts)
    print((out / "summary" / "LOCALITY_SUMMARY.md").read_text())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", default="results/circuit_locality")
    r.add_argument("--codes", nargs="+", default=DEFAULT_CODES)
    r.add_argument("--p", type=float, default=0.003)
    r.add_argument("--rounds-sc", type=int, default=0, help="0 = rounds equal to d")
    r.add_argument("--rounds-bb", type=int, default=6)
    r.add_argument("--workers", type=int, default=64)
    r.add_argument("--workers-big", type=int, default=96, help="workers for sc d >= --big-from and for bb (memory)")
    r.add_argument("--big-from", type=int, default=11)
    r.add_argument("--completion-sc", choices=["full", "fast"], default="full")
    r.add_argument("--completion-bb", choices=["full", "fast"], default="fast",
                   help="fast: global elements are lightened greedily only (their weights are upper bounds); full ~2 h on bb72")
    s = sub.add_parser("summary")
    s.add_argument("--out", default="results/circuit_locality")
    a = ap.parse_args()
    {"run": cmd_run, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()
