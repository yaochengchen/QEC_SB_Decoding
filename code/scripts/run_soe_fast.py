"""Sparse search kernels + O(m) local SOE: identity check, timing and implementation-independent work counts.

What changed (none of it changes a single decision, so no decoding is rerun):
  * qec_sb/cp_sb.py, qec_sb/cp_sa.py: the single-generator flip energy dE_j is computed from the <= 4 qubits
    generator j touches ("sparse", now the default) instead of the whole 2n-bit error ("dense", QECSB_DELTA=dense).
  * qec_sb/fast_soe.py: SOE r = 1, 2 from the weight histogram in O(m) instead of enumerating the ball.

Subcommands (all output under the new folder results_soe_fast/; results/ is not touched)
  check    on the same shots: dense vs sparse search pools must be identical (CP-SB 8x100 all d, 8x400 d>=11,
           CP-SA 8x100 all d), and old enumeration vs fast SOE (r1, r2) must give the same scores and decisions.
           d = 3..13, p = .06/.10/.14/.18. Parallel, resumable.
           -> results_soe_fast/check_v2/SOE_FAST_CHECK.md   (first line ALL_OK / FAIL / INCOMPLETE)
  timing   single process on an idle exclusive node; same syndromes, solver seeds and pools as
           `run_paper1_extra.py timing` (TIMING_SEED), d = 9/11/13; every budget config (7 SB + 6 SA) with the
           sparse kernels and fast SOE; on sbR8S100 and saC8W100 also the old dense search and the old SOE
           enumeration (before/after on identical inputs). Also writes the work counts.
           -> results_soe_fast/timing_v2_p0.10/{timing_budget.csv, work_counts.csv, TIMING_V2.md}
  counts   only the work counts (instant, no decoding) -> results_soe_fast/work_counts.csv

Work counts (exact, from the algorithms, independent of implementation and hardware), per decode:
  dE evaluations   SB R x S : 4 classes x 2 restarts x R x S x m      (every step evaluates all m flips per replica)
                   SA C x W : 4 classes x 2 modes    x C x W x m      (one evaluation per proposed flip)
  serial rounds    SB : S         (all m variables and all replicas are updated synchronously in one step)
                   SA : W x m     (Metropolis flips within one chain are sequential)

  python scripts/run_soe_fast.py check  --workers 500
  python scripts/run_soe_fast.py timing --timing-p 0.10
"""

from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paper1_extra as X  # noqa: E402  (also puts the package root on sys.path)
from qec_sb import cp_sb  # noqa: E402
from qec_sb.adaptive_soe import fixed_soe_prediction  # noqa: E402
from qec_sb.fast_soe import fast_soe_prediction, get_structure, SOEStructure  # noqa: E402

CHECK_SEED = 3_171_000_000                     # own namespace, disjoint from grid / extra / timing seeds
CHECK_D = (3, 5, 7, 9, 11, 13)
CHECK_P = (.06, .10, .14, .18)
CHECK_SHOTS = {3: 1000, 5: 1000, 7: 1000, 9: 1000, 11: 500, 13: 400}
CHECK_CHUNK = {3: 500, 5: 250, 7: 200, 9: 100, 11: 40, 13: 20}
TIMING_CONFIGS = tuple(X.ALL_BUDGET)           # 7 SB + 6 SA budget configs
OLD_TIMED = ("sbR8S100", "saC8W100")           # also timed with the old dense search and the old SOE enumeration
TAG = "v2"


def cpu_model():
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as f:
            for line in f:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def n_spins(d):
    return d * d - 1


def work_counts(d, cfg):
    """(dE evaluations per decode, serial rounds) for a budget config name like sbR8S400 / saC8W100."""
    m = n_spins(d)
    if cfg.startswith("sbR"):
        r, s = map(int, cfg[3:].split("S"))
        f = X.sb_cfg(r, s)
        runs = (int(f.warm_restarts) + int(f.random_restarts)) * len(tuple(f.coupling_multipliers))
        return 4 * runs * r * s * m, s
    c, w = map(int, cfg[3:].split("W"))
    return 4 * 2 * c * w * m, w * m


def counts_table(ds=X.BUDGET_D, configs=TIMING_CONFIGS):
    rows = []
    for d in ds:
        for cfg in configs:
            ev, sr = work_counts(d, cfg)
            rows.append({"d": d, "m": n_spins(d), "config": cfg, "solver": "SB" if cfg.startswith("sb") else "SA",
                         "dE_evaluations": ev, "serial_rounds": sr})
    return pd.DataFrame(rows)


def search(runner, sx, sz, key, solver_seed):
    """Seed + search exactly as run_paper1_extra.BudgetRunner.decode (same seed arithmetic), without inference."""
    code, m = runner.code, runner.method
    t0 = time.perf_counter()
    ex0, ez0 = runner.mwpm.decode(sx, sz)
    t_seed = time.perf_counter() - t0
    t0 = time.perf_counter()
    if m.startswith("sbR"):
        pools, hams, _ = X.collect_cp_class_pools_no_oracle(code, ex0, ez0, solver_seed + key * 1009
                                                            + runner.rep * X.REP_SEED_STRIDE, runner.cfg)
    else:
        seed = solver_seed + X.SA_SEED_OFFSET + key * 1009 + runner.rep * X.REP_SEED_STRIDE
        pools, hams = {}, {}
        for c in range(4):
            cex, cez = X.logical_class_seed(code, ex0, ez0, c)
            ham = X.build_stabilizer_hamiltonian(code, cex, cez)
            res = X.run_constraint_preserving_sa(ham, runner.sa, int(seed + c * 100003))
            u = np.unique(np.vstack([np.zeros((1, ham.n_spins), np.uint8), res.terminal_bits]), axis=0)
            pools[c], hams[c] = u[np.argsort(ham.energy_bits(u), kind="stable")], ham
    t_search = time.perf_counter() - t0
    return pools, hams, t_seed, t_search


def search_mode(mode, *args):
    old = cp_sb.DELTA_MODE
    cp_sb.DELTA_MODE = mode
    try:
        return search(*args)
    finally:
        cp_sb.DELTA_MODE = old


# ------------------------------------------------------------------------------------ check
def check_plan():
    tasks = []
    for d in CHECK_D:
        n, ch = CHECK_SHOTS[d], CHECK_CHUNK[d]
        for p in CHECK_P:
            for cfg in ("sbR8S100", "saC8W100") + (("sbR8S400",) if d >= 11 else ()):
                for lo in range(0, n, ch):
                    tasks.append({"d": d, "p": p, "config": cfg, "lo": lo, "hi": min(n, lo + ch)})
    cost = {3: .01, 5: .05, 7: .2, 9: .8, 11: 3.0, 13: 9.0}
    return sorted(tasks, key=lambda t: -cost[t["d"]] * (t["hi"] - t["lo"]))


def task_file(out, t):
    return Path(out) / "tasks" / f"{TAG}_d{t['d']}_p{t['p']:.2f}_{t['config']}_{t['lo']:05d}.csv"


def _pools_equal(a, b):
    return all(np.array_equal(a[c], b[c]) for c in range(4))


def run_check_task(out, t):
    path = task_file(out, t)
    if path.exists():
        return str(path), 0.0
    t_start = time.perf_counter()
    d, p = t["d"], float(t["p"])
    code = X.build_rotated_surface_code(d)
    s = X.sample_depolarizing(code, p, CHECK_SHOTS[d], CHECK_SEED + 1000 * d + int(round(p * 1000)))
    runner = X.BudgetRunner(code, p, t["config"])
    rows = []
    for i in range(t["lo"], t["hi"]):
        sx, sz = s["sx"][i], s["sz"][i]
        key = X._syn_key(sx, sz)
        pd_, _, _, _ = search_mode("dense", runner, sx, sz, key, CHECK_SEED + 77 + d)
        pools, hams, _, _ = search_mode("sparse", runner, sx, sz, key, CHECK_SEED + 77 + d)
        same_pool = int(_pools_equal(pd_, pools))
        for r in (1, 2):
            po, so, co, mo = fixed_soe_prediction(pools, hams, p, 1, r)
            pn, sn, cn, mn, hist = fast_soe_prediction(pools, hams, p, 1, r, return_hist=True)
            so, sn = np.asarray(so, float), np.asarray(sn, float)
            rel = float(np.max(np.abs(sn - so) / np.maximum(np.abs(so), 1e-300)))
            exact_tie = bool(po != pn and hist[int(po)] == hist[int(pn)])
            srt = np.sort(sn)[::-1]
            rows.append({"d": d, "p": p, "config": t["config"], "shot": i, "radius": r, "pools_identical": same_pool,
                         "pred_old": int(po), "pred_fast": int(pn), "agree": int(po == pn),
                         "exact_tie": int(exact_tie), "max_rel_diff": rel,
                         "counts_equal": int(np.array_equal(co, cn)), "minima_equal": int(np.array_equal(mo, mn)),
                         "top2_ratio": float(srt[0] / max(srt[1], 1e-300))})
    path.parent.mkdir(parents=True, exist_ok=True)
    X.write_atomic(pd.DataFrame(rows), path)
    return str(path), time.perf_counter() - t_start


def summarize_check(out):
    out = Path(out)
    files = sorted((out / "tasks").glob(f"{TAG}_*.csv"))
    if not files:
        print("no task files yet")
        return
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    g = df.groupby(["d", "p", "config", "radius"])
    S = g.agg(shots=("shot", "size"), pools_identical=("pools_identical", "min"),
              disagree=("agree", lambda x: int((x == 0).sum())), exact_tie=("exact_tie", "sum"),
              max_rel_diff=("max_rel_diff", "max"),
              counts_equal=("counts_equal", "min"), minima_equal=("minima_equal", "min")).reset_index()
    S["real_disagree"] = S.disagree - S.exact_tie
    S.to_csv(out / "check_summary.csv", index=False)
    expected = len(check_plan())
    complete = len(files) == expected
    real = int(S.real_disagree.sum())
    pools_ok = bool(S.pools_identical.min() == 1)
    ok = (complete and pools_ok and real == 0 and float(S.max_rel_diff.max()) < 1e-9
          and S.counts_equal.min() == 1 and S.minima_equal.min() == 1)
    verdict = "ALL_OK" if ok else ("INCOMPLETE" if not complete else "FAIL")
    n_pool = int(df[df.radius == 1].shape[0])
    n_pool_bad = int((df[df.radius == 1].pools_identical == 0).sum())
    tot = df.groupby("radius").agg(n=("shot", "size"), dis=("agree", lambda x: int((x == 0).sum())), tie=("exact_tie", "sum"))
    lines = [f"# {verdict}", "",
             f"Task files {len(files)}/{expected}.",
             f"Sparse vs dense search kernels: pools identical on {n_pool - n_pool_bad}/{n_pool} decodes "
             f"(CP-SB 8x100, 8x400 at d>=11, CP-SA 8x100).",
             f"Fast O(m) SOE vs enumeration on identical pools: real disagreements (not exact ties) {real}; "
             f"largest relative difference of any class score {float(S.max_rel_diff.max()):.2e}; "
             f"ball sizes and minimum weights identical: {bool(S.counts_equal.min() == 1 and S.minima_equal.min() == 1)}.", ""]
    for r, row in tot.iterrows():
        lines.append(f"- radius {r}: {int(row.n)} decodes, {int(row.dis)} decisions differ, of which {int(row.tie)} "
                     f"are exact ties (identical weight histograms, Z_a == Z_b exactly; the fast code returns the "
                     f"lowest class index, the enumeration's choice there depends on float summation order).")
    lines += ["", S.to_markdown(index=False, floatfmt=".3g"), ""]
    (out / "SOE_FAST_CHECK.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[:8]))


def cmd_check(out, workers):
    out = Path(out); (out / "tasks").mkdir(parents=True, exist_ok=True)
    tasks = [t for t in check_plan() if not task_file(out, t).exists()]
    print(f"check: {len(tasks)} tasks to run ({len(check_plan())} total), workers={workers}", flush=True)
    done = 0
    if tasks:
        with ProcessPoolExecutor(max_workers=max(1, int(workers))) as ex:
            futs = [ex.submit(run_check_task, str(out), t) for t in tasks]
            for f in as_completed(futs):
                path, dt = f.result()
                done += 1
                if done % 20 == 0 or done == len(tasks):
                    print(f"  {done}/{len(tasks)}  last {Path(path).name} {dt:.0f}s", flush=True)
    summarize_check(out)


# ----------------------------------------------------------------------------------- timing
def cmd_timing(out, shots, ps, ds, configs):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rows, agree_rows, struct_rows = [], [], []
    W = counts_table(ds, configs)
    W.to_csv(out / "work_counts.csv", index=False)
    for p in ps:
        for d in ds:
            code = X.build_rotated_surface_code(d)
            s = X.sample_depolarizing(code, p, shots, X.TIMING_SEED + d)          # same syndromes as the extra timing
            keys, seen = [], set()
            for i in range(shots):
                key = X._syn_key(s["sx"][i], s["sz"][i])
                if key in seen or key == 0:
                    continue
                seen.add(key); keys.append(i)
            solver_seed = X.TIMING_SEED + 500_000 + d
            # one-time per-code structures (generator supports, overlapping pairs): timed separately, then warm-up
            i0 = keys[0]
            warm = X.BudgetRunner(code, p, "sbR8S100")
            pools, hams, _, _ = search(warm, s["sx"][i0], s["sz"][i0], X._syn_key(s["sx"][i0], s["sz"][i0]), solver_seed)
            t0 = time.perf_counter(); S = SOEStructure(hams[0].generators, hams[0].n_data); t_struct = time.perf_counter() - t0
            t0 = time.perf_counter(); cp_sb.GeneratorSupport(hams[0].generators, hams[0].n_data); t_sup = time.perf_counter() - t0
            get_structure(hams[0]); cp_sb.generator_support(hams[0])
            fast_soe_prediction(pools, hams, p, 1, 2)
            struct_rows.append({"d": d, "p": p, "m": S.m, "overlapping_pairs": S.n_overlap_pairs,
                                "ball_r2": 1 + S.m + S.m * (S.m - 1) // 2, "soe_structure_s": t_struct,
                                "generator_support_s": t_sup})
            for cfg in configs:
                runner = X.BudgetRunner(code, p, cfg)
                rec = []
                for i in keys:
                    sx, sz = s["sx"][i], s["sz"][i]
                    key = X._syn_key(sx, sz)
                    pools, hams, t_seed, t_search = search_mode("sparse", runner, sx, sz, key, solver_seed)
                    t0 = time.perf_counter(); f1, *_ = fast_soe_prediction(pools, hams, p, 1, 1); t_f1 = time.perf_counter() - t0
                    t0 = time.perf_counter(); f2, *_ = fast_soe_prediction(pools, hams, p, 1, 2); t_f2 = time.perf_counter() - t0
                    r = {"seed": t_seed, "search": t_search, "soe_r1": t_f1, "soe_r2": t_f2}
                    if cfg in OLD_TIMED:
                        pd_, _, _, t_dense = search_mode("dense", runner, sx, sz, key, solver_seed)
                        r["search_dense"] = t_dense
                        t0 = time.perf_counter(); o1, *_ = fixed_soe_prediction(pools, hams, p, 1, 1); r["soe_r1_old"] = time.perf_counter() - t0
                        t0 = time.perf_counter(); o2, *_ = fixed_soe_prediction(pools, hams, p, 1, 2); r["soe_r2_old"] = time.perf_counter() - t0
                        agree_rows.append({"d": d, "p": p, "config": cfg, "pools_identical": int(_pools_equal(pd_, pools)),
                                           "r1_agree": int(o1 == f1), "r2_agree": int(o2 == f2)})
                    rec.append(r)
                df = pd.DataFrame(rec)
                stages = {"seed": ["seed"], "search": ["search"], "soe_r1": ["soe_r1"], "soe_r2": ["soe_r2"],
                          "end_to_end_r1": ["seed", "search", "soe_r1"], "end_to_end_r2": ["seed", "search", "soe_r2"]}
                if "soe_r2_old" in df:
                    stages.update({"search_dense": ["search_dense"], "soe_r1_old": ["soe_r1_old"],
                                   "soe_r2_old": ["soe_r2_old"],
                                   "end_to_end_r2_old": ["seed", "search_dense", "soe_r2_old"]})
                for st, cols in stages.items():
                    x = df[cols].sum(axis=1).to_numpy()
                    rows.append({"d": d, "p": p, "config": cfg, "stage": st, "syndromes": len(x),
                                 "median_s": float(np.median(x)), "p90_s": float(np.percentile(x, 90))})
                med = df.median()
                print(f"p={p} d={d} {cfg:10s} search {med['search']:.4f}s"
                      + (f" (dense {med['search_dense']:.3f}s)" if "search_dense" in med else "")
                      + f"  SOE r1 {med['soe_r1']*1e3:.2f} ms  r2 {med['soe_r2']*1e3:.2f} ms"
                      + (f" (enumeration {med['soe_r2_old']:.3f}s)" if "soe_r2_old" in med else "")
                      + f"  ({len(keys)} syndromes)", flush=True)
                pd.DataFrame(rows).to_csv(out / "timing_budget.csv", index=False)     # checkpoint
    T = pd.DataFrame(rows)
    T.to_csv(out / "timing_budget.csv", index=False)
    A, St = pd.DataFrame(agree_rows), pd.DataFrame(struct_rows)
    St.to_csv(out / "structures.csv", index=False)
    env = X.env_manifest(); env["cpu_model"] = cpu_model(); env["hostname"] = platform.node()
    lines = ["# Sparse search kernels + fast SOE: timing", "", "```json", json.dumps(env, indent=1), "```", "",
             "Same syndromes / solver seeds / pools as `run_paper1_extra.py timing`. Single process. "
             "Medians over unique non-trivial syndromes (s). Per-code structures are built once per d "
             "(`structures.csv`) and are not part of the per-syndrome time.", ""]
    piv = T[T.config.isin(OLD_TIMED)].pivot_table(index=["config", "d"], columns="stage", values="median_s")
    if len(piv):
        cols = [c for c in ("search_dense", "search", "soe_r2_old", "soe_r2", "soe_r1", "end_to_end_r2_old", "end_to_end_r2")
                if c in piv.columns]
        P = piv[cols].copy()
        if {"search_dense", "search"} <= set(P.columns):
            P["speedup_search"] = P.search_dense / P.search
        if {"soe_r2_old", "soe_r2"} <= set(P.columns):
            P["speedup_soe_r2"] = P.soe_r2_old / P.soe_r2
        lines += ["## Before (dense search, SOE enumeration) vs after (sparse, fast SOE)", "",
                  P.reset_index().to_markdown(index=False, floatfmt=".4g"), ""]
    if len(A):
        lines += [f"On these syndromes: pools identical {int(A.pools_identical.sum())}/{len(A)}, "
                  f"SOE decisions identical r1 {int(A.r1_agree.sum())}/{len(A)}, r2 {int(A.r2_agree.sum())}/{len(A)}.", ""]
    sv = T[T.stage == "search"].pivot_table(index="config", columns="d", values="median_s")
    Wd = W.pivot_table(index="config", columns="d", values="dE_evaluations")
    Wr = W.pivot_table(index="config", columns="d", values="serial_rounds")
    lines += ["## Search time (sparse) and work counts per decode", "",
              "### search median (s)", "", sv.reset_index().to_markdown(index=False, floatfmt=".4g"), "",
              "### dE evaluations", "", Wd.reset_index().to_markdown(index=False, floatfmt=".4g"), "",
              "### serial rounds", "", Wr.reset_index().to_markdown(index=False), "",
              "## Per-code structures", "", St.to_markdown(index=False, floatfmt=".4g"), "",
              "## All stages", "", T.to_markdown(index=False, floatfmt=".4g"), ""]
    (out / "TIMING_V2.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out / 'TIMING_V2.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "summarize", "timing", "counts"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--timing-shots", type=int, default=100)
    ap.add_argument("--timing-p", type=float, nargs="+", default=[.10])
    ap.add_argument("--timing-d", type=int, nargs="+", default=list(X.BUDGET_D))
    ap.add_argument("--configs", nargs="+", default=list(TIMING_CONFIGS))
    a = ap.parse_args()
    if a.command in ("check", "summarize"):
        out = a.out or "results_soe_fast/check_v2"
        cmd_check(out, a.workers) if a.command == "check" else summarize_check(out)
    elif a.command == "counts":
        out = Path(a.out or "results_soe_fast"); out.mkdir(parents=True, exist_ok=True)
        counts_table(a.timing_d, a.configs).to_csv(out / "work_counts.csv", index=False)
        print(f"wrote {out / 'work_counts.csv'}")
    else:
        out = a.out or f"results_soe_fast/timing_v2_p{a.timing_p[0]:.2f}"
        cmd_timing(out, a.timing_shots, a.timing_p, a.timing_d, a.configs)


if __name__ == "__main__":
    main()
