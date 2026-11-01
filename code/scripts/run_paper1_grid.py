"""Paper 1 (v2) p x d grid: code-capacity depolarizing, rotated surface code, CPU only.

Held-out evaluation. Run it only AFTER ``run_paper1_v2.py analyze`` has written
``main_method.json`` (the main decoder is chosen on the confirmation set by a
pre-declared rule and is not changed afterwards).

New, independent shots (seed namespace 2029000000+, disjoint from the confirmation
family). Every method decodes the same shots, so all comparisons are paired.

  d = 3,5,7,9,11,13
  p = .05,.06,.08,.10,.12,.14,.16,.18   (d>=11: p>=.08 only; low-p cells would have ~0 failures)
  methods: frozen QEC-SB stages (MWPM, CP-SB MAP, local SOE r1, importance, Bennett),
           local SOE r2, PyMatching, BP+OSD, CP-SA (MAP, local SOE), MPS chi=16 (all shots),
           MPS chi=32 (all shots for d<=9; hard shots only for d>=11), exact MLD (d=3)
  reference: exact MLD (d=3); MPS chi=32 (d<=9); MPS chi=16 with chi=32 on hard shots (d>=11)

Subcommands
  estimate   time a few shots of every (d, method) and print the projected cost
  run        all decoding tasks (resumable; finished chunks are skipped)
  hardchi    d>=11: re-decode the hard shots (any decoder wrong or disagreeing) with MPS chi=32
  analyze    LER, paired tests, gap closure, effective distance, chi check, figures
  timing     single-process latency on an otherwise idle machine (paper timing numbers)
  all        run -> hardchi -> analyze

  python scripts/run_paper1_grid.py estimate --workers 1024
  python scripts/run_paper1_grid.py all --out results/paper1_grid --workers 1024
  python scripts/run_paper1_grid.py timing --out results/paper1_timing      # on an idle node
"""

from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_paper1_v2 import (MethodRunner, decode_samples, env_manifest, write_atomic, gap_closure,  # noqa: E402
                           build_rotated_surface_code, sample_depolarizing, sample_sha256,
                           build_frozen_p42a_config, run_p42a, exact_paired_statistics,
                           wilson_interval, adjust_p_values, effective_distance_table, _syn_key)

DISTANCES = (3, 5, 7, 9, 11, 13)
ORIG_PROBS = (.05, .06, .08, .10, .12, .14, .16, .18)   # seeds of these cells never change
LOWP_EXT = os.environ.get("QECSB_LOWP_EXT", "0") == "1"  # overnight low-p extension (see shots_for)
PROBS = ORIG_PROBS + ((.04,) if LOWP_EXT else ())
GRID_BASE = 2_029_000_000
TIMING_SEED = 2_929_000_000
CHUNK = {3: 5000, 5: 2000, 7: 1000, 9: 500, 11: 250, 13: 200}
STAGES = ("MWPM", "CP_SB_MAP", "CP_local_SOE", "P4_1_importance", "P4_2_Bennett")
HARD_CAP = 3000          # max hard shots per cell re-decoded with chi=32 (d>=11)
SCALE = 1.0              # --scale, smoke tests only


def is_ext_cell(d, p):
    """Cells added by the low-p extension: p=.04 for every d, and p=.05/.06 for d>=11."""
    return p == .04 or (d >= 11 and p < .08)


def shots_for(d, p):
    if LOWP_EXT and is_ext_cell(d, p):
        return 1_000_000 if (p == .04 and d >= 11) else 400_000
    base = {.05: 200_000, .06: 100_000, .08: 50_000, .10: 20_000, .12: 20_000}.get(p, 10_000)
    if d >= 11:
        return 0 if p < .08 else min(base, 50_000)
    return base


def chunk_for(d, p):
    if LOWP_EXT and is_ext_cell(d, p):       # bigger chunks keep the file count manageable
        return {3: 5000, 5: 5000, 7: 4000, 9: 2000, 11: 1000, 13: 1000}[d]
    return CHUNK[d]


def cell_seed(d, p):
    if p in ORIG_PROBS:                      # unchanged formula for the original 8 probabilities
        return GRID_BASE + 1_000_000 * (DISTANCES.index(d) * len(ORIG_PROBS) + ORIG_PROBS.index(p))
    ext = (.04,)                             # separate namespace, disjoint from every other seed in use
    return GRID_BASE + 100_000_000 + 1_000_000 * (DISTANCES.index(d) * 10 + ext.index(p))


def methods_for(d):
    m = ["ours", "r2", "pym", "bposd", "sa", "mps16"]
    if d <= 9:
        m.append("mps32")
    if d == 3:
        m.append("mld")
    return m


def cell_id(d, p):
    return f"d{d}_p{p:.3f}"


def plan():
    tasks = []
    for d in DISTANCES:
        for p in PROBS:
            n = shots_for(d, p)
            if n == 0:
                continue
            n = max(1, int(n * SCALE))
            ch = chunk_for(d, p)
            for k, lo in enumerate(range(0, n, ch)):
                hi = min(lo + ch, n)
                for m in methods_for(d):
                    tasks.append({"d": d, "p": p, "chunk": k, "lo": lo, "hi": hi, "method": m,
                                  "sample_seed": cell_seed(d, p) + k, "solver_seed": cell_seed(d, p) + 500_000 + k})
    return tasks


def chunk_path(out, t):
    return Path(out) / "chunks" / cell_id(t["d"], t["p"]) / f"c{t['chunk']:05d}__{t['method']}.csv"


def run_task(out, t, limit=None):
    """Decode one chunk with one method. Shots of a chunk are sampled from its own seed."""
    path = chunk_path(out, t)
    if path.exists() and limit is None:
        return f"{path.parent.name}/{path.name} cached"
    d, p, n = t["d"], t["p"], t["hi"] - t["lo"]
    if limit:
        n = min(n, limit)
    code = build_rotated_surface_code(d)
    samples = sample_depolarizing(code, p, n, t["sample_seed"])
    extra = {"d": d, "p": p, "chunk": t["chunk"], "sample_sha256": sample_sha256(samples)}
    t0 = time.time()
    if t["method"] == "ours":
        tmp = Path(out) / "tmp" / f"{cell_id(d, p)}_c{t['chunk']}_{os.getpid()}"
        shutil.rmtree(tmp, ignore_errors=True)
        run_p42a(tmp, build_frozen_p42a_config(d, p, n, t["sample_seed"], t["solver_seed"]))
        shot = pd.read_csv(tmp / "shot_level.csv", dtype={"syndrome_key": str})
        syn = pd.read_csv(tmp / "syndrome_level.csv", dtype={"syndrome_key": str}).set_index("syndrome_key")
        t_full = syn[["time_mwpm_s", "time_cp_s", "time_local_s", "time_p41_s", "time_bar_s"]].sum(axis=1)
        t_soe = syn[["time_mwpm_s", "time_cp_s", "time_local_s"]].sum(axis=1)   # MWPM seed + CP-SB + local SOE r1
        df = shot[["shot", "syndrome_key", "truth_class"] + [f"pred_{s}" for s in STAGES]].copy()
        df["error_weight"] = np.count_nonzero(samples["pauli"], axis=1)
        df["time_ours_s"] = df.syndrome_key.map(t_full).to_numpy()
        df["time_ours_soe_s"] = df.syndrome_key.map(t_soe).to_numpy()
        for k, v in extra.items():
            df[k] = v
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        runner = MethodRunner(code, p, t["method"])
        df = pd.DataFrame(decode_samples(samples, 0, n, runner, t["solver_seed"], extra))
    if limit:
        return df
    path.parent.mkdir(parents=True, exist_ok=True)
    write_atomic(df, path)
    return f"{path.parent.name}/{path.name} {time.time() - t0:.0f}s"


def _cost_key(t):
    return -(t["d"] ** 3) * (t["hi"] - t["lo"]) * {"sa": 3, "ours": 2, "r2": 2, "mps32": 2}.get(t["method"], 1)


def run_parallel(fn, jobs, workers, label):
    print(f"{label}: {len(jobs)} tasks with {workers} workers", flush=True)
    if not jobs:
        return
    failures, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, *args): args for args in jobs}
        for k, f in enumerate(as_completed(futs), 1):
            try:
                msg = f.result()
                if k % 50 == 0 or k == len(futs):
                    print(f"[{k}/{len(futs)}] {time.time() - t0:.0f}s  {msg}", flush=True)
            except Exception as e:
                failures.append((futs[f][1:], repr(e)))
                print(f"FAILED {futs[f][1:]}: {e!r}", flush=True)
    if failures:
        raise SystemExit(f"{len(failures)} task(s) failed, e.g. {failures[0]}")


SHARD, NSHARDS = 0, 1    # --shard/--nshards: split tasks over Slurm array jobs (round-robin on the cost-sorted list)


def cmd_run(out, workers):
    tasks = sorted(plan(), key=_cost_key)[SHARD::NSHARDS]
    tasks = [t for t in tasks if not chunk_path(out, t).exists()]
    run_parallel(run_task, [(out, t) for t in tasks], workers, f"run shard {SHARD}/{NSHARDS}")


# ------------------------------------------------------------------ hard-shot chi check
def load_cell(cdir):
    by_chunk = {}
    for f in cdir.glob("c*__*.csv"):
        by_chunk.setdefault(int(f.name[1:6]), []).append(pd.read_csv(f, dtype={"syndrome_key": str}))
    frames = []
    for k, dfs in sorted(by_chunk.items()):
        base = dfs[0][["d", "p", "chunk", "shot", "syndrome_key", "truth_class", "error_weight", "sample_sha256"]].copy()
        for df in dfs:
            if not (np.array_equal(df.truth_class, base.truth_class) and df.sample_sha256.iloc[0] == base.sample_sha256.iloc[0]):
                raise ValueError(f"{cdir.name} chunk {k}: methods disagree on shots")
            cols = [c for c in df.columns if c.startswith(("pred_", "time_"))]
            base = pd.concat([base, df[cols].reset_index(drop=True)], axis=1)
        frames.append(base)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def hard_mask(g):
    """Shots where chi=16 MPS is wrong, or where any other decoder disagrees with it."""
    preds = [c for c in g.columns if c.startswith("pred_") and c not in
             ("pred_MPS_chi16", "pred_MPS_chi32_hard", "pred_MPS_ref")]
    m = g.pred_MPS_chi16 != g.truth_class
    for c in preds:
        m |= g[c] != g.pred_MPS_chi16
    return m


def hardchi_task(out, d, p, chunk, shots):
    path = Path(out) / "hard" / cell_id(d, p) / f"c{chunk:05d}__mps32.csv"
    if path.exists():
        return f"{path.name} cached"
    t = next(t for t in plan() if t["d"] == d and t["p"] == p and t["chunk"] == chunk and t["method"] == "mps16")
    code = build_rotated_surface_code(d)
    s = sample_depolarizing(code, p, t["hi"] - t["lo"], t["sample_seed"])
    runner = MethodRunner(code, p, "mps32")
    rows = []
    for i in shots:
        pred = runner.decode(s["sx"][i], s["sz"][i], _syn_key(s["sx"][i], s["sz"][i]), 0)
        rows.append({"chunk": chunk, "shot": int(i), **pred})
    path.parent.mkdir(parents=True, exist_ok=True)
    write_atomic(pd.DataFrame(rows), path)
    return f"{cell_id(d, p)}/{path.name} {len(shots)} hard shots"


def cmd_hardchi(out, workers):
    jobs = []
    for cdir in sorted((Path(out) / "chunks").glob("d*_p*")):
        g = load_cell(cdir)
        if not len(g) or int(g.d.iloc[0]) < 11:
            continue
        d, p = int(g.d.iloc[0]), float(g.p.iloc[0])
        hard = g[hard_mask(g)]
        if len(hard) > HARD_CAP:
            hard = hard.sample(HARD_CAP, random_state=0)
        for chunk, h in hard.groupby("chunk"):
            jobs.append((out, d, p, int(chunk), sorted(h.shot.tolist())))
    run_parallel(hardchi_task, jobs[SHARD::NSHARDS], workers, f"hardchi shard {SHARD}/{NSHARDS}")


# ---------------------------------------------------------------------------- analyze
def load_grid(out):
    frames = []
    for cdir in sorted((Path(out) / "chunks").glob("d*_p*")):
        g = load_cell(cdir)
        if not len(g):
            continue
        hdir = Path(out) / "hard" / cdir.name
        hf = sorted(hdir.glob("c*__mps32.csv")) if hdir.exists() else []
        if hf:
            g["hard_shot"] = hard_mask(g)
            h = pd.concat([pd.read_csv(f) for f in hf], ignore_index=True)
            h = h.rename(columns={"pred_MPS_chi32": "pred_MPS_chi32_hard"})[["chunk", "shot", "pred_MPS_chi32_hard"]]
            g = g.merge(h, on=["chunk", "shot"], how="left")
            g["hard_checked"] = g.pred_MPS_chi32_hard.notna()
            g["pred_MPS_ref"] = g.pred_MPS_chi32_hard.fillna(g.pred_MPS_chi16).astype(int)
        frames.append(g)
    df = pd.concat(frames, ignore_index=True)
    for c in [c for c in df.columns if c.startswith("pred_")]:
        df["correct_" + c[5:]] = np.where(df[c].isna(), np.nan, (df[c] == df.truth_class).astype(float))
    return df


def reference_for(d, g):
    if d == 3:
        return "MLD_exact"
    if "correct_MPS_chi32" in g and g["correct_MPS_chi32"].notna().all():
        return "MPS_chi32"
    if "correct_MPS_ref" in g and g["correct_MPS_ref"].notna().all():
        return "MPS_ref"
    return "MPS_chi16"


def load_main(main, main_from):
    if main:
        return main
    f = Path(main_from)
    if not f.exists():
        raise SystemExit(f"{f} not found: run `run_paper1_v2.py analyze` first, or pass --main")
    return json.loads(f.read_text(encoding="utf-8"))["main"]


def cmd_analyze(out, main):
    out = Path(out)
    df = load_grid(out)
    methods = ["MWPM", "PyMatching", "BP_OSD", "CP_SA_MAP", "CP_SB_MAP", "CP_SA_local_SOE", "CP_local_SOE",
               "CP_local_SOE_r2", "P4_1_importance", "P4_2_Bennett", "MPS_chi16", "MPS_chi32", "MPS_ref", "MLD_exact"]
    ler, pairs, gaps, chi, repro = [], [], [], [], []
    for (d, p), g in df.groupby(["d", "p"]):
        ref = reference_for(d, g)
        for m in methods:
            if f"correct_{m}" not in g or g[f"correct_{m}"].isna().all():
                continue
            c = g[f"correct_{m}"].dropna()
            k, n = int(len(c) - c.sum()), len(c)
            lo, hi = wilson_interval(k, n)
            ler.append({"d": d, "p": p, "method": m, "shots": n, "failures": k, "LER": k / n,
                        "LER_ci95_lo": lo, "LER_ci95_hi": hi, "reference": ref, "is_main": m == main})
        comps = [(main, "MWPM"), (main, "PyMatching"), (main, "BP_OSD"), (main, ref),
                 ("CP_SA_local_SOE", "CP_local_SOE"), ("CP_SA_MAP", "CP_SB_MAP"),
                 ("CP_local_SOE_r2", "CP_local_SOE"), ("CP_local_SOE", "CP_SB_MAP"),
                 ("BP_OSD", ref), ("PyMatching", ref), ("CP_SA_local_SOE", ref)]
        for a, b in dict.fromkeys(comps):
            if a == b or f"correct_{a}" not in g or f"correct_{b}" not in g:
                continue
            st = exact_paired_statistics(g[f"correct_{a}"].to_numpy(np.uint8), g[f"correct_{b}"].to_numpy(np.uint8))
            pairs.append({"d": d, "p": p, "a": a, "b": b, "family": f"{a} vs {'REF' if b == ref else b}",
                          "LER_a_minus_LER_b": -st["paired_difference"],
                          "ci95_lo": -st["paired_ci95_hi"], "ci95_hi": -st["paired_ci95_lo"],
                          "a_right_b_wrong": st["discordant_wins"], "a_wrong_b_right": st["discordant_losses"],
                          "exact_p": st["exact_p_value"]})
        for a in dict.fromkeys((main, "CP_local_SOE", "CP_local_SOE_r2", "CP_SA_local_SOE", "BP_OSD", "PyMatching")):
            val, lo, hi = gap_closure(g.correct_MWPM.to_numpy(), g[f"correct_{a}"].to_numpy(),
                                      g[f"correct_{ref}"].to_numpy(), reps=1000)
            gaps.append({"d": d, "p": p, "method": a, "reference": ref, "gap_closed": val, "ci95_lo": lo, "ci95_hi": hi})
        if "pred_MPS_chi32" in g and g.pred_MPS_chi32.notna().any():
            sub = g[g.pred_MPS_chi32.notna()]
            chi.append({"d": d, "p": p, "subset": "all shots", "shots": len(sub),
                        "chi16_vs_chi32_disagreements": int((sub.pred_MPS_chi16 != sub.pred_MPS_chi32).sum())})
        if "hard_checked" in g and g.hard_checked.notna().any():
            sub = g[g.hard_checked == True]  # noqa: E712
            chi.append({"d": d, "p": p, "subset": f"hard shots ({int(g.hard_shot.sum())} hard, cap {HARD_CAP})",
                        "shots": len(sub),
                        "chi16_vs_chi32_disagreements": int((sub.pred_MPS_chi16 != sub.pred_MPS_chi32_hard).sum())})
        if "pred_CP_SB_MAP_rebuilt" in g:
            repro.append({"d": d, "p": p, "r2_rebuilt_pools_match_frozen":
                          bool((g.pred_CP_SB_MAP_rebuilt == g.pred_CP_SB_MAP).all()
                               and (g.pred_CP_local_SOE_rebuilt == g.pred_CP_local_SOE).all())})
    pairs = pd.DataFrame(pairs)
    for fam, idx in pairs.groupby("family").groups.items():
        pairs.loc[idx, "holm_p"] = adjust_p_values(pairs.loc[idx, "exact_p"].tolist(), "holm")
    present = [m for m in methods if f"correct_{m}" in df]
    tables = {"grid_ler": pd.DataFrame(ler), "grid_paired": pairs, "grid_gap_closure": pd.DataFrame(gaps),
              "grid_mps_chi_check": pd.DataFrame(chi), "grid_r2_pool_check": pd.DataFrame(repro),
              "grid_effective_distance": effective_distance_table(df, present, "d")}
    for name, t in tables.items():
        t.to_csv(out / f"{name}.csv", index=False)
    md = ["# Paper 1 grid results (code capacity, depolarizing)", "",
          f"Main method (fixed on the confirmation set): **{main}**.",
          "Reference: exact MLD (d=3); MPS chi=32 (d<=9); MPS chi=16 with chi=32 on hard shots (d>=11).",
          "Holm correction within each comparison family across all (d, p) cells.",
          "Timing is NOT reported here (server under full load); see the `timing` subcommand.", ""]
    for name, t in tables.items():
        md += [f"## {name}", "", t.to_markdown(index=False, floatfmt=".4g") if len(t) else "(empty)", ""]
    (out / "RESULTS_GRID.md").write_text("\n".join(md), encoding="utf-8")
    plot_vs_p(tables["grid_ler"], main, out / "grid_ler_vs_p.pdf")
    plot_vs_d(tables["grid_ler"], main, out / "grid_ler_vs_d.pdf")
    print(f"wrote {out / 'RESULTS_GRID.md'} and figures")


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _err(ax, x, s, **kw):
    y = s.LER.clip(lower=1e-7)
    ax.errorbar(x, y, yerr=[(y - s.LER_ci95_lo.clip(lower=1e-7)).clip(lower=0), (s.LER_ci95_hi.clip(lower=1e-7) - y).clip(lower=0)],
                marker="o", ms=3, capsize=2, **kw)


def plot_vs_p(ler, main, path):
    plt = _plt()
    show = ["PyMatching", "BP_OSD", "CP_SA_local_SOE", main, "reference"]
    ds = sorted(ler.d.unique())
    fig, axes = plt.subplots(1, len(ds), figsize=(3.2 * len(ds), 3.2), sharey=True)
    for ax, d in zip(np.atleast_1d(axes), ds):
        sub = ler[ler.d == d]
        for m in show:
            if m == "reference":
                s = sub[sub.method == sub.reference].sort_values("p")
                lab = "reference (MLD / MPS)"
            else:
                s, lab = sub[sub.method == m].sort_values("p"), m
            if len(s):
                _err(ax, s.p, s, label=lab)
        ax.set_yscale("log"); ax.set_xlabel("p"); ax.set_title(f"d={d}")
    np.atleast_1d(axes)[0].set_ylabel("logical error rate")
    np.atleast_1d(axes)[-1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=150)


def plot_vs_d(ler, main, path):
    """LER vs d at fixed p: threshold crossing and sub-threshold suppression."""
    plt = _plt()
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.3), sharey=True)
    for ax, m in zip(axes, [main, "reference", "PyMatching"]):
        for p, s in ler.groupby("p"):
            s = (s[s.method == s.reference] if m == "reference" else s[s.method == m]).sort_values("d")
            if len(s):
                _err(ax, s.d, s, label=f"p={p:.2f}")
        ax.set_yscale("log"); ax.set_xlabel("d"); ax.set_title(m)
    axes[0].set_ylabel("logical error rate"); axes[-1].legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=150)


# ---------------------------------------------------------------------------- timing
def cmd_timing(out, shots, p, methods):
    """Single process, sequential, per unique non-trivial syndrome. Run on an idle node."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for d in DISTANCES:
        for m in methods:
            if m == "mld" and d != 3:
                continue
            t = {"d": d, "p": p, "chunk": 0, "lo": 0, "hi": shots, "method": m,
                 "sample_seed": TIMING_SEED + d, "solver_seed": TIMING_SEED + 500_000 + d}
            df = run_task(str(out), t, limit=shots)
            u = df.drop_duplicates("syndrome_key")
            u = u[u.syndrome_key.astype(str) != "0"]
            cols = ["time_ours_s", "time_ours_soe_s"] if m == "ours" else [f"time_{m}_s"]
            for tc in cols:
                x = u[tc].to_numpy(float)
                rows.append({"d": d, "p": p, "method": tc[5:-2], "unique_nontrivial_syndromes": len(x),
                             "median_s": float(np.median(x)), "p90_s": float(np.percentile(x, 90)),
                             "mean_s": float(x.mean())})
                print(f"d={d:2d} {tc[5:-2]:9s} median {rows[-1]['median_s']*1e3:9.2f} ms  "
                      f"p90 {rows[-1]['p90_s']*1e3:9.2f} ms", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / "timing.csv", index=False)
    (out / "TIMING.md").write_text(
        "# Single-process CPU wall time per unique non-trivial syndrome\n\n"
        "`ours` = full frozen pipeline; `ours_soe` = MWPM seed + CP-SB + local SOE r1; "
        "`r2` = MWPM seed + CP-SB pool rebuild + SOE r1 and r2.\n\n```json\n"
        + json.dumps(env_manifest(), indent=1) + "\n```\n\n" + t.to_markdown(index=False, floatfmt=".4g") + "\n",
        encoding="utf-8")


def cmd_estimate(out, workers, shots=12):
    rows = []
    for d in DISTANCES:
        for m in methods_for(d):
            t = {"d": d, "p": .12, "chunk": 0, "lo": 0, "hi": shots, "method": m, "sample_seed": 1, "solver_seed": 2}
            t0 = time.time()
            run_task(out, t, limit=shots)
            per = (time.time() - t0) / shots
            total = sum(t2["hi"] - t2["lo"] for t2 in plan() if t2["d"] == d and t2["method"] == m)
            rows.append({"d": d, "method": m, "s_per_shot": per, "planned_shots": total, "core_hours": per * total / 3600})
            print(f"d={d:2d} {m:6s} {per:8.4f} s/shot  x {total:>9,d} shots = {per * total / 3600:7.1f} core-h", flush=True)
    df = pd.DataFrame(rows)
    ch = df.core_hours.sum()
    print(f"\nTotal ~{ch:.0f} core-hours (timed at p=.12; p>=.14 cells are slower, low-p cells faster).")
    print(f"With {workers} workers: ~{ch / workers * 60:.0f} min if perfectly parallel; "
          f"longest single task ~{df['s_per_shot'].max() * max(CHUNK.values()) / 60:.0f} min. hardchi adds little.")


def main():
    global SCALE, SHARD, NSHARDS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["estimate", "run", "hardchi", "analyze", "timing", "all"])
    ap.add_argument("--out", default="results/paper1_grid")
    ap.add_argument("--workers", type=int, default=min(1024, max(1, (os.cpu_count() or 2) - 2)))
    ap.add_argument("--scale", type=float, default=1.0, help="multiply all shot counts (smoke test only)")
    ap.add_argument("--main", default=None, help="override the main method (default: read --main-from)")
    ap.add_argument("--main-from", default="results/paper1_v2/main_method.json")
    ap.add_argument("--timing-shots", type=int, default=300)
    ap.add_argument("--timing-p", type=float, default=.10)
    ap.add_argument("--timing-methods", nargs="+", default=["ours", "r2", "pym", "bposd", "sa", "mps16"])
    ap.add_argument("--shard", type=int, default=0, help="this job's index (Slurm array: $SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--nshards", type=int, default=1, help="number of array jobs sharing the tasks")
    a = ap.parse_args()
    SCALE, SHARD, NSHARDS = a.scale, a.shard, a.nshards
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.command == "timing":
        cmd_timing(out, a.timing_shots, a.timing_p, a.timing_methods)
        return
    main_method = load_main(a.main, a.main_from) if a.command in ("analyze", "all") else None
    m = out / "run_manifest.json"
    if not m.exists():
        m.write_text(json.dumps({"env": env_manifest(), "grid_base_seed": GRID_BASE, "distances": DISTANCES,
                                 "probabilities": PROBS, "scale": a.scale, "main_method": main_method,
                                 "shots": {cell_id(d, p): shots_for(d, p) for d in DISTANCES for p in PROBS}},
                                indent=2), encoding="utf-8")
    if a.command == "estimate":
        cmd_estimate(out, a.workers)
    if a.command in ("run", "all"):
        cmd_run(out, a.workers)
    if a.command in ("hardchi", "all"):
        cmd_hardchi(out, a.workers)
    if a.command in ("analyze", "all"):
        cmd_analyze(out, main_method)


if __name__ == "__main__":
    main()
