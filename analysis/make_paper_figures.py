"""Paper 1 figures (PDF, IEEE two-column widths) from the summary CSVs only. Runs in seconds.

Usage (from the QEC_SB_paper1 folder):
    python analysis/make_paper_figures.py --results results --out paper/figures

Inputs (all produced by code/scripts/run_paper1_{grid,extra}.py analyze / timing):
    results/paper1_grid/grid_ler.csv
    results/paper1_extra/{bposd_corr_ler,budget_ler,budget_scaled_ler,ler_vs_d_budget_scaled}.csv
    results/paper1_extra_timing_p0.10|p0.14/timing_budget.csv   (merged; p=0.10 preferred)
    results/paper1_timing/timing.csv                  (fallback for MPS / PyMatching times)
    results_soe_fast/timing_v2_p*/timing_budget.csv   (optional, code/scripts/run_soe_fast.py timing): sparse search kernels
                                                      + O(m) SOE for every budget config; replaces their seed/search/SOE/
                                                      end-to-end times (MPS, PyMatching, BP+OSD stay from the extra timing)
Outputs:
    fig_ratio_vs_p.pdf      LER/REF vs p, one panel per d (main comparison)
    fig_ler_vs_d.pdf        LER vs d: frozen budget / budget-scaled / reference
    fig_decomposition.pdf   MAP -> SOE r1 -> SOE r2, frozen (solid) and 8x400 (dashed, d>=11)
    fig_budget.pdf          d=13: accuracy vs search time (top) and vs end-to-end time (bottom)
The LaTeX source picks these up automatically (paper/main.tex falls back to the older PNGs if absent).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Fixed categorical order (same palette as the run_paper1_* analysis figures); marker shape doubles identity.
COL = {"MWPM": "#52514e", "BPc": "#1baf7a", "SA": "#eb6834", "SB": "#2a78d6", "SB400": "#6250d6",
       "REF": "#0b0b0b", "MAP": "#6250d6", "r1": "#eb6834", "r2": "#2a78d6"}
MK = {"MWPM": "s", "BPc": "D", "SA": "^", "SB": "o", "SB400": "v", "MAP": "s", "r1": "^", "r2": "o"}
W1, W2 = 3.5, 7.16  # IEEE column / page width in inches

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 6.5,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "lines.linewidth": 1.4, "lines.markersize": 3.5,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5, "pdf.fonttype": 42,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.02})


def ref_ler(grid):
    """REF LER per (d, p): the row whose method equals the cell's reference column."""
    r = grid[grid.method == grid.reference][["d", "p", "LER"]]
    return r.set_index(["d", "p"]).LER


def series(grid, method):
    return grid[grid.method == method].set_index(["d", "p"]).LER


def fig_ratio_vs_p(R, out):
    grid = pd.read_csv(R / "paper1_grid/grid_ler.csv")
    ref = ref_ler(grid)
    bpc = pd.read_csv(R / "paper1_extra/bposd_corr_ler.csv")
    bpc = bpc[bpc.method == "BP_OSD_corr"].set_index(["d", "p"]).LER
    bs = pd.read_csv(R / "paper1_extra/budget_scaled_ler.csv")
    bs = bs[bs.method == "sbR8S400_r2"].set_index(["d", "p"]).LER_over_REF   # ratio on the decoded subset
    lines = [("MWPM", series(grid, "MWPM") / ref, "MWPM / PyMatching"),
             ("BPc", bpc / ref, "BP+OSD (correlated)"),
             ("SA", series(grid, "CP_SA_local_SOE") / ref, "CP-SA 8x100 + SOE r1"),
             ("SB", series(grid, "CP_local_SOE_r2") / ref, "CP-SB 8x100 + SOE r2 (frozen)"),
             ("SB400", bs, "CP-SB 8x400 + SOE r2 (d$\\geq$11)")]
    ds = sorted(grid.d.unique())
    fig, axes = plt.subplots(2, 3, figsize=(W2, 3.9), sharex=True, sharey=True)
    for ax, d in zip(axes.flat, ds):
        for key, s, lab in lines:
            if d not in s.index.get_level_values(0):
                continue
            v = s.xs(d, level=0).sort_index()
            ls = "--" if key == "SB400" else "-"
            ax.plot(v.index, v.to_numpy(), ls, color=COL[key], marker=MK[key], label=lab)
        ax.axhline(1.0, color=COL["REF"], lw=0.8, ls=":")
        ax.set_title(f"d = {d}")
        ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("physical error rate p")
    for ax in axes[:, 0]:
        ax.set_ylabel("LER / LER$_{\\mathrm{REF}}$")
    axes[0, 0].set_yticks([1, 1.5, 2, 3, 5, 10])
    axes[0, 0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[0, 0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    h, lab = axes[1, 2].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "fig_ratio_vs_p.pdf")
    plt.close(fig)


def fig_ler_vs_d(R, out):
    F = pd.read_csv(R / "paper1_extra/ler_vs_d_budget_scaled.csv")
    cmap = plt.get_cmap("viridis")
    ps = sorted(F.p.unique())
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.5), sharey=True)
    for ax, (col, title) in zip(axes, (("frozen_main_LER", "(a) CP-SB 8x100 + SOE r2, frozen"),
                                       ("budget_scaled_LER", "(b) steps 400 for d$\\geq$11"),
                                       ("REF_LER", "(c) reference (MLD / MPS)"))):
        for k, p in enumerate(ps):
            s = F[F.p == p].sort_values("d")
            ax.plot(s.d, s[col].clip(lower=1e-7), "-o", ms=2.5, lw=1.1, color=cmap(k / max(len(ps) - 1, 1)),
                    label=f"{p:.2f}")
        ax.set_yscale("log")
        ax.set_xlabel("code distance d")
        ax.set_title(title)
        ax.set_xticks(sorted(F.d.unique()))
    axes[0].set_ylabel("logical error rate")
    h, lab = axes[2].get_legend_handles_labels()
    fig.legend(h, lab, title="p", loc="center left", bbox_to_anchor=(1.0, 0.55), fontsize=6.5, title_fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_ler_vs_d.pdf")
    plt.close(fig)


def fig_decomposition(R, out):
    grid = pd.read_csv(R / "paper1_grid/grid_ler.csv")
    ref = ref_ler(grid)
    E = pd.read_csv(R / "paper1_extra/budget_scaled_ler.csv")
    ds = sorted(grid.d.unique())
    fig, axes = plt.subplots(2, 3, figsize=(W2, 3.7), sharex=True, sharey=True)
    for ax, d in zip(axes.flat, ds):
        for m, key, lab in (("CP_SB_MAP", "MAP", "MAP (search only)"), ("CP_local_SOE", "r1", "+ SOE r1"),
                            ("CP_local_SOE_r2", "r2", "+ SOE r2")):
            v = (series(grid, m) / ref).xs(d, level=0).sort_index()
            ax.plot(v.index, v.to_numpy(), "-", color=COL[key], marker=MK[key], label=lab + ", 8x100")
        e = E[E.d == d]
        for inf, key in (("MAP", "MAP"), ("r1", "r1"), ("r2", "r2")):
            s = e[e.method == f"sbR8S400_{inf}"].sort_values("p")
            if len(s):
                ax.plot(s.p, s.LER_over_REF, "--", color=COL[key], marker=MK[key], mfc="white",
                        label=f"{'MAP' if inf == 'MAP' else '+ SOE ' + inf}, 8x400")
        ax.axhline(1.0, color=COL["REF"], lw=0.8, ls=":")
        ax.set_title(f"d = {d}")
        ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("physical error rate p")
    for ax in axes[:, 0]:
        ax.set_ylabel("LER / LER$_{\\mathrm{REF}}$")
    axes[0, 0].set_yticks([1, 1.5, 2, 3])
    axes[0, 0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[0, 0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    h, lab = axes[1, 2].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "fig_decomposition.pdf")
    plt.close(fig)


def _prefer_p010(T):
    # prefer p=0.10 (both timing campaigns have it); otherwise the nearest p available per (d, config, stage)
    T = T.copy()
    T["pri"] = (T.p - 0.10).abs()
    return T.sort_values("pri", kind="stable").drop_duplicates(["d", "config", "stage"], keep="first")


def load_timing(R, F=None):
    frames = [pd.read_csv(f) for f in sorted(R.glob("paper1_extra_timing_p*/timing_budget.csv"))]
    T = _prefer_p010(pd.concat(frames, ignore_index=True)).set_index(["d", "config", "stage"]).median_s.to_dict()
    # v1 (timing_p*): fast SOE only -> replace SOE / end-to-end stages.
    # v2 (timing_v2_p*): sparse search kernels + fast SOE for every budget config -> replace all CP stages.
    for pattern, stages in (("timing_p*/timing_budget.csv", ["soe_r1", "soe_r2", "end_to_end_r1", "end_to_end_r2"]),
                            ("timing_v2_p*/timing_budget.csv", ["seed", "search", "soe_r1", "soe_r2",
                                                                "end_to_end_r1", "end_to_end_r2"])):
        fast = sorted(Path(F).glob(pattern)) if F is not None else []
        if not fast:
            continue
        N = _prefer_p010(pd.concat([pd.read_csv(f) for f in fast], ignore_index=True))
        N = N[N.stage.isin(stages)]
        for d, c, st, v in zip(N.d, N.config, N.stage, N.median_s):
            T[(int(d), c, st)] = float(v)
        print(f"     {pattern.split('/')[0]}: {len(N.config.unique())} configs, stages {stages}")
    return pd.Series(T)


def baseline_time(T, R, d, names):
    for n in names:
        for st in ("total", "decode"):
            if (d, n, st) in T.index:
                return float(T[(d, n, st)])
    t = pd.read_csv(R / "paper1_timing/timing.csv")
    t = t[(t.d == d) & (t.method.isin(names))]
    return float(t.median_s.iloc[0]) if len(t) else np.nan


def fig_budget(R, out, d=13, ps=(0.12, 0.16), F=None):
    B = pd.read_csv(R / "paper1_extra/budget_ler.csv")
    T = load_timing(R, F)
    bpc = pd.read_csv(R / "paper1_extra/bposd_corr_ler.csv")
    groups = [("SB", "SB: steps 100$\\to$400 (8 rep.)", ["sbR8S100", "sbR8S200", "sbR8S400"], "-", "o"),
              ("SB", "SB: replicas 8$\\to$64 (100 steps)", ["sbR8S100", "sbR16S100", "sbR32S100", "sbR64S100"], ":", "s"),
              ("SB", "SB 32x400", ["sbR32S400"], "", "*"),
              ("SA", "SA: sweeps 25$\\to$200 (8 chains)", ["saC8W25", "saC8W50", "saC8W100", "saC8W200"], "-", "^"),
              ("SA", "SA 16x100, 32x200", ["saC16W100", "saC32W200"], "", "P")]
    fig, axes = plt.subplots(2, len(ps), figsize=(W2, 4.6), sharey="row")
    for j, p in enumerate(ps):
        b = B[(B.d == d) & np.isclose(B.p, p)].set_index(["config", "inference"]).LER_over_REF
        ax = axes[0, j]
        for key, lab, cfgs, ls, mk in groups:
            x = [T.get((d, c, "search"), np.nan) for c in cfgs]
            y = [b.get((c, "r1"), np.nan) for c in cfgs]
            ax.plot(x, y, ls if ls else "none", marker=mk, ms=5 if mk == "*" else 3.5, color=COL[key], label=lab)
        ax.axhline(1.0, color=COL["REF"], lw=0.8, ls=":")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"(a{j + 1}) d={d}, p={p:.2f}: SOE r1 vs search time")
        ax.set_xlabel("search time per syndrome (s, 1 core)")
        ax = axes[1, j]
        for c, key, mk, lab in (("sbR8S400", "SB", "o", "CP-SB 8x400"), ("sbR32S400", "SB", "*", "CP-SB 32x400"),
                                ("saC8W100", "SA", "^", "CP-SA 8x100"), ("saC32W200", "SA", "P", "CP-SA 32x200")):
            t1 = T.get((d, c, "end_to_end_r1"), np.nan)
            if not np.isfinite(t1):   # older timing runs had no end_to_end_r1 stage
                t1 = T.get((d, c, "seed"), 0.0) + T.get((d, c, "search"), np.nan) + T.get((d, c, "soe_r1"), np.nan)
            t2 = T.get((d, c, "end_to_end_r2"), np.nan)
            ax.plot([t1, t2], [b.get((c, "r1"), np.nan), b.get((c, "r2"), np.nan)], "-", marker=mk,
                    ms=5 if mk == "*" else 3.5, color=COL[key], label=lab + " (r1$\\to$r2)")
        tm = baseline_time(T, R, d, ["mps16"])
        ax.plot([tm], [1.0], "k*", ms=7, label="MPS $\\chi$=16 (reference class)")
        tb = baseline_time(T, R, d, ["bposd3", "bposd_corr"])
        cell = bpc[(bpc.d == d) & np.isclose(bpc.p, p)].set_index("method").LER
        if "BP_OSD_corr" in cell.index and "REF" in cell.index:
            ax.plot([tb], [cell["BP_OSD_corr"] / cell["REF"]], "D", color=COL["BPc"], label="BP+OSD (correlated)")
        ax.axhline(1.0, color=COL["REF"], lw=0.8, ls=":")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"(b{j + 1}) d={d}, p={p:.2f}: end-to-end time")
        ax.set_xlabel("end-to-end time per syndrome (s, 1 core)")
    for ax in axes[:, 0]:
        ax.set_ylabel("LER / LER$_{\\mathrm{REF}}$")
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[0, 0].set_yticks([1, 1.25, 1.5, 2, 3])
    axes[1, 0].set_yticks([1, 1.25, 1.5, 2, 3])
    axes[0, -1].legend(loc="upper right", fontsize=6)
    axes[1, -1].legend(loc="upper right", fontsize=6)
    fig.tight_layout()
    fig.savefig(out / "fig_budget.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="paper/figures")
    ap.add_argument("--soe-fast", default="results_soe_fast", help="folder from run_soe_fast.py (optional)")
    a = ap.parse_args()
    R, out = Path(a.results), Path(a.out)
    F = Path(a.soe_fast) if Path(a.soe_fast).is_dir() else None
    out.mkdir(parents=True, exist_ok=True)
    for f in (fig_ratio_vs_p, fig_ler_vs_d, fig_decomposition, fig_budget):
        try:
            f(R, out, F=F) if f is fig_budget else f(R, out)
            print("ok  ", f.__name__)
        except Exception as e:  # one broken input must not block the other figures
            print("FAIL", f.__name__, repr(e))


if __name__ == "__main__":
    main()
