"""Paper 1 post-run analysis (reads only the summary CSVs; runs in seconds).

Usage (from paper1/results):
    python analysis/analyze_paper1.py --root . --out analysis
Inputs : paper1_v2/{ler_table,paired_comparisons,gap_closure,effective_distance,mps_chi_convergence}.csv,
         paper1_v2/main_method.json, paper1_grid/{grid_ler,grid_paired,grid_gap_closure,
         grid_effective_distance,grid_mps_chi_check}.csv
Outputs: analysis/ANALYSIS.md, analysis/*.csv, analysis/fig_*.pdf/png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# reference categorical palette (dataviz skill), fixed order, never cycled
C = {"main": "#2a78d6", "CP_local_SOE": "#eb6834", "CP_SA_local_SOE": "#1baf7a", "BP_OSD": "#eda100",
     "MWPM": "#52514e", "REF": "#0b0b0b"}
LABEL = {"CP_local_SOE_r2": "CP-SB + SOE r2 (main)", "CP_local_SOE": "CP-SB + SOE r1",
         "CP_SA_local_SOE": "CP-SA + SOE r1", "BP_OSD": "BP+OSD (indep. X/Z)", "MWPM": "MWPM",
         "REF": "MLD / MPS reference"}


def md(df, fmt=".4g"):
    return df.to_markdown(index=False, floatfmt=fmt) if len(df) else "(empty)"


def family_summary(P):
    """Count Holm-significant wins/losses per comparison family."""
    rows = []
    for fam, g in P.groupby("family"):
        sig = g.holm_p < .05
        rows.append({"comparison": fam, "cells": len(g),
                     "a_better_sig": int((sig & (g.LER_a_minus_LER_b < 0)).sum()),
                     "a_worse_sig": int((sig & (g.LER_a_minus_LER_b > 0)).sum()),
                     "median_LER_diff_a_minus_b": g.LER_a_minus_LER_b.median()})
    return pd.DataFrame(rows)


def ler_wide(L, main):
    L = L.copy()
    L["m"] = np.where(L.method == L.reference, "REF", L.method)
    keep = ["MWPM", "BP_OSD", "CP_SA_local_SOE", "CP_local_SOE", main, "P4_2_Bennett", "REF"]
    w = L[L.m.isin(keep)].pivot_table(index=["d", "p"], columns="m", values="LER").reset_index()
    w["shots"] = L.groupby(["d", "p"]).shots.max().to_numpy()
    w["main_over_REF"] = w[main] / w["REF"]
    w["SA_over_REF"] = w["CP_SA_local_SOE"] / w["REF"]
    return w[["d", "p", "shots"] + keep + ["main_over_REF", "SA_over_REF"]]


def crossings(L, methods):
    """p where LER(d) and LER(d') curves cross (log-LER linear interpolation between grid p's)."""
    L = L.copy()
    L["m"] = np.where(L.method == L.reference, "REF", L.method)
    rows = []
    ds = sorted(L.d.unique())
    for m in methods:
        for d1, d2 in zip(ds[:-1], ds[1:]):
            a = L[(L.m == m) & (L.d == d1)].set_index("p").LER
            b = L[(L.m == m) & (L.d == d2)].set_index("p").LER
            ps = sorted(set(a.index) & set(b.index))
            diff = [np.log(b[p]) - np.log(a[p]) for p in ps]  # <0: larger d is better (below threshold)
            cross = None
            for (p0, f0), (p1, f1) in zip(zip(ps, diff), zip(ps[1:], diff[1:])):
                if f0 < 0 <= f1:
                    cross = p0 + (p1 - p0) * (-f0) / (f1 - f0)
                    break
            rows.append({"method": m, "d_pair": f"{d1}-{d2}", "p_range": f"{min(ps):.2f}-{max(ps):.2f}",
                         "crossing_p": cross if cross is not None else
                         ("none (larger d better everywhere)" if max(diff) < 0 else "none (larger d worse everywhere)")})
    return pd.DataFrame(rows)


def plot_gap(G, main, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps = [.06, .10, .14, .18]
    series = [(main, "main"), ("CP_local_SOE", "CP_local_SOE"), ("CP_SA_local_SOE", "CP_SA_local_SOE"),
              ("BP_OSD", "BP_OSD")]
    fig, axes = plt.subplots(1, len(ps), figsize=(12, 3.2), sharey=True)
    for ax, p in zip(axes, ps):
        for m, ck in series:
            s = G[(G.method == m) & np.isclose(G.p, p)].sort_values("d")
            if not len(s):
                continue
            ax.errorbar(s.d, s.gap_closed, yerr=[(s.gap_closed - s.ci95_lo).clip(lower=0), (s.ci95_hi - s.gap_closed).clip(lower=0)],
                        color=C[ck], marker="o", ms=4, lw=2, capsize=2, label=LABEL[m])
        ax.axhline(1, color="#9a9994", lw=1, ls="--"); ax.axhline(0, color="#9a9994", lw=1)
        ax.set_title(f"p = {p:.2f}", fontsize=10); ax.set_xlabel("code distance d")
        ax.grid(alpha=.25)
    axes[0].set_ylabel("MWPM-to-reference gap closed")
    axes[-1].legend(fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=160)


def plot_ratio(W, main, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3), sharey=True)
    for ax, (col, ck, name) in zip(axes, [("main_over_REF", "main", LABEL[main]),
                                         ("SA_over_REF", "CP_SA_local_SOE", LABEL["CP_SA_local_SOE"])]):
        for i, (p, s) in enumerate(W.groupby("p")):
            s = s.sort_values("d")
            shade = 0.25 + 0.75 * i / max(1, W.p.nunique() - 1)   # sequential: light -> dark = low -> high p
            ax.plot(s.d, s[col], marker="o", ms=4, lw=2, color=C[ck], alpha=shade, label=f"p={p:.2f}")
        ax.axhline(1, color="#9a9994", lw=1, ls="--")
        ax.set_title(name, fontsize=10); ax.set_xlabel("code distance d"); ax.grid(alpha=.25)
    axes[0].set_ylabel("LER / LER(reference)")
    axes[1].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(path); fig.savefig(str(path).replace(".pdf", ".png"), dpi=160)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="analysis")
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    v2, gr = root / "paper1_v2", root / "paper1_grid"
    sel = json.loads((v2 / "main_method.json").read_text(encoding="utf-8"))
    mainm = sel["main"]

    # ---- stage 1: 12,000-shot set
    L1 = pd.read_csv(v2 / "ler_table.csv")
    t1 = L1.pivot_table(index="method", columns=["d", "p"], values="LER")
    t1.columns = [f"d{d} p{p:.2f}" for d, p in t1.columns]
    t1 = t1.reset_index()
    P1 = pd.read_csv(v2 / "paired_comparisons.csv")
    P1["family"] = P1.a + " vs " + P1.b
    f1 = family_summary(P1)
    G1 = pd.read_csv(v2 / "gap_closure.csv")
    g1 = G1.pivot_table(index="method", columns=["d", "p"], values="gap_closed")
    g1.columns = [f"d{d} p{p:.2f}" for d, p in g1.columns]
    g1 = g1.reset_index()
    chi1 = pd.read_csv(v2 / "mps_chi_convergence.csv")
    E1 = pd.read_csv(v2 / "effective_distance.csv")

    # ---- stage 2: grid
    L = pd.read_csv(gr / "grid_ler.csv")
    P = pd.read_csv(gr / "grid_paired.csv")
    G = pd.read_csv(gr / "grid_gap_closure.csv")
    E = pd.read_csv(gr / "grid_effective_distance.csv")
    chi = pd.read_csv(gr / "grid_mps_chi_check.csv")
    W = ler_wide(L, mainm)
    fam = family_summary(P)
    gw = G.pivot_table(index=["d", "p"], columns="method", values="gap_closed").reset_index()
    cr = crossings(L, ["MWPM", "BP_OSD", "CP_SA_local_SOE", mainm, "REF"])
    sa_vs_sb = P[P.family == "CP_SA_local_SOE vs CP_local_SOE"][["d", "p", "LER_a_minus_LER_b", "a_right_b_wrong",
                                                                  "a_wrong_b_right", "holm_p"]]
    lowd = E[E.failures_below_ceil_d_half > 0]

    for name, t in {"stage1_ler": t1, "stage1_family_summary": f1, "stage1_gap_closure": g1,
                    "grid_ler_wide": W, "grid_family_summary": fam, "grid_gap_closure_wide": gw,
                    "grid_crossings": cr, "grid_SA_vs_SB_SOE": sa_vs_sb}.items():
        t.to_csv(out / f"{name}.csv", index=False)
    plot_gap(G, mainm, out / "fig_gap_closure_vs_d.pdf")
    plot_ratio(W, mainm, out / "fig_ler_ratio_to_reference.pdf")

    txt = [f"# Paper 1 analysis (code capacity, depolarizing, rotated surface code)", "",
           f"Main method (rule fixed before the grid): **{mainm}** "
           f"(r2 vs r1 on 12,000 shots: {sel['r2_wins']}/{sel['r2_losses']}, exact p={sel['exact_p']:.3g}).", "",
           "## Stage 1 — 12,000-shot set (d=5,7,9; p=.10,.12)", "", "### LER", "", md(t1), "",
           "### Paired comparisons (Holm within family over 6 cells)", "", md(f1), "",
           "### Fraction of the MWPM→MPS gap closed", "", md(g1, ".3f"), "",
           "### MPS chi=16 vs chi=32 class disagreements", "", md(chi1), "",
           f"Failures with error weight < ceil(d/2): {int(E1.failures_below_ceil_d_half.sum())} (all methods).", "",
           "## Stage 2 — p×d grid (held out; d=3..13, p=.05..18)", "",
           "### Paired comparisons (Holm within family over all cells)", "", md(fam), "",
           "### LER (REF = exact MLD at d=3, MPS chi=32 at d<=9, MPS chi=16 + chi=32 on hard shots at d>=11)", "",
           md(W, ".3e"), "",
           "### Fraction of the MWPM→reference gap closed", "", md(gw, ".3f"), "",
           "### Threshold crossings between consecutive distances", "", md(cr, ".3f"), "",
           "### CP-SA + SOE vs CP-SB + SOE (r1), per cell", "", md(sa_vs_sb, ".3g"), "",
           "### MPS chi check", "", md(chi), "",
           "### Effective distance: failures with error weight < ceil(d/2)", "",
           md(lowd) if len(lowd) else "none for any method", ""]
    (out / "ANALYSIS.md").write_text("\n".join(txt), encoding="utf-8")
    print(f"wrote {out / 'ANALYSIS.md'}")


if __name__ == "__main__":
    main()
