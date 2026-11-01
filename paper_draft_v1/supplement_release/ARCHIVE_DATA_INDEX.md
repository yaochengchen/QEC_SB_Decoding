# Archived data and figure index

Repository-local source index for the current manuscript and Supplemental Material. These paths do not by themselves establish a public release or stable external identifier.

The [budget-figure manifest](../analysis_outputs/evidence_budget/manifest.json) records input hashes and analysis conventions. [Resource plot data](../analysis_outputs/evidence_budget/resource_plot_data.csv), [paired plot data](../analysis_outputs/evidence_budget/paired_plot_data.csv), and [plotting script](../analysis/plot_evidence_budget.py) support the resource and paired-effect figures. The manifest does not capture a complete software environment or output hashes.

| Evidence | Archived file |
|---|---|
| Main-decoder LER and Wilson intervals | [results/paper1_grid/grid_ler.csv](../results/paper1_grid/grid_ler.csv) |
| Main-decoder paired comparisons | [results/paper1_grid/grid_paired.csv](../results/paper1_grid/grid_paired.csv) |
| Correlated BP+OSD comparisons | [results/paper1_extra/bposd_corr_paired.csv](../results/paper1_extra/bposd_corr_paired.csv) |
| Budget LER, all inference rules | [results/paper1_extra/budget_ler.csv](../results/paper1_extra/budget_ler.csv) |
| Selected configuration-pair effects | [results/paper1_extra/budget_paired_matched.csv](../results/paper1_extra/budget_paired_matched.csv) |
| Extended-budget LER | [results/paper1_extra/budget_scaled_ler.csv](../results/paper1_extra/budget_scaled_ler.csv) |
| Extended-budget paired tests | [results/paper1_extra/budget_scaled_paired.csv](../results/paper1_extra/budget_scaled_paired.csv) |
| Historical baseline timings used in Fig. 4; historical CP timings retained in the project archive | [results/paper1_extra_timing_p0.10/timing_budget.csv](../results/paper1_extra_timing_p0.10/timing_budget.csv) |
| Optimized implementation timings used in Table 2 | [results_soe_fast/timing_v2_p0.10/timing_budget.csv](../results_soe_fast/timing_v2_p0.10/timing_budget.csv) |
| Defined search operation counts | [results_soe_fast/timing_v2_p0.10/work_counts.csv](../results_soe_fast/timing_v2_p0.10/work_counts.csv) |

Additional resource slices: [p = 0.08](../analysis_outputs/evidence_budget/resources_p0.08.pdf), [0.10](../analysis_outputs/evidence_budget/resources_p0.10.pdf), [0.14](../analysis_outputs/evidence_budget/resources_p0.14.pdf), [0.16](../analysis_outputs/evidence_budget/resources_p0.16.pdf), and [0.18](../analysis_outputs/evidence_budget/resources_p0.18.pdf). The [archived LER-versus-distance view](../paper/figures/fig_ler_vs_d.pdf) is additional context. All views reuse existing observations.
