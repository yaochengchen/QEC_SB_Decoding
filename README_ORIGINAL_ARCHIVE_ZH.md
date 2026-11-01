# QEC_SB_paper1：第一篇论文的整理版（代码 + 汇总结果 + 分析 + 论文）

更新：2026-10-19。这里只放运行代码、分析/画图代码、能画出论文全部图表的汇总结果和论文本身。逐 shot 的原始数据（chunks/、hard/、DEM、npz 等）不在这里，完整结果你自己保存（`paper1/QEC_SB_paper1_results_logs.tar.gz`）。

这个文件夹用来归档和写论文，**不要解压覆盖到服务器上**。服务器上继续用 `QEC_SB_paper1_v3`；`code/` 与 `QEC_SB_paper1_v3.zip`（含 circuit_level）逐字节一致。

## 目录

```
QEC_SB_paper1/
├── code/        = QEC_SB_paper1_v3（qec_sb/、qec_sb/circuit_level/、scripts/、slurm/、run.sh）
├── results/     汇总结果
│   ├── paper1_v2/                   选择集 12,000 shot：主方法选择 main_method.json、MPS d=3 校验
│   ├── paper1_grid/                 留出网格 54 格：LER、配对检验、差距闭合、有效距离、χ 校验
│   ├── paper1_extra/                相关 BP+OSD、预算扫描（含 p=.18）、预算扩展行、χ=64 复核、种子复现、run_manifest_v4.json
│   ├── paper1_timing/               网格各解码器计时（p=.10）
│   ├── paper1_extra_timing_p0.10/   预算扫描各配置 + MPS χ=16/32 + PyMatching + 相关 BP+OSD 计时，100 syndrome
│   ├── paper1_extra_timing_p0.14/   同上，p=.14
│   └── circuit_locality/            电路级 DEM 稳定子群局域性：summary/（csv、md、tex 表）、每个实例的 T1_REPORT/T1_INSPECT、locality.log
├── results_soe_fast/   稀疏内核 + 快速 SOE：一致性检查 check_v2/（ALL_OK，tasks/ 未放）和计时 timing_v2_p0.10/（表 II、fig_budget 用）
├── analysis/
│   ├── make_paper_figures.py        从 results/ 生成论文的 4 张 PDF 图
│   └── analyze_paper1_summary.py    选择集 + 网格的汇总分析
├── paper/
│   ├── main.tex, refs.bib           论文（IEEEtran，投 IEEE TQE）
│   ├── figures/fallback/            PDF 图生成前的占位 PNG
│   └── main_preview.pdf             通读用预览（用 article 替代类编译，正式排版以 Overleaf 为准）
└── notes/FINDINGS_ZH.md             当前结论汇总
```

## 生成论文图

在本文件夹根目录运行（几秒钟）：

```bash
python analysis/make_paper_figures.py --results results --out paper/figures
```

若 `results_soe_fast/timing_v2_p*/` 存在，`fig_budget` 的 CP-SB/CP-SA 搜索、SOE 和端到端时间自动改用新计时（`--soe-fast` 可指定别的路径）。

生成 `fig_ratio_vs_p.pdf`、`fig_ler_vs_d.pdf`、`fig_decomposition.pdf`、`fig_budget.pdf`。`main.tex` 会优先用 `paper/figures/*.pdf`，没有时退回 `figures/fallback/` 的 PNG。脚本只做过语法检查，某张图报错时会打印 FAIL 并继续生成其他图。

## 编译论文

Overleaf 直接编译，或本地 `pdflatex main && bibtex main && pdflatex main && pdflatex main`。

## 还需要你填的（正文红色 TODO）

- JiaPei 的姓和单位；是否加其他作者
- FAPESP 基金号
- 代码/数据仓库链接
- `refs.bib` 里标了 `VERIFY` 的 5 条 2024–2026 年 arXiv 文献，需要对照 arXiv 核对作者和标题

## 本地可以删掉的旧东西（我没有删任何文件）

- `QEC_SB_paper1/results/paper1_extra_timing/`：旧的 30-syndrome 计时，已被 `_p0.10`、`_p0.14` 取代（画图脚本已不再读它）
- `paper1/QEC_SB_paper1_v3.zip`（10-17 旧版）、`paper1_v2/`、`QEC_SB_paper1_v2/`、`QEC_SB_paper1_v2.zip`
- `paper1/analysis_2026-10-18/`、`paper1/paper1_extra/`、`paper1/results/`：内容都已并入本文件夹
