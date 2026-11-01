# QEC-SB 第一篇补实验（v3，含补跑）

## 新增（2026-10-19，v2）：稀疏搜索内核 + SOE r=2 的 O(m) 精确算法 + 工作量计数（只重测计时，不重跑任何解码）

```bash
bash slurm/submit_soe_fast.sh     # 两个独占节点作业并行：一致性检查（约 20-40 分钟）+ 计时 v2（单进程，约 2-3 小时）
```

如果上一版（`results_soe_fast/check/`、`timing_p0.10/`）的作业还在排队或运行，直接 `scancel`，已被 v2 取代。
结果全部写到**新文件夹 `results_soe_fast/`**，不动 `results/`：

- `results_soe_fast/check_v2/SOE_FAST_CHECK.md`：第一行应为 **ALL_OK**。同一批 shot 上检查两件事：
  (1) 稀疏与原稠密 ΔE 内核得到的候选池逐位相同（CP-SB 8×100、d≥11 另加 8×400、CP-SA 8×100）；
  (2) 快速 SOE 与原枚举的每类得分、球大小、最小权重、判决相同（判决不同只允许出现在 Z 严格相等的"精确平局"上）。
  d=3–13，p=.06/.10/.14/.18。
- `results_soe_fast/timing_v2_p0.10/TIMING_V2.md`、`timing_budget.csv`、`work_counts.csv`：与原计时相同的 syndrome、
  solver seed；13 个预算配置全部用稀疏内核 + 快速 SOE 重新计时；sbR8S100、saC8W100 同时计原稠密搜索和原 SOE 枚举
  （前后对比、加速倍数）；与实现无关的工作量计数（每次解码的单翻转 ΔE 评估次数、串行轮数）；CPU 型号。

跑完把整个 `results_soe_fast/` 发回（`check_v2/tasks/` 可不带）。

改动（都不改变任何判决，已有结果全部有效）：

- `qec_sb/cp_sb.py`、`qec_sb/cp_sa.py`：单个生成元翻转的 ΔE 只读该生成元碰到的 ≤4 个比特（默认 `sparse`），
  原来对整条 2n 位错误向量重算（`QECSB_DELTA=dense` 可切回）。整数完全相同、随机数调用顺序不变 → 轨迹逐位相同。
- `qec_sb/fast_soe.py`：SOE 半径 1/2 由权重直方图精确算出，O(m)，与 `fixed_soe_prediction` 同接口。
- `scripts/run_soe_fast.py`（check / timing / counts）、`slurm/60_soe_fast_check.sbatch`、`slurm/61_soe_fast_timing.sbatch`、
  `slurm/submit_soe_fast.sh`。

工作量计数（每次解码，m=d²−1）：SB R×S：ΔE 评估 8·R·S·m，串行轮数 S；SA C×W：ΔE 评估 8·C·W·m，串行轮数 W·m。

## 新增（2026-10-19）：电路级 DEM 稳定子群的局部性（附录表），只在 n013 上跑

```bash
bash slurm/submit_locality.sh            # 一个作业，n013，512 CPU；可续跑（重新提交即可）
WBIG=64 bash slurm/submit_locality.sh    # d>=11 和 bb72 用的 worker 数（默认 96），内存不够时调小
```

只用 CPU，不解码。实例：表面码 d=3,5,7,9,11,13（rounds = d，Stim `rotated_memory_z`，四项噪声都 = p）和
[[72,12,6]] BB 码（Bravyi depth-8 电路，6 轮），p=0.003。每个实例：建 DEM → 自解码生成局部生成元 → 分析全局元素；
每做完一个实例就重写一次汇总。预计：d≤9 几分钟，d=11 约 10 分钟，d=13 约 20–60 分钟，bb72 约 10 分钟。
进度 `results/circuit_locality/locality.log`。跑完发回 `results/circuit_locality/summary/`（3 个小文件）和
`results/circuit_locality/*/T1_INSPECT.md`。
代码：`qec_sb/circuit_level/`（gf2, dem, kernel_basis, codes, circuits；来自第二篇 Phase 0 包，ldpc 导入改为不依赖 sinter）、
`scripts/run_circuit_locality.py`、`slurm/50_circuit_locality.sbatch`、`slurm/submit_locality.sh`。

全部是 CPU 程序，只做 code-capacity depolarizing 噪声。分两个阶段，顺序不能换：

1. **确认集（12,000 shot，原 seed）**：复现原表，在同一批 shot 上补 PyMatching、BP+OSD、CP-SA、MPS（χ=16/32）和聚合半径 2 消融。分析结束时按**事先写死的规则**选出主方法，写进 `results/paper1_v2/main_method.json`：默认是 CP-SB + local SOE（半径 1）；只有当半径 2 在全部 12,000 shot 上失败更少、且配对精确检验 p<0.05 时才改用半径 2。importance 和 Bennett 只作消融，不参与选择。
2. **p×d 网格（全新 seed）**：主方法已固定，网格是独立的留出检验。d=3…13；p=.05,.06,.08,.10,.12,.14,.16,.18（d≥11 只跑 p≥.08）。

## 补跑（网格、补充实验、过夜延伸都跑完以后）

```bash
QECSB_LOWP_EXT=1 python scripts/run_paper1_extra.py estimate --out results/paper1_extra   # 先看剩余工作量
bash slurm/submit_followup.sh
```

一条命令提交三步：补充实验数组（所有可用节点）→ 两个计时作业（p=.10、p=.14，各独占一个节点，并行）→ 补充分析。
已有的块全部跳过，中断后重新执行同一命令即可续跑。新增内容：

- **预算扫描延到 p=.18**（d=9/11/13，前 1 万 shot）：只跑 `sbR8S100`、`sbR8S400`、`sbR32S400`、`saC8W100`、`saC32W200`，
  用来实际测出预算随 d 扩展后 d=11→13 的交叉点（之前只能说"≥.16"）。
- **d=11/13 的预算扩展行**：按事先声明的规则 CP-SB 8×400（d≥11 步数 400，副本不变）在网格格点本身上解码，
  输出 MAP / SOE r1 / SOE r2，与冻结主方法、CP-SA、参照逐 shot 配对（`budget_scaled_ler.csv`、`budget_scaled_paired.csv`）。
  p≥.08 用全部网格 shot；低 p 延伸格（p≤.06）默认只取前 10 万 shot（环境变量 `QECSB_FULL_LOWP_CAP`）。
  扫描里已有的 8×400 文件（p=.08–.16 前 1 万 shot）原样复用。规则写在 `run_manifest_v4.json`。
- **SB 对 SA 同半径配对检验** `budget_paired_matched.csv`：例如 8×400 对 8×100、32×400 对 32×200，
  MAP/r1/r2 各自比较，另含"步数对副本"的 8×400 对 32×100 / 64×100。
- **计时**改为每个 p 一个作业、默认 100 个 syndrome，同时计 MPS χ=16/32 和 PyMatching，
  结果在 `results/paper1_extra_timing_p0.10/`、`..._p0.14/`；分析会把所有计时目录合并。
- 三张新图：`fig_budget_accuracy_vs_cost`（横轴改为端到端时间，SB/SA 各画 r1、r2，并标出 MPS χ=16 和相关 BP+OSD）、
  `fig_ler_vs_d_budget_scaled`（LER 随 d：冻结 / 预算扩展 / 参照）、`fig_decomposition`（MAP → r1 → r2 的差距分解）。
- 工作量（空闲节点时间估计）：预算扩展行约 740 core-h，p=.18 延伸约 160 core-h，计时每个 p 约 4–5 小时（单核）。

### 节点选择（所有 submit 脚本通用，见 `slurm/common.sh`）

- 默认**不使用 n013**（`#SBATCH --exclude=n013` 已写进每个 sbatch 文件；`EXCLUDE=n013,n014` 可多排，`EXCLUDE=` 置空则不排除）。
- 数组大小 `N` 默认 = 分区里除排除节点外的全部可用节点数（`sinfo` 查到的 idle/mixed/allocated，不含 down/drain），
  上限 `MAXN=32`；也可以直接 `N=6 bash slurm/submit_followup.sh`。每个数组任务 `--exclusive` 独占一个节点。
- 老脚本（`submit_all.sh`、`submit_extra.sh`、`submit_overnight.sh`、`resubmit_failed.sh`）同样改用这套设置，
  原来的 `GRID_JOBS`/`JOBS` 仍可覆盖。

## 原说明

## 运行

```bash
conda activate <你的环境>     # Python ≥ 3.10
bash run.sh 1024             # 1024 = 并行进程数，每个进程约 150 MB 内存
```

`run.sh` 依次做：装依赖 → 冒烟测试 → 阶段 1 → 复现闸门（不通过就停）→ 阶段 2。中断后重跑同一命令会续跑。

先估算网格耗时：`python scripts/run_paper1_grid.py estimate --workers 1024`

**计时**：服务器满载时测的时间不能用。全部跑完后，在一台没有其他任务的机器上单独跑：

```bash
python scripts/run_paper1_grid.py timing --out results/paper1_timing
```

跑完把 `results/paper1_v2/`、`results/paper1_grid/`、`results/paper1_timing/` 打包发回来（`chunks/` 很大的话可以不带，汇总 csv 和 md 足够）。

## 跑完先看

1. `results/paper1_v2/reproduction_check.csv`：`match` 和 `r2_pools_rebuilt_identically` 全为 True（run.sh 会自动检查）。
2. `results/paper1_v2/MPS_D3_CHECK.md`：非平局 shot 上 MPS 与穷举 MLD 100% 一致。
3. `results/paper1_v2/main_method.json`：选出的主方法。
4. `results/paper1_grid/grid_mps_chi_check.csv`：χ=16 与 32 的分歧（d≥11 只在"难 shot"上检查：MPS χ=16 判错或任何解码器与它不一致的 shot）。
5. `grid_effective_distance.csv`：`failures_below_ceil_d_half` 应为 0；不为 0 说明解码器丢了有效距离。
6. 主要结果：`RESULTS_V2.md`、`RESULTS_GRID.md`、`grid_ler_vs_p.pdf`、`grid_ler_vs_d.pdf`。

## 和上一版（v2）相比的改动

- 网格加入半径 2、p=.16/.18（用于阈值交叉）；d≥11 去掉 p=.05/.06（几乎没有失败，没有信息量）。
- 网格新增主方法对 MWPM / PyMatching / BP+OSD / 参照的配对检验，以及 CP-SB 对 CP-SA；Holm 按比较族、跨所有格点校正。
- d≥11 的参照改为：χ=16，难 shot 上用 χ=32 覆盖。
- 每个 shot 记录物理错误权重，输出有效距离检查表。
- 计时单独成 `timing` 子命令：单进程、只计非平凡 unique syndrome，报中位数和 p90。
- 新增 LER 随 d 变化图。
- ldpc 安装失败（sinter 编译问题）时自动退回不带依赖安装。

## 分步命令

```bash
python scripts/run_paper1_v2.py {check|mps-d3|ours|baselines|analyze|all} --out results/paper1_v2 --workers N
python scripts/run_paper1_grid.py {estimate|run|hardchi|analyze|timing|all} --out results/paper1_grid --workers N
```

## 在 Slurm 集群上运行（sbatch）

```bash
# 1) 登录节点，建 conda 环境（一次）
bash slurm/setup_env.sh
# 2) 冒烟测试 + 耗时估算（约 30 分钟），看 logs/qecsb_smoke_*.log 末尾是否有 SMOKE_OK
mkdir -p logs && sbatch slurm/00_smoke.sbatch
# 3) 一次提交全部正式作业（按依赖关系自动排队）
bash slurm/submit_all.sh          # 默认 4 节点 × 512 核
```

默认 CPUS=512、GRID_JOBS=4（每个数组任务独占一个节点，--exclusive --mem=0）。
需要指定分区/账号时加 `EXTRA="--partition=xxx --account=yyy"`。
某一步超时或被杀：`STEP=grid bash slurm/resubmit_failed.sh`，已完成的块会跳过。
复现闸门不通过时，后面的作业会一直挂着（DependencyNeverSatisfied），用 `scancel` 取消。

## 补充实验（网格跑完以后）

```bash
bash slurm/submit_extra.sh
```

在网格已有的 shot 上逐 shot 配对加跑两项，结果写到 `results/paper1_extra/`：

- **相关 BP+OSD**（全部 44 个格点）：在 [Hz 0 Hz; 0 Hx Hx] 上解 X/Z/Y 三列，每列先验 p/3，替换原来 X/Z 分开解的弱基线。
- **搜索预算扫描**（d=9/11/13，p=.08–.16，每格前 1 万 shot）：
  - CP-SB 的 (replicas, steps)：8×100（冻结值）、16×100、32×100、64×100、8×200、8×400、32×400
  - CP-SA 的 (每种模式链数, sweeps)：8×25、8×50、8×100（网格值）、8×200、16×100、32×200
  - 每个配置都算 MAP、SOE 半径 1、SOE 半径 2。

8×100 这两个基准配置必须和网格结果逐 shot 一致（见 `budget_consistency.csv`）。
解码结束后，会在独占节点上为每个配置单独计时（`results/paper1_extra_timing/`），
图 `fig_budget_accuracy_vs_cost` 的横轴用的就是这次计时。

## 过夜补充（`AFTER=<补充实验作业号> bash slurm/submit_overnight.sh`）

用环境变量 `QECSB_LOWP_EXT=1` 打开低 p 延伸。原有格点的 seed 和结果完全不变，新格点使用独立的 seed 命名空间。

- 新格点：p=.04 覆盖全部 d（d≤9 用 40 万 shot，d≥11 用 100 万 shot），另加 d=11/13 的 p=.05/.06（各 40 万 shot），所有方法逐 shot 配对。
- 新格点上的相关 BP+OSD。
- 基准 CP-SB 和 CP-SA 各做两组独立 solver seed 的重复（d=9/11/13，p=.10/.14），结果在 `seed_replicates.csv`。
- d≥11 的难 shot 上，用 MPS χ=64 复核参照，结果在 `mps_chi64_hard_check.csv`。
- 最后把网格分析和补充分析整体重做一遍，覆盖原有结果。
