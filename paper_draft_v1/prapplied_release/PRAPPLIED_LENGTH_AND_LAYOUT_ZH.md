# Physical Review Applied 篇幅调研与本稿双栏排版

核查日期：2026-10-31。

## 对本稿的判断

建议按 **Regular Article** 准备。本次在不压缩、删减或改写已定稿科学内容的前提下，用 REVTeX 4.2 的 `aps,prapplied,reprint` 选项排成双栏，得到 **15页**。全文仍包括摘要、六节正文、三个附录、七张图和35条参考文献。此前的19页是单栏阅读版，两种版式的页数不能直接当作内容多少或期刊限制来比较。

当前正文及结论延伸至第10页；附录从第10页开始，参考文献位于第14–15页。起止页有共享，因此不能将这些页数简单相加。图片较多且采用跨栏显示，是总页数的重要组成部分。本稿15页不构成 Regular Article 的形式性超长，是否还应精简应由信息必要性和可读性决定。

## 1. 官方要求：普通研究论文没有硬性长度上限

Physical Review Applied 的作者说明明确区分：Regular Articles 无长度上限；Letters 限4,500词；Review Articles 限30,000词。无上限并非鼓励冗长，普通论文仍应与所呈现信息相称。本文属于原创方法与评估研究，适合 Regular Article；不能把 Letter 的限制套用到它。见[期刊作者说明](https://journals.aps.org/prapplied/authors)与[编辑政策](https://journals.aps.org/authors/editorial-policies-practices-physical-review-applied)。

APS 的篇幅规则采用正文、展示公式、图片和表格的等效词数，而非单纯数 PDF 页数。对于 Letters/短文，摘要、作者信息和参考文献等不计入所述正文词数，但图表及公式有对应的计数规则。因此此前用 Markdown 空白分词得到的压缩比例，只能用于版本比较，不能当作期刊正式字数。见[APS 篇幅计算指南](https://journals.aps.org/authors/length-guide)。

## 2. 已核实的相关论文长度

下表统计的是正式发表排版 PDF 的总页数，包含文件内的附录与参考文献，不包含另行下载的 Supplemental Material。优先选择2025–2026年的模拟分岔、物理启发式优化与解码论文，补充一篇2024年紧密相关的实验解码工作。不是随机抽样，也没有把 arXiv 预印本页数混入正式版。

| 论文 | 年份与卷号 | 正式PDF总页数 | 与本稿的关联 |
|---|---|---:|---|
| [Improved update rule for probabilistic computers](https://journals.aps.org/prapplied/abstract/10.1103/t1wc-j7yp) | 2025，24，054009 | 11 | Physics-Inspired Computing 专题中的算法与FPGA资源比较；页数来自[作者实验室提供的正式发表PDF](https://bpb-us-w2.wpmucdn.com/u.osu.edu/dist/7/38882/files/2025/11/PhysRevAppl24_054009_2025.pdf) |
| [Edge-of-chaos-enhanced quantum-inspired algorithm for combinatorial optimization](https://journals.aps.org/prapplied/pdf/10.1103/2qd9-x6v8) | 2026，25，044011 | 13 | 直接相关的模拟分岔改进与性能评估 |
| [Reducing the error rate of a superconducting logical qubit using analog readout information](https://journals.aps.org/prapplied/pdf/10.1103/PhysRevApplied.22.044031) | 2024，22，044031 | 19 | 量子纠错解码、实验数据与性能分析 |
| [Practical hybrid decoding scheme for parity-encoded spin systems](https://journals.aps.org/prapplied/pdf/10.1103/drry-blqs) | 2026，25，014046 | 20 | 编码自旋系统、混合解码与仿真；注意不是表面码QEC解码论文 |
| [Enhancing decoding performance using efficient error learning](https://journals.aps.org/prapplied/pdf/10.1103/vtj7-nh6l) | 2026，25，054065 | 22 | 量子纠错中的噪声学习与最大似然解码性能 |

这五篇的范围为11–22页，中位数为19页。这个中位数只是**所选相关样本**的描述，不能代表全刊，也不能用作建议的目标页数。样本偏向方法、推导及评估较充分的研究，且受到公开可获取PDF的影响。两个优化算法例子为11和13页，三个解码相关例子为19、20和22页，显示相关研究的篇幅跨度较大。

因此，不能凭此次调研断言“PRApplied 通常只有8–10页”或“全刊典型长度为19页”。更稳妥的实际判断是：本稿15页位于已核实相关论文的范围内，不需要为了一个不存在的普通论文页数上限再次删减科学内容。若以后主动追求更短，可另行评估附录B的文件索引及附录C是否移至补充材料；本次没有进行这种迁移。

## 3. 双栏版本采取的排版措施

- 使用 APS 官方推荐的 REVTeX 4.2，期刊选项为 `prapplied`，`reprint` 生成双栏预览；这属于作者排版稿，并非出版社最终校样。官方格式入口见[作者说明](https://journals.aps.org/prapplied/authors)和[REVTeX主页](https://journals.aps.org/revtex)。
- 标题、作者、单位、通讯邮箱和摘要使用 REVTeX 原生结构。
- 主文节、小节及附录使用模板的层次样式；保留既有正文与附录编号含义。
- 图1–5及S1–S2跨双栏，以保持多面板坐标和图例可读；表1、表C1和无编号复现文件索引跨栏，表2单栏。没有改变原图像像素或图中数据。
- 长公式使用多行排版，保留所有数学项、关系、符号和(1)–(10)、(A1)、(A2)编号。
- 参考文献改用原生 bibliography/citation 结构；35条著录内容及编号顺序不变。为保持现有编号，表号继续采用1、2、C1，没有擅自改成罗马数字。

## 4. 核对结果与限制

当前正式Markdown与最新图1是输入。排版过程中发现图1图注和对应中文译文被工作区中的其他操作更新，因此没有继续使用过时的单栏LaTeX作为内容来源；双栏构建先从当前Markdown重新转换，再施加格式变换。本任务没有回写正式正文、中文译文或旧单栏PDF。

从当前输入生成的240个行内/展示数学项，在仅忽略换行、对齐和间距命令后逐项一致；10条图表标题及说明完整保留；七幅图片哈希一致；35条参考文献及正文引用序列一致；公式和内部交叉引用目标有效。

全部15页已渲染并检查，重点复核长公式、表C1、多面板图与参考文献。编译无越栏、缺字或未定义引用错误。日志仍有双栏两端对齐造成的 underfull 提示，以及 REVTeX/hyperref 的标签兼容性提示；未出现文字截断或遮挡，不将这些提示隐瞒为“零警告”。详细机器核查见本目录 verification.json。

本文科学内容未修改；本次也没有新增实验、重新做文献新颖性审计或替作者补写投稿声明。项目内的复现文件链接仍是项目相对路径，独立PDF不能替代完整数据归档。
