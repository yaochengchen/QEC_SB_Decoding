# 51项语言修订合并与一致性检查

## 合并范围与批准依据

依据用户“这51项修改全部采用”，将 language_review_remaining/edits.json 中 III-01—III-10、IV-01—IV-12、V-01—V-09、VI-01—VI-02、A-01—A-10、SM-01—SM-08 全部合并。正式英文正文和补充材料与已批准 CLEAN 文件逐字一致。第一轮26项保持不变。

合并前版本保存在 before_approved_language_51/；两轮红绿对比目录均保留。第二轮对比仍可打开 [正文](language_review_remaining/manuscript.pdf)、[补充材料](language_review_remaining/SUPPLEMENTAL_MATERIAL.pdf) 和 [逐项理由](language_review_remaining/EDIT_REASONS_ZH.md)。

## 同步与检查结果

- 已逐项同步 III—VI、附录A及补充材料中文译文。SM-07修复英文双谓语，原中文没有这一语病，保持等义；SM-08同步明确切片指代。
- 正式稿与批准候选完全相同；合并前后全部行内/独立公式、数值序列、图片路径保持不变。摘要、引言及方法全文未改变。
- 8个正式中译章节通过公式、图注编号、图路径、引用与版本指纹检查；文献管理的当前引用位置索引一致，33条文献不变。对修订句另作人工语义对照，指纹检查本身不证明翻译等义。
- 数据检查核对54个主评估单元、7,120,000 shots及现有配对比较汇总；原始论文11个受保护文件保持不变。此次不重新运行解码实验或外部文献审计。
- 单栏源码的207项数学表达与Markdown对应；公式编号1—10、A1—A2，5张正文图、33个文献目标和交叉链接通过检查。全部发布图文件与备份逐字节一致。
- 双栏正文12页、单栏16页、补充材料7页。三个构建日志均无Overfull、未定义引用或缺字错误。单栏保留1条Underfull行宽提示，不影响内容。
- 已检查双栏全部页面缩略图及补充材料页面，另放大检查第2页S2/S3衔接；图表无可见截断。PDF文本边界检查通过。
- 原补充构建脚本省略S3引导文字，现纳入已批准文字，位于第2页；后续五张切片及图号不变。
- 双栏包、单栏包、独立补充包与完整投稿准备包已重建。独立补充包解压到临时目录后成功重建7页PDF；完整包中的正文及补充PDF与正式发布PDF逐字节一致。

## 发布入口

- [双栏正式PDF](prapplied_release/manuscript.pdf)
- [单栏校对PDF](latex_release/manuscript.pdf)
- [补充材料PDF](SUPPLEMENTAL_MATERIAL.pdf)
- [完整发布包](SUBMISSION_PREP_BUNDLE.zip)
- [独立补充包](SUPPLEMENT_SOURCE_AND_PDF.zip)

公开数据归档URL/DOI仍待补充，这是既有发布事项。此次未增加科学内容、修改数据或更换图源。
