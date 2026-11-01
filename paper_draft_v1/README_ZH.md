# 当前正式稿

已按用户批准合并全文真实性审计修订及Abstract、Introduction、Methods的26项语言润色，以及其余章节与补充材料的51项语言润色。图4为优化实现的准确率—端到端时间图；正文保留附录A，资源图S1和配对图S2及五个扩展切片位于补充材料。附录B/C的既有安排不变，参考文献32条。

- [英文正式稿](MANUSCRIPT.md)
- [Physical Review Applied双栏PDF](prapplied_release/manuscript.pdf)
- [单栏校对PDF](latex_release/manuscript.pdf)
- [补充材料PDF](SUPPLEMENTAL_MATERIAL.pdf)；[补充材料中译](SUPPLEMENTAL_MATERIAL_ZH.md)
- [独立补充源码与PDF包](SUPPLEMENT_SOURCE_AND_PDF.zip)；[构建脚本](supplement_release/build.py)
- [完整投稿准备包](SUBMISSION_PREP_BUNDLE.zip)
- [摘要中译](ABSTRACT_ZH.md)、[引言](INTRODUCTION_ZH.md)、[方法](SECTION_II_ZH.md)、[评估](SECTION_III_ZH.md)、[结果](SECTION_IV_ZH.md)、[讨论](SECTION_V_ZH.md)、[结论](SECTION_VI_ZH.md)、[附录A](APPENDIX_A_ZH.md)
- [文献管理](REFERENCES_MANAGEMENT_ZH.md)；[数据索引](ARCHIVE_DATA_INDEX.md)

## 本次采用的修订

修正补充PDF旧图4描述；SOE新旧时间统一采用优化批次，旧值为3.53979秒（表2为3.540）；明确42,800次实现检查覆盖范围；补齐独立补充包图源和拼接脚本；自由能术语改为effective free energy；精简重复计时说明。公式、图片、参考文献及其他科学结果不变。

完整合并前备份：before_approved_factual_audit/。获批候选及红绿对比保留于factual_audit_candidate/，其中“待批准”措辞为当时状态记录，当前已采用。旧图4仍保留在figures/preview_fig4.png，不用于正文或新增补充页。

补充材料可独立构建；数据索引仍指向项目文件，公开归档URL/DOI尚待补充。独立构建需要Python及pymupdf、Pandoc、XeLaTeX和Times New Roman字体。


26项语言修订合并前备份：before_approved_language_26/；完整红绿对比及理由保留于language_review_abs_intro_methods/。后续51项修订也已获用户全部批准并合并。

## 51项语言修订发布记录

合并前备份：`before_approved_language_51/`。红绿对比及逐项理由完整保留于 [正文对比](language_review_remaining/manuscript.pdf)、[补充材料对比](language_review_remaining/SUPPLEMENTAL_MATERIAL.pdf) 和 [修改理由](language_review_remaining/EDIT_REASONS_ZH.md)。这些目录中的待审标记记录历史状态，当前51项已全部采用。

已同步中文译文、双栏和单栏LaTeX/PDF、补充材料与发布包。补充材料S3引导文字现在纳入PDF，仍为7页。见 [最终一致性检查](LANGUAGE_51_MERGE_CHECK_ZH.md)。

## 2026-11-01 引文删除

按作者要求删除原文献[28]及引言tabu介绍；原[29–33]顺延为[28–32]。中译、LaTeX及发布包同步。备份：`before_remove_reference28_20261101/`。其他公式、数据和图片不变；历史审阅报告使用旧编号。
