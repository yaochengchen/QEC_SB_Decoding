# 标题与 References 正式格式整理记录

日期：2026-10-31。

## 已完成

1. 采用用户确认的 A 方案：**Quantum-Inspired Search Meets Quantum Error Decoding: Constraint-Preserving Simulated Bifurcation with Local Degeneracy Aggregation**。
2. References 全部改为英文正式著录，去除“出版社”“论文全文存档”“作者引用元数据”等链接说明、中文标点和重复来源链接。
3. 按 APS/Physical Review 风格使用作者姓名首字母、期刊缩写、粗体卷号、起始页或文章号、年份；PDF 文献标签为 [1]–[35]。DOI 或公开来源链接保留在题名上，便于直接跳转。
4. [15]、[17]、[21] 补全不超过十人的作者名单；较长作者组使用首位作者 et al.。保留 Google Quantum AI and Collaborators 的正式团体署名。
5. [11] 由 arXiv 预印本更新为同一研究的正式期刊版本：F. Martínez-García, F. R. F. Pereira, and P. Parrado-Rodríguez, Entropy 28, 91 (2026), DOI 10.3390/e28010091。[27] 题名连字符按正式版本校正。
6. [30]、[31] 根据作者仓库的 CITATION.cff 保留为 unpublished manuscript (2026)，不把“已投稿”“审稿中”写成“已发表”或“已接受”。其余尚无核实期刊信息的条目保留 arXiv 编号。
7. 同步英文 Markdown、参考文献管理中的正式清单及当前位置索引、LaTeX、PDF 和源码包。中文管理文档中的介绍、引用理由和历史核验说明保留。

## 核对结果

- 35 条编号与正文引用对应关系不变；正文科学内容（排除首行标题与 References）与修改前逐字一致。
- 10 个正式中文译文的版本指纹与结构同步检查通过；正文和中文译文无需因本轮著录格式而重写。
- 269 处数学表达式、公式 (1)–(10)/(A1)/(A2) 和七张图片均通过核对。
- 新 PDF 为 21 页，标题页和两页参考文献已重新渲染检查；不存在中文著录残留、重复逗号、缺字、溢出或未解析引用警告。
- 一致性检查及其 6 个回归测试通过；原有 paper/ 下登记的 11 个文件保持不变。
- 修改前初稿、PDF、LaTeX 和相关文件保存在 `../before_reference_format/`。

## 格式与著录来源

- [APS Style Basics](https://journals.aps.org/authors/style-basics)：作者名单、数字编号和已发表版本优先等规则。
- [APS Journals Style Guide](https://res.cloudinary.com/apsphysics/image/upload/v1736779890/APS_Journals_Style_Guide_Authors_Nov2024_ua91lv.pdf)：期刊缩写、卷号、起始页等著录格式。
- [第 11 条正式出版页面](https://www.mdpi.com/1099-4300/28/1/91)，其出版信息同时由 Crossref DOI 元数据核对。
- [clock-dSB 作者引用信息](https://github.com/yaochengchen/clock-dsb/blob/main/CITATION.cff)和 [PdSB 作者引用信息](https://github.com/yaochengchen/PdSB/blob/main/CITATION.cff)：未发表稿件状态及作者顺序。

Crossref 核对结果保存在 `reference_format_audit/crossref.json`。期刊年份按卷期版本保留，不用较早的在线发表日期覆盖。该轮是参考文献著录整理和相关信息核验，不是重新评估全部文献的科学结论，也不是宣称整篇稿件已经满足特定期刊的全部投稿要求。
