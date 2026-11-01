# 优化图4正式合并记录

用户批准采用优化实现图作为正文图4，原历史图只留项目归档。现已合并，不改变S1/S2和附录安排。

左面板展示全部归档预算的半径二错误率与端到端时间，固定χ=16 MPS使用自身准确率记录。右面板改为逐配置连接半径一、二的错误率点估计，并单列端到端中位数之差，避免微小时间差在对数坐标中重叠。右面板不提供配对显著性或计时置信区间；原候选图和数据继续保留。

同步范围：英文图注与相关正文、III/IV中文译文、LaTeX、单双栏PDF、补充材料中引用正文图4的一句话、README、索引和发布ZIP。没有改动其他附录结构、公式和参考文献编号，也没有新增实验。不同准确率样本量、非平凡综合征计时及一次性准备排除条件写入图注。

旧图仍在figures/preview_fig4.png；完整合并前备份在before_approved_fig4_optimized/。新绘图脚本为analysis/plot_figure4_optimized_final.py，数值导出在analysis_outputs/figure4_optimized_final/plot_data.csv。

全文verify_draft.py通过；公式与参考文献保持不变。PDF编译无未定义引用或横向溢出；页码和页面预览见fig4_optimized_merge_audit/。结果只支持具体实现下的准确率—时间取舍，不建立与MPS准确率等效，也不是硬件加速证明。
