__version__ = "2.10.6"
# 2023.08.11. v1.1.0	adding a option for not normalizing the bw file
# 2023.10.26. v1.9.0	prerelease for v2.0
# 2023.10.27. v2.0.0	大更新，还没微调
# 2023.10.28. v2.1.0	修复bug，增加计算信号长度的功能
# 2023.10.28. v2.2.0	修复bug，改变计算信号长度的算法
# 2023.10.29. v2.3.0	增加 overall signal 计算
# 2023.11.01. v2.3.1	增加 signal_only 选项
# 2023.11.02. v2.3.2	修改 sample signal 和 group mean 的计算顺序
# 2023.11.04. v2.3.3	修复 overall score 标准化时排序错误的问题
# 2023.11.05. v2.3.4	修复判断单边溢出信号时的列名选取错误
# 2023.11.13. v2.3.5	微调 track score
# 2023.12.05. v2.3.6	candidates 增加 cleavage site，修正 alignment 有 deletion 会错位的 bug
# 2023.12.05. v2.3.7	用 cleavage site 代替 midpoint # 还没改完
# 2023.12.07. v2.3.8	df_score 增加 df_exp, df_ctr 各自列。修复没 df_ctr 时的 bug。track score 用 proximal
# 2023.12.09. v2.4.0	为了兼顾 proximal 和 overall，当 normalized overall signal 高于 2 时，增加 overall signal 的加分
# 2023.12.09. v2.5.0	尝试新的加权位置
# 2023.12.10. v2.6.0	加入 trackseq v4 的计算分支，即考虑 Region 内的 positive_pct，避免短而尖锐的信号
# 2023.12.10. v2.6.1	有些非特异信号数值很大，如果在 control 组是大负数，可能导致减 control 后假高信号，因此给负数一个 clip
# 2023.12.30. v2.7.0	增加 X_offplot 模块，用于绘图
# 2023.12.31. v2.7.1	control 的负数值 clip 由 -5 改为 -1，进一步减少假阳性。另外不加 overall 了
# 2024.01.01. v2.7.2	权重改为 proximal + pct = 1 + 1. 防信号外溢假阳性标准由<0改为<=0
# 2024.01.02. v2.7.3	flank regions 默认值改为 1000 2000 3000 5000。之前 control 的负数值 clip 相当于直接在 final score，现在改为每个单独 clip 后重新算 score，默认值为 CtrClip=-0.5
# 2024.01.03. v2.7.4	更新了 blacklist.bed
# 2024.01.04. v2.7.5	更新了 hg38 blacklist.bed
# 2024.01.12. v2.7.6	修复小bug，输出 fdr 改为 <0.05。
# 2024.01.23. v2.7.7	Snakefile_offtracker: add --fixedStep to bigwigCompare for not merging neighbouring bins with equal values.
# 2024.02.01. v2.7.8	逐步添加 X_offplot.py 功能
# 2024.06.02. v2.7.9	添加 offtracker_plot.py
# 2024.06.03. v2.7.10	修复 bugs，offtable 添加 threshold = 2 的分界
# 2024.06.04. v2.7.11	readme 修改
# 2024.11.19. v2.7.12	offtracker_candidates.py 新增 --pam_location 参数指定 upstream 或 downstream，用于非 Cas9 情况
# 2025.04.25. v2.8.0	修复了 offtracker candidates 会把小写序列转换成 N 的 bug
# 2025.05.22. v2.9.0	翻新部分代码结构
# 2025.06.05. v2.10.0	增加了QC模块。保留了负数score的记录，并在plot时显示为红字。增加了 "--ignore_chr" 用于跳过common chr过滤。
# 2025.06.17. v2.10.6   修复翻新代码结构导致的bug