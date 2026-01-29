# -*- coding: utf-8 -*-
"""
对比分析脚本：SupportCut方法 vs RelativeBalance求解器直接求解
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
supportcut_df = pd.read_csv('D:/NJU/dro理论学习/机会约束区域划分/TDP-pre/TDP-pre/output/csv/GG20+SupportCut.csv')
d1_df = pd.read_csv('D:/NJU/dro理论学习/机会约束区域划分/TDP-pre/TDP-pre/output/csv/GG20+RelativeBalance+D1.csv')
d2_df = pd.read_csv('D:/NJU/dro理论学习/机会约束区域划分/TDP-pre/TDP-pre/output/csv/GG20+RelativeBalance+D2.csv')

# 合并D1和D2数据
solver_df = pd.concat([d1_df, d2_df], ignore_index=True)

print("SupportCut数据行数:", len(supportcut_df))
print("RelativeBalance合并数据行数:", len(solver_df))

# 定义实验的唯一标识列（用于匹配对应实验）
match_cols = ['InstanceName', 'NumUnits', 'NumRegions', 'RSD', 'r', 'gamma', 'Scenarios', 'UseD1', 'Delta1', 'Delta2']

# 创建唯一标识
supportcut_df['ExperimentKey'] = supportcut_df[match_cols].astype(str).agg('_'.join, axis=1)
solver_df['ExperimentKey'] = solver_df[match_cols].astype(str).agg('_'.join, axis=1)

# 对于solver数据，选择UseExactMethod=false的（直接求解器求解）
solver_direct = solver_df[solver_df['UseExactMethod'] == False].copy()

print("\nSolver直接求解数据行数:", len(solver_direct))

# 筛选两个数据集中都成功求解的实验
supportcut_success = supportcut_df[supportcut_df['SolveSuccess'] == True].copy()
solver_success = solver_direct[solver_direct['SolveSuccess'] == True].copy()

print("\nSupportCut成功求解数:", len(supportcut_success))
print("Solver直接求解成功数:", len(solver_success))

# 找到共同成功的实验
common_keys = set(supportcut_success['ExperimentKey']) & set(solver_success['ExperimentKey'])
print("\n两种方法都成功求解的实验数:", len(common_keys))

# 筛选共同成功的实验
supportcut_common = supportcut_success[supportcut_success['ExperimentKey'].isin(common_keys)].copy()
solver_common = solver_success[solver_success['ExperimentKey'].isin(common_keys)].copy()

# 按ExperimentKey排序确保对应
supportcut_common = supportcut_common.sort_values('ExperimentKey').reset_index(drop=True)
solver_common = solver_common.sort_values('ExperimentKey').reset_index(drop=True)

# 验证匹配
assert len(supportcut_common) == len(solver_common), "数据行数不匹配！"
assert (supportcut_common['ExperimentKey'].values == solver_common['ExperimentKey'].values).all(), "实验顺序不匹配！"

print("\n最终用于对比的实验数:", len(supportcut_common))

# 准备绘图数据
experiment_ids = range(1, len(supportcut_common) + 1)
objective_diff = supportcut_common['Objective'].values - solver_common['Objective'].values
oos_cut = supportcut_common['OutOfSamplePerformance'].values
oos_solver = solver_common['OutOfSamplePerformance'].values
time_cut = supportcut_common['Runtime(s)'].values
time_solver = solver_common['Runtime(s)'].values

# 创建图形 - 三个子图
fig, axes = plt.subplots(3, 1, figsize=(14, 12))
fig.suptitle('SupportCut方法 vs 求解器直接求解 对比分析', fontsize=16, fontweight='bold')

# 子图1: 目标函数值之差 (折线图)
ax1 = axes[0]
ax1.plot(experiment_ids, objective_diff, 'b-o', linewidth=1.5, markersize=4, label='目标函数差值 (Cut - Solver)')
ax1.axhline(y=0, color='r', linestyle='--', alpha=0.5, label='零线')
ax1.set_xlabel('实验序号', fontsize=11)
ax1.set_ylabel('目标函数差值', fontsize=11)
ax1.set_title('目标函数值对比 (Cut方法 - Solver方法)', fontsize=12)
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0.5, len(experiment_ids) + 0.5)

# 添加统计信息
mean_diff = np.mean(objective_diff)
ax1.text(0.02, 0.95, f'平均差值: {mean_diff:.4f}', transform=ax1.transAxes, 
         fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 子图2: Out of Sample对比 (折线图)
ax2 = axes[1]
ax2.plot(experiment_ids, oos_cut, 'b-o', linewidth=1.5, markersize=4, label='SupportCut', alpha=0.8)
ax2.plot(experiment_ids, oos_solver, 'r-s', linewidth=1.5, markersize=4, label='求解器直接求解', alpha=0.8)
ax2.set_xlabel('实验序号', fontsize=11)
ax2.set_ylabel('Out of Sample Performance', fontsize=11)
ax2.set_title('Out of Sample性能对比', fontsize=12)
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0.5, len(experiment_ids) + 0.5)
ax2.set_ylim(0.9, 1.01)

# 添加平均值信息
mean_oos_cut = np.mean(oos_cut)
mean_oos_solver = np.mean(oos_solver)
ax2.text(0.02, 0.05, f'Cut平均OOS: {mean_oos_cut:.4f}\nSolver平均OOS: {mean_oos_solver:.4f}', 
         transform=ax2.transAxes, fontsize=10, verticalalignment='bottom', 
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 子图3: 计算时间对比 (双折线图)
ax3 = axes[2]
ax3.plot(experiment_ids, time_cut, 'b-o', linewidth=1.5, markersize=4, label='SupportCut', alpha=0.8)
ax3.plot(experiment_ids, time_solver, 'r-s', linewidth=1.5, markersize=4, label='求解器直接求解', alpha=0.8)
ax3.set_xlabel('实验序号', fontsize=11)
ax3.set_ylabel('计算时间 (秒)', fontsize=11)
ax3.set_title('计算时间对比', fontsize=12)
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0.5, len(experiment_ids) + 0.5)

# 添加平均时间信息
mean_time_cut = np.mean(time_cut)
mean_time_solver = np.mean(time_solver)
ax3.text(0.02, 0.95, f'Cut平均时间: {mean_time_cut:.3f}s\nSolver平均时间: {mean_time_solver:.3f}s', 
         transform=ax3.transAxes, fontsize=10, verticalalignment='top', 
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('comparison_analysis.png', dpi=150, bbox_inches='tight')
plt.savefig('comparison_analysis.pdf', bbox_inches='tight')
print("\n图表已保存为 comparison_analysis.png 和 comparison_analysis.pdf")

# 输出详细统计
print("\n" + "="*60)
print("详细统计信息")
print("="*60)
print(f"\n1. 目标函数差值 (Cut - Solver):")
print(f"   平均值: {np.mean(objective_diff):.4f}")
print(f"   标准差: {np.std(objective_diff):.4f}")
print(f"   最大值: {np.max(objective_diff):.4f}")
print(f"   最小值: {np.min(objective_diff):.4f}")

print(f"\n2. Out of Sample Performance:")
print(f"   Cut方法平均: {mean_oos_cut:.4f}")
print(f"   Solver方法平均: {mean_oos_solver:.4f}")

print(f"\n3. 计算时间:")
print(f"   Cut方法平均: {mean_time_cut:.3f}s")
print(f"   Solver方法平均: {mean_time_solver:.3f}s")
print(f"   Cut方法总时间: {np.sum(time_cut):.3f}s")
print(f"   Solver方法总时间: {np.sum(time_solver):.3f}s")

plt.show()

