# -*- coding: utf-8 -*-
"""
对比分析脚本：CC vs 支撑超平面Cut方法 Out of Sample 性能对比
1. 样本近似算法 (CC) - 来自CSV文件
2. 支撑超平面Cut求解DRCC模型 (DRCC+Solver) - 来自xlsx文件
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 获取脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 数据文件路径
base_path = os.path.join(os.path.dirname(script_dir), 'csv', 'csv_plot')

print("=" * 70)
print("CC vs 支撑超平面Cut方法 Out of Sample 性能对比分析")
print("=" * 70)
print(f"数据目录: {base_path}")

# ============================================================
# 读取两组数据
# ============================================================
print("\n" + "-" * 50)
print("读取数据文件...")
print("-" * 50)

# 1. 样本近似算法 (CC)
cc_path = os.path.join(base_path, 'GG20+GG50+GG80+CC.csv')
cc_df = pd.read_csv(cc_path)
print(f"\n1. 样本近似算法 (CC):")
print(f"   文件: GG20+GG50+GG80+CC.csv")
print(f"   总行数: {len(cc_df)}")
print(f"   列名: {cc_df.columns.tolist()}")

# 2. 支撑超平面Cut求解DRCC模型
solver_path = os.path.join(base_path, 'GG20+GG50+GG80+DRCC+Solver.xlsx')
solver_df = pd.read_excel(solver_path)
print(f"\n2. 支撑超平面Cut方法 (DRCC+Solver):")
print(f"   文件: GG20+GG50+GG80+DRCC+Solver.xlsx")
print(f"   总行数: {len(solver_df)}")
print(f"   列名: {solver_df.columns.tolist()}")

# ============================================================
# 数据处理 - 创建匹配键
# ============================================================
print("\n" + "-" * 50)
print("创建匹配键...")
print("-" * 50)

# 基础匹配列（所有数据集都应该有的列）
base_match_cols = ['InstanceName', 'NumUnits', 'NumRegions', 'RSD', 'r', 'gamma']

def create_key(df, cols):
    """创建用于匹配的键"""
    available_cols = [col for col in cols if col in df.columns]
    if len(available_cols) == 0:
        raise ValueError(f"没有找到匹配列！可用列: {df.columns.tolist()}")
    key_parts = [df[col].astype(str) for col in available_cols]
    return pd.Series(['_'.join(x) for x in zip(*key_parts)])

# 创建匹配键
cc_df['ExperimentKey'] = create_key(cc_df, base_match_cols)
solver_df['ExperimentKey'] = create_key(solver_df, base_match_cols)

print(f"使用匹配列: {base_match_cols}")
print(f"CC 键示例: {cc_df['ExperimentKey'].iloc[0]}")
print(f"支撑超平面Cut 键示例: {solver_df['ExperimentKey'].iloc[0]}")

# ============================================================
# 筛选成功求解的实验
# ============================================================
print("\n" + "-" * 50)
print("筛选成功求解的实验...")
print("-" * 50)

def filter_success(df, name):
    """筛选成功求解的数据"""
    # 检查可能的成功标志列
    if 'SolveSuccess' in df.columns:
        success_df = df[df['SolveSuccess'] == True].copy()
        print(f"  {name}: 使用 'SolveSuccess' 列筛选, 成功数 = {len(success_df)}")
    elif 'Status' in df.columns:
        # Status=2 通常表示OPTIMAL
        success_df = df[df['Status'] == 2].copy()
        print(f"  {name}: 使用 'Status==2' 筛选, 成功数 = {len(success_df)}")
    else:
        # 假设所有数据都成功
        success_df = df.copy()
        print(f"  {name}: 未找到成功标志列, 假设全部成功, 数量 = {len(success_df)}")
    return success_df

cc_success = filter_success(cc_df, 'CC')
solver_success = filter_success(solver_df, '支撑超平面Cut')

# ============================================================
# 找到两组数据共同成功的实验
# ============================================================
print("\n" + "-" * 50)
print("寻找共同成功的实验...")
print("-" * 50)

cc_keys = set(cc_success['ExperimentKey'])
solver_keys = set(solver_success['ExperimentKey'])

print(f"CC 唯一实验数: {len(cc_keys)}")
print(f"支撑超平面Cut 唯一实验数: {len(solver_keys)}")

common_keys = cc_keys & solver_keys
print(f"\n两组共同的实验数: {len(common_keys)}")

if len(common_keys) == 0:
    print("\n警告：没有找到两组数据都成功的实验！")
    print("\n尝试查看数据差异...")
    print(f"  CC独有键示例: {list(cc_keys - solver_keys)[:3]}")
    print(f"  支撑超平面Cut独有键示例: {list(solver_keys - cc_keys)[:3]}")
    exit(1)

# 筛选共同成功的实验
cc_common = cc_success[cc_success['ExperimentKey'].isin(common_keys)].copy()
solver_common = solver_success[solver_success['ExperimentKey'].isin(common_keys)].copy()

# 处理重复键 - 取每个键的第一个
cc_common = cc_common.drop_duplicates(subset='ExperimentKey', keep='first')
solver_common = solver_common.drop_duplicates(subset='ExperimentKey', keep='first')

# 按键排序确保对应
cc_common = cc_common.sort_values('ExperimentKey').reset_index(drop=True)
solver_common = solver_common.sort_values('ExperimentKey').reset_index(drop=True)

print(f"去重后用于对比的实验数: {len(cc_common)}")

# 验证匹配
assert len(cc_common) == len(solver_common), "数据行数不匹配！"
assert (cc_common['ExperimentKey'].values == solver_common['ExperimentKey'].values).all(), "CC和支撑超平面Cut实验顺序不匹配！"

# ============================================================
# 提取 Out of Sample 数据
# ============================================================
print("\n" + "-" * 50)
print("提取 Out of Sample 性能数据...")
print("-" * 50)

def get_oos_column(df, name):
    """获取Out of Sample列"""
    oos_cols = ['OutOfSamplePerformance', 'OOS', 'OutOfSample', 'out_of_sample']
    for col in oos_cols:
        if col in df.columns:
            print(f"  {name}: 使用列 '{col}'")
            return col
    raise ValueError(f"{name}: 未找到Out of Sample列！可用列: {df.columns.tolist()}")

oos_cc = cc_common[get_oos_column(cc_common, 'CC')].values
oos_solver = solver_common[get_oos_column(solver_common, '支撑超平面Cut')].values

print(f"\n数据统计:")
print(f"  CC           - 平均: {np.mean(oos_cc):.4f}, 最小: {np.min(oos_cc):.4f}, 最大: {np.max(oos_cc):.4f}")
print(f"  支撑超平面Cut - 平均: {np.mean(oos_solver):.4f}, 最小: {np.min(oos_solver):.4f}, 最大: {np.max(oos_solver):.4f}")

# ============================================================
# 绘图 - CC vs 支撑超平面Cut
# ============================================================
print("\n" + "-" * 50)
print("生成对比图表...")
print("-" * 50)

experiment_ids = range(1, len(cc_common) + 1)

# 创建单张图
fig, ax = plt.subplots(figsize=(16, 8))

ax.plot(experiment_ids, oos_cc, 'b-o', linewidth=1.5, markersize=4, 
        label='样本近似算法 (CC)', alpha=0.8)
ax.plot(experiment_ids, oos_solver, 'r-s', linewidth=1.5, markersize=4, 
        label='支撑超平面Cut方法 (DRCC)', alpha=0.8)

ax.set_xlabel('实验序号', fontsize=11)
ax.set_ylabel('Out of Sample Performance', fontsize=11)
# ax.set_title('CC vs 支撑超平面Cut方法 Out of Sample 性能对比\n(仅包含两种方法都成功求解的实验)',
             # fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.5, len(experiment_ids) + 0.5)

# 设置y轴范围
y_min = min(np.min(oos_cc), np.min(oos_solver))
y_max = max(np.max(oos_cc), np.max(oos_solver))
ax.set_ylim(max(0, y_min - 0.05), min(1.05, y_max + 0.02))
# ax.set_ylim(0.99, 1)

# 添加统计信息
stats_text = (f'CC: 平均={np.mean(oos_cc):.4f}\n'
              f'支撑超平面Cut: 平均={np.mean(oos_solver):.4f}\n'
              f'差异: {np.mean(oos_cc) - np.mean(oos_solver):+.4f}')
ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='bottom', horizontalalignment='left', 
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()

# 保存图表
output_path_png = os.path.join(script_dir, 'cc_vs_supportcut_oos_comparison.png')
output_path_pdf = os.path.join(script_dir, 'cc_vs_supportcut_oos_comparison.pdf')
plt.savefig(output_path_png, dpi=150, bbox_inches='tight')
plt.savefig(output_path_pdf, bbox_inches='tight')
print(f"\n图表已保存:")
print(f"  {output_path_png}")
print(f"  {output_path_pdf}")

# ============================================================
# 输出详细统计
# ============================================================
print("\n" + "=" * 70)
print("详细统计信息")
print("=" * 70)

print(f"\n1. Out of Sample Performance 统计:")
print(f"   样本近似算法 (CC):")
print(f"     平均值: {np.mean(oos_cc):.4f}")
print(f"     标准差: {np.std(oos_cc):.4f}")
print(f"     最小值: {np.min(oos_cc):.4f}")
print(f"     最大值: {np.max(oos_cc):.4f}")

print(f"\n   支撑超平面Cut方法 (DRCC+Solver):")
print(f"     平均值: {np.mean(oos_solver):.4f}")
print(f"     标准差: {np.std(oos_solver):.4f}")
print(f"     最小值: {np.min(oos_solver):.4f}")
print(f"     最大值: {np.max(oos_solver):.4f}")

print(f"\n2. CC vs 支撑超平面Cut方法对比:")
print(f"   平均值差异: {np.mean(oos_cc) - np.mean(oos_solver):+.4f}")
cc_better = np.sum(oos_cc > oos_solver)
solver_better = np.sum(oos_solver > oos_cc)
tie = np.sum(oos_cc == oos_solver)
print(f"   CC更优: {cc_better} ({cc_better/len(cc_common)*100:.1f}%)")
print(f"   支撑超平面Cut更优: {solver_better} ({solver_better/len(cc_common)*100:.1f}%)")
print(f"   相同: {tie} ({tie/len(cc_common)*100:.1f}%)")

plt.show()
