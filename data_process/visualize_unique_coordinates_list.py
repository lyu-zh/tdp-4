import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 获取脚本所在目录，确保无论从哪里运行都能找到数据文件
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

print("=" * 80)
print("可视化独特坐标列表")
print("=" * 80)

# 步骤1: 读取CSV文件
print("\n步骤1: 读取CSV文件...")
print("-" * 80)

# 使用相对于脚本目录的正确路径（向上返回到项目根目录）
csv_file = os.path.join(script_dir, '..', 'data', 'test', 'unique_coordinates_list_filtered_new.csv')
csv_file = os.path.normpath(csv_file)  # 规范化路径，处理 .. 和 . 符号

if not os.path.exists(csv_file):
    print(f"错误: 文件 {csv_file} 不存在！")
    exit(1)

try:
    print(f"正在读取 {csv_file}...")
    df = pd.read_csv(csv_file, encoding='utf-8-sig')
    print(f"成功读取 {len(df):,} 行数据")
    print(f"列名: {df.columns.tolist()}")
    
    # 检查经纬度列是否存在
    if '经度' not in df.columns or '纬度' not in df.columns:
        print("错误: 未找到'经度'或'纬度'列")
        print(f"可用列: {df.columns.tolist()}")
        exit(1)
    
    # 显示前几行数据
    print("\n前5行数据:")
    print(df.head())
    
except Exception as e:
    print(f"读取文件时出错: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# 步骤2: 数据清洗和验证
print("\n步骤2: 数据清洗和验证...")
print("-" * 80)

# 删除缺失值
df_valid = df[['序号', '经度', '纬度']].dropna()
print(f"有效数据行数: {len(df_valid):,}")

# 转换为数值类型
df_valid['经度'] = pd.to_numeric(df_valid['经度'], errors='coerce')
df_valid['纬度'] = pd.to_numeric(df_valid['纬度'], errors='coerce')

# 再次删除转换后的缺失值
df_valid = df_valid.dropna()
print(f"转换后有效数据行数: {len(df_valid):,}")

# 验证经纬度范围（中国大致范围）
df_valid = df_valid[
    (df_valid['纬度'] >= 15) & (df_valid['纬度'] <= 40) &
    (df_valid['经度'] >= 70) & (df_valid['经度'] <= 140)
]
print(f"范围验证后有效数据行数: {len(df_valid):,}")

# 提取经纬度列表
lons = df_valid['经度'].tolist()
lats = df_valid['纬度'].tolist()

# 步骤3: 数据统计
print("\n步骤3: 数据统计信息...")
print("-" * 80)
print(f"总坐标点数: {len(df_valid):,}")
print(f"纬度范围: {min(lats):.6f} 到 {max(lats):.6f}")
print(f"经度范围: {min(lons):.6f} 到 {max(lons):.6f}")
print(f"纬度跨度: {max(lats) - min(lats):.6f} 度")
print(f"经度跨度: {max(lons) - min(lons):.6f} 度")

# 检查是否有重复点
unique_coords = df_valid.drop_duplicates(subset=['经度', '纬度'])
duplicate_count = len(df_valid) - len(unique_coords)
if duplicate_count > 0:
    print(f"重复坐标点数: {duplicate_count:,}")
    print(f"独特坐标点数: {len(unique_coords):,}")
else:
    print("所有坐标点都是独特的")

# 步骤4: 可视化
print("\n步骤4: 生成可视化图表...")
print("-" * 80)

# 创建图表
fig, ax = plt.subplots(figsize=(16, 12))

# 绘制散点图
print("正在绘制散点图...")
scatter = ax.scatter(lons, lats, c='blue', s=15, alpha=0.6, edgecolors='none')

# 设置标签和标题
ax.set_xlabel('经度 (Longitude)', fontsize=14)
ax.set_ylabel('纬度 (Latitude)', fontsize=14)
ax.set_title(f'独特坐标列表分布图\n(共 {len(df_valid):,} 个坐标点)', fontsize=16, fontweight='bold')

# 添加网格
ax.grid(True, alpha=0.3, linestyle='--')

# 设置合适的比例尺
lat_range = max(lats) - min(lats)
lon_range = max(lons) - min(lons)

# 添加适当的边距（5-10%）
lat_margin = max(lat_range * 0.05, 0.01)
lon_margin = max(lon_range * 0.05, 0.01)

# 设置坐标范围
x_min = min(lons) - lon_margin
x_max = max(lons) + lon_margin
y_min = min(lats) - lat_margin
y_max = max(lats) + lat_margin

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# 不强制等比例，让matplotlib自动调整
ax.set_aspect('auto')

# 添加统计信息文本框
info_text = f'总坐标点数: {len(df_valid):,}\n'
if duplicate_count > 0:
    info_text += f'独特坐标点数: {len(unique_coords):,}\n'
info_text += f'纬度范围: [{min(lats):.4f}, {max(lats):.4f}]\n'
info_text += f'经度范围: [{min(lons):.4f}, {max(lons):.4f}]'
ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
        fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 保存图片
output_file = 'pic/unique_coordinates_list_visualization_center.png'
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n图表已保存为: {output_file}")

# 显示图表
plt.show()

print("\n" + "=" * 80)
print("可视化完成！")
print("=" * 80)
print(f"✓ 处理了 {len(df_valid):,} 个坐标点")
print(f"✓ 可视化图表已保存: {output_file}")
