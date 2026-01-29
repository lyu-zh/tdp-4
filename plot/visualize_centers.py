"""
可视化中心点和区域划分
- 读取 unique_coordinates_list 和 selected_centers_p200
- 中心点有边框
- 为每个中心点计算最近的 n/k 个点，并画透明圆
- 根据 demand_matrix 的 2022-07-17 列，用不同颜色区分0和1的点
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os

# 获取脚本所在目录的绝对路径
script_dir = os.path.dirname(os.path.abspath(__file__))
# 项目根目录（plot的上一级目录）
project_root = os.path.dirname(script_dir)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 数据文件路径（基于项目根目录）
data_dir = os.path.join(project_root, 'data', 'test')
coords_file = os.path.join(data_dir, 'unique_coordinates_list_filtered_new.csv')
centers_file = os.path.join(data_dir, 'selected_centers_filtered_new_p3.csv')
demand_file = os.path.join(data_dir, 'demand_matrix_filtered_new.csv')
output_dir = os.path.join(script_dir, 'pic')
target_date = '2023-09-03'

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

print("正在读取数据...")
print(f"数据目录: {data_dir}")
print(f"坐标文件: {coords_file}")
print(f"文件是否存在: {os.path.exists(coords_file)}")

# 读取坐标数据
if not os.path.exists(coords_file):
    raise FileNotFoundError(f"找不到坐标文件: {coords_file}")
coords_df = pd.read_csv(coords_file)
print(f"总点数: {len(coords_df)}")

# 特殊标记的点序号
special_point_ids = []
print(f"特殊标记的点序号: {special_point_ids}")

# 找到特殊点的坐标
special_point_coords = []
print(f"坐标文件列名: {coords_df.columns.tolist()}")

for point_id in special_point_ids:
    # 尝试不同的列名来查找序号
    found = False
    idx = None
    if '序号' in coords_df.columns:
        idx = coords_df[coords_df['序号'] == point_id].index
        found = len(idx) > 0
    elif coords_df.columns[0] in ['序号', 'id', 'ID', 'index', 'Index'] or coords_df.columns[0].lower() == 'id':
        idx = coords_df[coords_df.iloc[:, 0] == point_id].index
        found = len(idx) > 0
    else:
        # 尝试第一列（通常是序号）
        idx = coords_df[coords_df.iloc[:, 0] == point_id].index
        found = len(idx) > 0
    
    if found:
        row = coords_df.iloc[idx[0]]
        # 尝试不同的列名来获取经度纬度
        if '经度' in row.index and '纬度' in row.index:
            lon, lat = float(row['经度']), float(row['纬度'])
        elif 'longitude' in row.index and 'latitude' in row.index:
            lon, lat = float(row['longitude']), float(row['latitude'])
        else:
            # 假设第二列是经度，第三列是纬度
            lon, lat = float(row.iloc[1]), float(row.iloc[2])
        
        special_point_coords.append((lon, lat))
        print(f"  找到点 {point_id}: 索引={idx[0]}, 坐标=({lon:.6f}, {lat:.6f})")
    else:
        print(f"  警告: 未找到点 {point_id}")

# 读取中心点数据
centers_df = pd.read_csv(centers_file)
print(f"中心点数量: {len(centers_df)}")

# 读取需求矩阵
demand_df = pd.read_csv(demand_file)
print(f"需求矩阵行数: {len(demand_df)}")

# 计算每个中心点应包含的点数
n = len(coords_df)  # 总点数
k = len(centers_df)  # 中心点数
points_per_center = n // k
# points_per_center = 0
print(f"每个中心点应包含的点数: {points_per_center}")

# 提取坐标
all_points = coords_df[['经度', '纬度']].values
center_points = centers_df[['经度', '纬度']].values

# 检查目标日期列是否存在
if target_date not in demand_df.columns:
    print(f"警告: 找不到日期列 '{target_date}'，将使用第一列日期数据")
    # 找到第一个日期列（跳过经度、纬度列）
    date_cols = [col for col in demand_df.columns if col not in ['经度', '纬度']]
    if date_cols:
        target_date = date_cols[0]
        print(f"使用日期列: {target_date}")
    else:
        print("错误: 找不到任何日期列")
        exit(1)

# 创建需求字典：根据经度纬度匹配需求值
demand_dict = {}
for idx, row in demand_df.iterrows():
    lon = row['经度']
    lat = row['纬度']
    demand_value = row[target_date]
    # 使用经度纬度作为键（转换为字符串以便匹配）
    key = (float(lon), float(lat))
    demand_dict[key] = int(demand_value) if pd.notna(demand_value) else 0

print(f"需求数据点数量: {len(demand_dict)}")

# 为每个点匹配需求值
point_demands = []
for point in all_points:
    lon, lat = point[0], point[1]
    # 尝试精确匹配
    key = (float(lon), float(lat))
    if key in demand_dict:
        point_demands.append(demand_dict[key])
    else:
        # 如果精确匹配失败，尝试找最近的点（容差0.0001度，约11米）
        found = False
        for dkey, dval in demand_dict.items():
            if abs(dkey[0] - lon) < 0.0001 and abs(dkey[1] - lat) < 0.0001:
                point_demands.append(dval)
                found = True
                break
        if not found:
            point_demands.append(0)  # 默认值为0

point_demands = np.array(point_demands)
print(f"匹配到需求值的点数: {np.sum(point_demands == 1)} 个为1, {np.sum(point_demands == 0)} 个为0")

# 为每个中心点找到对应的demand值
print("正在匹配中心点的需求值...")
center_demands = []
for center in center_points:
    lon, lat = center[0], center[1]
    # 尝试精确匹配
    key = (float(lon), float(lat))
    if key in demand_dict:
        center_demands.append(demand_dict[key])
    else:
        # 如果精确匹配失败，尝试找最近的点（容差0.0001度，约11米）
        found = False
        for dkey, dval in demand_dict.items():
            if abs(dkey[0] - lon) < 0.0001 and abs(dkey[1] - lat) < 0.0001:
                center_demands.append(dval)
                found = True
                break
        if not found:
            center_demands.append(0)  # 默认值为0

center_demands = np.array(center_demands)
print(f"中心点需求值: {np.sum(center_demands == 1)} 个为1, {np.sum(center_demands == 0)} 个为0")

# 计算每个中心点到所有点的距离，并找到最近的 n/k 个点
print("正在计算每个中心点的最近邻点...")
center_regions = []  # 存储每个中心点的区域信息

for i, center in enumerate(center_points):
    # 计算到所有点的距离（欧氏距离）
    distances = np.sqrt(np.sum((all_points - center) ** 2, axis=1))
    
    # 找到最近的 points_per_center 个点的索引
    nearest_indices = np.argsort(distances)[:points_per_center]
    nearest_distances = distances[nearest_indices]
    
    # 找到最远点的距离（用于画圆）
    max_distance = np.max(nearest_distances)
    
    center_regions.append({
        'center': center,
        'nearest_indices': nearest_indices,
        'max_distance': max_distance
    })
    
    if (i + 1) % 20 == 0:
        print(f"  已处理 {i + 1}/{len(center_points)} 个中心点")

print("正在绘制图形...")
# 创建图形
fig, ax = plt.subplots(figsize=(16, 12))

# 首先绘制所有点，根据需求值用不同颜色
# 值为0的点
mask_0 = point_demands == 0
ax.scatter(all_points[mask_0, 0], all_points[mask_0, 1], 
           c='lightblue', s=15, alpha=0.6, label=f'需求为0的点 ({np.sum(mask_0)})', zorder=1)

# 值为1的点
mask_1 = point_demands == 1
ax.scatter(all_points[mask_1, 0], all_points[mask_1, 1], 
           c='orange', s=15, alpha=0.6, label=f'需求为1的点 ({np.sum(mask_1)})', zorder=1)

# 为每个中心点画透明圆
for region in center_regions:
    center = region['center']
    max_dist = region['max_distance']
    
    # 画透明圆
    circle = Circle((center[0], center[1]), max_dist, 
                   fill=True, alpha=0.15, facecolor='gray', edgecolor='none', zorder=2)
    ax.add_patch(circle)

# 绘制中心点（带边框，根据需求值设置颜色，大小与普通点一致）
# 需求为0的中心点（淡蓝色）
center_mask_0 = center_demands == 0
if np.any(center_mask_0):
    ax.scatter(center_points[center_mask_0, 0], center_points[center_mask_0, 1], 
               c='lightblue', s=15, edgecolors='black', linewidths=0.5, 
               marker='o', label=f'中心点(需求=0) ({np.sum(center_mask_0)})', zorder=3)

# 需求为1的中心点（橙色）
center_mask_1 = center_demands == 1
if np.any(center_mask_1):
    ax.scatter(center_points[center_mask_1, 0], center_points[center_mask_1, 1], 
               c='orange', s=15, edgecolors='black', linewidths=0.5, 
               marker='o', label=f'中心点(需求=1) ({np.sum(center_mask_1)})', zorder=3)

# 为特殊标记点添加带箭头的文本标签（点本身会按照需求值正常显示，不需要单独绘制）
if special_point_coords:
    for point_id, (lon, lat) in zip(special_point_ids, special_point_coords):
        ax.annotate(f'点{point_id}', 
                   xy=(lon, lat), 
                   xytext=(15, 15), 
                   textcoords='offset points',
                   fontsize=11, 
                   fontweight='bold',
                   color='darkred',
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='yellow', alpha=0.8, edgecolor='darkred', linewidth=1.5),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', color='darkred', lw=2),
                   zorder=6)

# 设置标签和标题
ax.set_xlabel('经度', fontsize=12)
ax.set_ylabel('纬度', fontsize=12)
ax.set_title(f'中心点区域划分可视化 (日期: {target_date})', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3)

# 保存图片
output_file = os.path.join(output_dir, f'centers_visualization_{target_date.replace("-", "_")}.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"图片已保存到: {output_file}")

# 也保存一个PDF版本
output_file_pdf = os.path.join(output_dir, f'centers_visualization_{target_date.replace("-", "_")}.pdf')
plt.savefig(output_file_pdf, bbox_inches='tight')
print(f"PDF已保存到: {output_file_pdf}")

plt.close()
print("可视化完成！")
