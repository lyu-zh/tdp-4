"""
使用K-means方法对坐标点进行聚类，聚成3类
在每一类中挑选最接近聚类中心的点作为中心点
输出三个中心点的序号和坐标
"""

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import os

# 获取脚本所在目录和项目根目录
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)  # 项目根目录（data_process的父目录）
data_dir = os.path.join(project_root, 'data', 'test')

# 文件路径
input_file = os.path.join(data_dir, 'selected_low_ratio_points_top20.csv')
output_file = os.path.join(data_dir, 'selected_centers_low_ratio_p3.csv')

print("=" * 80)
print("K-means聚类分析")
print("=" * 80)

# 步骤1: 读取CSV文件
print("\n步骤1: 读取CSV文件...")
print("-" * 80)

try:
    df = pd.read_csv(input_file, encoding='utf-8-sig')
    print(f"成功读取 {len(df):,} 行数据")
    print(f"列名: {df.columns.tolist()}")
    
    # 检查必要的列是否存在
    if '经度' not in df.columns or '纬度' not in df.columns:
        print("错误: 未找到'经度'或'纬度'列")
        print(f"可用列: {df.columns.tolist()}")
        exit(1)
    
    if '序号' not in df.columns:
        print("警告: 未找到'序号'列，将自动创建")
        df.insert(0, '序号', range(1, len(df) + 1))
    
    # 显示前几行数据
    print("\n前5行数据:")
    print(df.head())
    
except Exception as e:
    print(f"读取文件时出错: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# 步骤2: 数据准备
print("\n步骤2: 数据准备...")
print("-" * 80)

# 提取经纬度数据
coordinates = df[['经度', '纬度']].values
print(f"坐标数据形状: {coordinates.shape}")

# 检查是否有缺失值
if pd.isna(coordinates).any():
    print("警告: 发现缺失值，将删除包含缺失值的行")
    valid_mask = ~pd.isna(coordinates).any(axis=1)
    df = df[valid_mask].reset_index(drop=True)
    coordinates = df[['经度', '纬度']].values
    print(f"清理后坐标数据形状: {coordinates.shape}")

# 步骤3: K-means聚类
print("\n步骤3: 执行K-means聚类（k=3）...")
print("-" * 80)

n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
cluster_labels = kmeans.fit_predict(coordinates)

# 将聚类标签添加到数据框
df['聚类标签'] = cluster_labels

print(f"聚类完成，共分为 {n_clusters} 类")
print(f"各类别点数统计:")
cluster_counts = pd.Series(cluster_labels).value_counts().sort_index()
for cluster_id, count in cluster_counts.items():
    print(f"  类别 {cluster_id}: {count} 个点")

# 步骤4: 在每一类中选择最接近聚类中心的点
print("\n步骤4: 选择各类别的中心点...")
print("-" * 80)

cluster_centers = kmeans.cluster_centers_
selected_centers = []

for cluster_id in range(n_clusters):
    # 获取当前类别的所有点
    cluster_mask = cluster_labels == cluster_id
    cluster_points = coordinates[cluster_mask]
    cluster_df = df[cluster_mask].copy()
    
    # 计算每个点到聚类中心的距离
    center = cluster_centers[cluster_id]
    distances = np.sqrt(np.sum((cluster_points - center) ** 2, axis=1))
    
    # 找到距离最小的点（最接近聚类中心的点）
    min_distance_idx = np.argmin(distances)
    selected_point = cluster_df.iloc[min_distance_idx]
    
    selected_centers.append({
        '序号': selected_point['序号'],
        '经度': selected_point['经度'],
        '纬度': selected_point['纬度'],
        '聚类标签': cluster_id,
        '到聚类中心距离': distances[min_distance_idx]
    })
    
    print(f"类别 {cluster_id}:")
    print(f"  聚类中心坐标: ({center[0]:.6f}, {center[1]:.6f})")
    print(f"  选择的中心点序号: {selected_point['序号']}")
    print(f"  选择的中心点坐标: ({selected_point['经度']:.6f}, {selected_point['纬度']:.6f})")
    print(f"  到聚类中心的距离: {distances[min_distance_idx]:.6f}")

# 步骤5: 保存结果
print("\n步骤5: 保存结果...")
print("-" * 80)

# 创建结果数据框（只包含序号、经度、纬度）
result_df = pd.DataFrame([
    {
        '序号': center['序号'],
        '经度': center['经度'],
        '纬度': center['纬度']
    }
    for center in selected_centers
])

# 按照聚类标签排序
result_df = result_df.sort_values(by='序号').reset_index(drop=True)

# 保存到CSV文件
result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
print(f"结果已保存到: {output_file}")
print(f"\n选中的3个中心点:")
print(result_df.to_string(index=False))

print("\n" + "=" * 80)
print("聚类分析完成！")
print("=" * 80)
print(f"✓ 处理了 {len(df):,} 个坐标点")
print(f"✓ 聚成 {n_clusters} 类")
print(f"✓ 选择了 {len(selected_centers)} 个中心点")
print(f"✓ 结果已保存: {output_file}")
