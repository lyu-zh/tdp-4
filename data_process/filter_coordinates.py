"""
筛选 unique_coordinates_list.csv 中经纬度在指定范围内的数据
保留：纬度 23-23.2，经度 113-113.5
筛选后重新排序并重新编号（从0开始连续编号）
"""

import pandas as pd
import os

# 获取脚本所在目录和项目根目录
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)  # 项目根目录（data_process的父目录）
data_dir = os.path.join(project_root, 'data', 'test')

# 文件路径
input_file = os.path.join(data_dir, 'unique_coordinates_list.csv')
output_file = os.path.join(data_dir, 'unique_coordinates_list_filtered_new.csv')

print("正在读取文件...")
df = pd.read_csv(input_file)

print(f"原始数据行数: {len(df)}")
print(f"原始数据列: {df.columns.tolist()}")

# 筛选条件：纬度 23-23.2，经度 113-113.5
lat_min, lat_max = 23.0, 23.06
lon_min, lon_max = 113.4, 113.48

print(f"\n筛选条件:")
print(f"  纬度范围: {lat_min} - {lat_max}")
print(f"  经度范围: {lon_min} - {lon_max}")

# 筛选数据
filtered_df = df[
    (df['纬度'] >= lat_min) & 
    (df['纬度'] <= lat_max) & 
    (df['经度'] >= lon_min) & 
    (df['经度'] <= lon_max)
].copy()

print(f"\n筛选后数据行数: {len(filtered_df)}")
print(f"筛选掉的数据行数: {len(df) - len(filtered_df)}")

# 对所有点重新排序（按照序号从小到大）
if '序号' in filtered_df.columns:
    filtered_df = filtered_df.sort_values(by='序号', ascending=True).reset_index(drop=True)
    print(f"\n已按照序号排序")
elif '经度' in filtered_df.columns and '纬度' in filtered_df.columns:
    # 如果没有序号列，按照经度、纬度排序
    filtered_df = filtered_df.sort_values(by=['经度', '纬度'], ascending=True).reset_index(drop=True)
    print(f"\n已按照经度、纬度排序")
else:
    print(f"\n警告: 未找到合适的排序列，保持原顺序")

# 重新编号（从0开始连续编号）
start_index = 1  # 如果需要从1开始，改为 1
if '序号' in filtered_df.columns:
    filtered_df['序号'] = range(start_index, start_index + len(filtered_df))
    print(f"已重新编号，从 {start_index} 开始，共 {len(filtered_df)} 个点")
else:
    # 如果没有序号列，创建新的序号列
    filtered_df.insert(0, '序号', range(start_index, start_index + len(filtered_df)))
    print(f"已创建序号列并重新编号，从 {start_index} 开始，共 {len(filtered_df)} 个点")

# 保存结果
filtered_df.to_csv(output_file, index=False, encoding='utf-8-sig')

print(f"\n筛选结果已保存到: {output_file}")
print(f"前5行数据预览:")
print(filtered_df.head())
