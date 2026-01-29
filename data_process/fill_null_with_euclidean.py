"""
填充CSV文件中的Null值和NaN值

- Null值：替换为对应点到区域中心点的球面距离（米），使用Haversine公式
- NaN值：替换为0

数据来源：
- 所有点坐标：data/unique_coordinates_list.csv
- 区域中心点坐标：data/selected_centers_p200.csv
- 待处理CSV文件目录：output/travel_dist_dual_values_by_date/

使用方法：
    python fill_null_with_euclidean.py
    
输出：
    处理后的文件保存到 output/travel_dist_dual_values_by_date_filled/ 目录
"""

import os
import pandas as pd
import numpy as np
from math import radians, cos, sin, sqrt, atan2
import glob
import time
from multiprocessing import Pool, Manager, cpu_count
import multiprocessing

# 配置路径
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # data_process 目录
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)  # 项目根目录 (TDP-pre-5)
# 数据文件存储在项目根目录下的 data/test 文件夹中
DATA_DIR = os.path.join(PROJECT_ROOT, 'data', 'test')
# 输出目录在项目根目录下的 output 文件夹中
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'output', 'travel_dist_dual_values_filtered_by_date_new')
OUTPUT_FILLED_DIR = os.path.join(PROJECT_ROOT, 'output', 'travel_dist_dual_values_filtered_by_date_new_filled')

COORDS_FILE = os.path.join(DATA_DIR, 'unique_coordinates_list_filtered_new.csv')
CENTERS_FILE = os.path.join(DATA_DIR, 'selected_centers_filtered_new_p3.csv')

# 全局变量，用于多进程共享
DIST_MATRIX = None


def haversine_distance_meters(lon1, lat1, lon2, lat2):
    """
    使用Haversine公式计算两个经纬度点之间的球面距离（米）
    
    Haversine公式考虑了地球曲率，适用于任意两点间的距离计算。
    
    Args:
        lon1, lat1: 第一个点的经度和纬度（度）
        lon2, lat2: 第二个点的经度和纬度（度）
    
    Returns:
        两点之间的球面距离（米）
    """
    R = 6371000  # 地球平均半径（米）
    
    # 转换为弧度
    lat1_rad = radians(lat1)
    lat2_rad = radians(lat2)
    delta_lat = radians(lat2 - lat1)
    delta_lon = radians(lon2 - lon1)
    
    # Haversine公式
    a = sin(delta_lat / 2) ** 2 + cos(lat1_rad) * cos(lat2_rad) * sin(delta_lon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    
    return R * c


def load_point_coordinates(coords_file):
    """加载所有点的坐标"""
    df = pd.read_csv(coords_file)
    coords = {}
    for _, row in df.iterrows():
        point_id = int(row['序号']) - 1
        lon = float(row['经度'])
        lat = float(row['纬度'])
        coords[point_id] = (lon, lat)
    return coords


def load_center_coordinates(centers_file):
    """加载区域中心点的坐标"""
    df = pd.read_csv(centers_file)
    center_coords = {}
    for region_id, row in df.iterrows():
        lon = float(row['经度'])
        lat = float(row['纬度'])
        center_coords[region_id] = (lon, lat)
    return center_coords


def precompute_distance_matrix(point_coords, center_coords):
    """预计算所有点到所有中心点的欧氏距离矩阵"""
    num_points = max(point_coords.keys()) + 1
    num_regions = len(center_coords)
    
    dist_matrix = np.zeros((num_points, num_regions))
    
    print(f"预计算距离矩阵: {num_points} 点 × {num_regions} 区域...")
    
    for point_id in range(num_points):
        if point_id not in point_coords:
            continue
        p_lon, p_lat = point_coords[point_id]
        for region_id in range(num_regions):
            if region_id not in center_coords:
                continue
            c_lon, c_lat = center_coords[region_id]
            dist_matrix[point_id, region_id] = haversine_distance_meters(p_lon, p_lat, c_lon, c_lat)
    
    print("距离矩阵预计算完成")
    return dist_matrix


def init_worker(dist_matrix):
    """初始化工作进程，设置全局距离矩阵"""
    global DIST_MATRIX
    DIST_MATRIX = dist_matrix


def process_single_csv_worker(args):
    """
    工作进程函数，处理单个CSV文件
    - Null值：替换为球面距离
    - NaN值：替换为0
    """
    csv_path, output_path, file_idx, total_files = args
    global DIST_MATRIX
    
    filename = os.path.basename(csv_path)
    
    try:
        # 读取CSV文件（保持字符串格式以便检测Null和NaN）
        df = pd.read_csv(csv_path, dtype=str)
        
        null_count = 0
        nan_count = 0
        
        # 获取区域列
        region_columns = [col for col in df.columns if col.startswith('Region_')]
        
        # 遍历每一行（每一行对应一个点）
        for idx, row in df.iterrows():
            point_id = int(row['PointID'])
            
            if point_id >= DIST_MATRIX.shape[0]:
                continue
            
            for region_col in region_columns:
                region_id = int(region_col.replace('Region_', ''))
                
                if region_id >= DIST_MATRIX.shape[1]:
                    continue
                
                value = row[region_col]
                
                # 检查是否为Null或NaN
                if isinstance(value, str):
                    value_lower = value.strip().lower()
                    if value_lower == 'null':
                        # Null值：替换为球面距离
                        dist = DIST_MATRIX[point_id, region_id]
                        df.at[idx, region_col] = f"{dist:.10f}"
                        null_count += 1
                    elif value_lower == 'nan':
                        # NaN值：替换为0
                        df.at[idx, region_col] = '0'
                        nan_count += 1
                elif pd.isna(value):
                    # pandas NA值：替换为0
                    df.at[idx, region_col] = '0'
                    nan_count += 1
        
        # 保存处理后的文件
        df.to_csv(output_path, index=False)
        
        return (file_idx, filename, null_count, nan_count, None)
    
    except Exception as e:
        return (file_idx, filename, 0, 0, str(e))


def main():
    print("=" * 60)
    print("CSV文件Null/NaN值填充工具 (并行版本)")
    print("- Null值：替换为球面距离（米），使用Haversine公式")
    print("- NaN值：替换为0")
    print("=" * 60)
    
    # 检测CPU核心数
    num_cores = cpu_count()
    num_workers = min(num_cores, 20)  # 最多使用20个线程
    print(f"\n检测到 {num_cores} 个CPU核心，将使用 {num_workers} 个工作进程")
    
    start_time = time.time()
    
    # 创建输出目录
    os.makedirs(OUTPUT_FILLED_DIR, exist_ok=True)
    print(f"输出目录: {OUTPUT_FILLED_DIR}")
    
    # 加载坐标数据
    print("\n加载坐标数据...")
    point_coords = load_point_coordinates(COORDS_FILE)
    print(f"  加载了 {len(point_coords)} 个点的坐标")
    
    center_coords = load_center_coordinates(CENTERS_FILE)
    print(f"  加载了 {len(center_coords)} 个区域中心点的坐标")
    
    # 预计算距离矩阵
    print("\n预计算欧氏距离矩阵...")
    dist_matrix = precompute_distance_matrix(point_coords, center_coords)
    
    # 获取所有CSV文件
    csv_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, '*.csv')))
    total_files = len(csv_files)
    print(f"\n找到 {total_files} 个CSV文件待处理")
    
    if total_files == 0:
        print("没有找到CSV文件，程序结束")
        return
    
    # 准备任务参数
    tasks = []
    for idx, csv_path in enumerate(csv_files):
        filename = os.path.basename(csv_path)
        output_path = os.path.join(OUTPUT_FILLED_DIR, filename)
        tasks.append((csv_path, output_path, idx + 1, total_files))
    
    # 并行处理
    print(f"\n开始并行处理 ({num_workers} 个工作进程)...")
    print("-" * 60)
    
    total_null = 0
    total_nan = 0
    processed = 0
    errors = 0
    
    # 使用进程池并行处理
    with Pool(processes=num_workers, initializer=init_worker, initargs=(dist_matrix,)) as pool:
        # 使用imap_unordered获取结果，以便实时输出进度
        for result in pool.imap_unordered(process_single_csv_worker, tasks):
            file_idx, filename, null_count, nan_count, error = result
            
            if error:
                print(f"  [错误] {filename}: {error}")
                errors += 1
            else:
                total_null += null_count
                total_nan += nan_count
                processed += 1
                
                elapsed = time.time() - start_time
                rate = processed / elapsed if elapsed > 0 else 0
                remaining = (total_files - processed) / rate if rate > 0 else 0
                
                print(f"  [{processed}/{total_files}] {filename} - "
                      f"Null→距离: {null_count}, NaN→0: {nan_count} (剩余: {remaining:.0f} 秒)")
    
    elapsed = time.time() - start_time
    
    print("-" * 60)
    print("\n" + "=" * 60)
    print("处理完成！")
    print(f"  成功处理: {processed} 个文件")
    print(f"  处理错误: {errors} 个文件")
    print(f"  Null值→球面距离: {total_null:,}")
    print(f"  NaN值→0: {total_nan:,}")
    print(f"  总耗时: {elapsed:.1f} 秒")
    print(f"  平均速度: {processed / elapsed:.1f} 文件/秒")
    print(f"  输出目录: {OUTPUT_FILLED_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    # Windows下需要这个保护
    multiprocessing.freeze_support()
    main()
