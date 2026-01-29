"""
检查 selected_centers_p200.csv 中的中心点序号和 unique_coordinates_list.csv 中对应序号的点的经纬度是否一致
"""

import pandas as pd
import numpy as np
import os

# 获取脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'data')

# 读取文件
centers_file = os.path.join(data_dir, 'selected_centers_p200.csv')
coords_file = os.path.join(data_dir, 'unique_coordinates_list.csv')

print("正在读取文件...")
centers_df = pd.read_csv(centers_file)
coords_df = pd.read_csv(coords_file)

print(f"中心点数量: {len(centers_df)}")
print(f"坐标点数量: {len(coords_df)}")
print()

# 创建坐标字典，以序号为键
coords_dict = {}
for idx, row in coords_df.iterrows():
    point_id = int(row['序号'])
    lon = float(row['经度'])
    lat = float(row['纬度'])
    coords_dict[point_id] = (lon, lat)

print("开始检查中心点坐标一致性...")
print("=" * 80)

mismatches = []
matches = 0
tolerance = 1e-10  # 浮点数比较容差（精确匹配）
coords_tolerance = 1e-6  # 经纬度匹配容差（约0.1米）

for idx, center_row in centers_df.iterrows():
    center_id = int(center_row['序号'])
    center_lon = float(center_row['经度'])
    center_lat = float(center_row['纬度'])
    
    # 在坐标文件中查找对应序号
    if center_id in coords_dict:
        coord_lon, coord_lat = coords_dict[center_id]
        
        # 比较经纬度
        lon_diff = abs(center_lon - coord_lon)
        lat_diff = abs(center_lat - coord_lat)
        
        if lon_diff < tolerance and lat_diff < tolerance:
            matches += 1
        else:
            # 尝试通过经纬度匹配（可能序号不对应，但坐标对应）
            found_by_coords = False
            for coord_id, (c_lon, c_lat) in coords_dict.items():
                if abs(center_lon - c_lon) < coords_tolerance and abs(center_lat - c_lat) < coords_tolerance:
                    found_by_coords = True
                    mismatches.append({
                        '序号': center_id,
                        '中心点经度': center_lon,
                        '中心点纬度': center_lat,
                        '坐标文件序号': coord_id,
                        '坐标文件经度': c_lon,
                        '坐标文件纬度': c_lat,
                        '经度差': abs(center_lon - c_lon),
                        '纬度差': abs(center_lat - c_lat),
                        '匹配方式': '通过坐标匹配（序号不匹配）'
                    })
                    break
            
            if not found_by_coords:
                mismatches.append({
                    '序号': center_id,
                    '中心点经度': center_lon,
                    '中心点纬度': center_lat,
                    '坐标文件序号': center_id,
                    '坐标文件经度': coord_lon,
                    '坐标文件纬度': coord_lat,
                    '经度差': lon_diff,
                    '纬度差': lat_diff,
                    '匹配方式': '序号匹配但坐标不一致'
                })
    else:
        mismatches.append({
            '序号': center_id,
            '中心点经度': center_lon,
            '中心点纬度': center_lat,
            '坐标文件经度': None,
            '坐标文件纬度': None,
            '经度差': None,
            '纬度差': None,
            '错误': '在坐标文件中未找到该序号'
        })

print(f"检查完成！")
print(f"一致的点数: {matches}/{len(centers_df)}")
print(f"不一致的点数: {len(mismatches)}")
print()

if mismatches:
    print("=" * 80)
    print("不一致的点详情:")
    print("=" * 80)
    for mismatch in mismatches:
        print(f"\n序号: {mismatch['序号']}")
        print(f"  中心点坐标: ({mismatch['中心点经度']:.15f}, {mismatch['中心点纬度']:.15f})")
        if mismatch['坐标文件经度'] is not None:
            print(f"  坐标文件序号: {mismatch.get('坐标文件序号', mismatch['序号'])}")
            print(f"  坐标文件坐标: ({mismatch['坐标文件经度']:.15f}, {mismatch['坐标文件纬度']:.15f})")
            print(f"  经度差: {mismatch['经度差']:.2e}")
            print(f"  纬度差: {mismatch['纬度差']:.2e}")
            if '匹配方式' in mismatch:
                print(f"  匹配方式: {mismatch['匹配方式']}")
        else:
            print(f"  错误: {mismatch.get('错误', '未知错误')}")
else:
    print("✓ 所有中心点的坐标都与坐标文件中的对应点一致！")

print()
print("=" * 80)
