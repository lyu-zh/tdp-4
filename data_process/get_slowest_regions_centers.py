#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
获取前10个最慢区域对应的中心点序号
"""

import csv
from pathlib import Path

def load_centers(centers_file):
    """
    加载中心点数据
    返回: 字典 {区域编号: 中心点序号}
    """
    centers = {}
    # 尝试多种编码方式
    encodings = ['utf-8-sig', 'utf-8', 'gbk', 'gb2312']
    content = None
    encoding_used = None
    
    for encoding in encodings:
        try:
            with open(centers_file, 'r', encoding=encoding) as f:
                content = f.read()
            encoding_used = encoding
            break
        except Exception:
            continue
    
    if content is None:
        raise ValueError(f"无法读取文件 {centers_file}")
    
    # 解析CSV内容
    lines = content.strip().split('\n')
    if len(lines) < 2:
        raise ValueError(f"文件 {centers_file} 格式不正确")
    
    # 跳过标题行，从第二行开始
    for region_idx, line in enumerate(lines[1:], 0):
        parts = line.split(',')
        if len(parts) >= 1:
            try:
                center_index = int(parts[0].strip())
                centers[region_idx] = center_index
            except ValueError:
                continue
    
    return centers

def main():
    # 前10个最慢区域
    slowest_regions = [107, 32, 69, 57, 46, 9, 186, 125, 123, 189]
    
    # 加载中心点数据
    centers_file = "data/selected_centers_p200.csv"
    centers = load_centers(centers_file)
    
    print("="*70)
    print("前 10 个最慢区域对应的中心点序号:")
    print("="*70)
    print(f"{'排名':<6} {'区域编号':<12} {'中心点序号':<15} {'经度':<18} {'纬度':<18}")
    print("-"*70)
    
    # 读取中心点的详细信息
    center_details = {}
    encodings = ['utf-8-sig', 'utf-8', 'gbk', 'gb2312']
    content = None
    
    for encoding in encodings:
        try:
            with open(centers_file, 'r', encoding=encoding) as f:
                content = f.read()
            break
        except Exception:
            continue
    
    if content:
        lines = content.strip().split('\n')
        # 跳过标题行
        for region_idx, line in enumerate(lines[1:], 0):
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 3:
                try:
                    center_index = int(parts[0])
                    longitude = float(parts[1])
                    latitude = float(parts[2])
                    center_details[region_idx] = {
                        'index': center_index,
                        'longitude': longitude,
                        'latitude': latitude
                    }
                except ValueError:
                    continue
    
    for rank, region in enumerate(slowest_regions, 1):
        if region in center_details:
            info = center_details[region]
            print(f"{rank:<6} {region:<12} {info['index']:<15} {info['longitude']:<18.10f} {info['latitude']:<18.10f}")
        else:
            print(f"{rank:<6} {region:<12} {'未找到':<15}")
    
    # 保存结果到文件
    output_file = "slowest_regions_centers.txt"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("前 10 个最慢区域对应的中心点序号\n")
        f.write("="*70 + "\n")
        f.write(f"{'排名':<6} {'区域编号':<12} {'中心点序号':<15} {'经度':<18} {'纬度':<18}\n")
        f.write("-"*70 + "\n")
        for rank, region in enumerate(slowest_regions, 1):
            if region in center_details:
                info = center_details[region]
                f.write(f"{rank:<6} {region:<12} {info['index']:<15} {info['longitude']:<18.10f} {info['latitude']:<18.10f}\n")
            else:
                f.write(f"{rank:<6} {region:<12} {'未找到':<15}\n")
    
    print(f"\n结果已保存到: {output_file}")
    
    # 输出中心点序号列表（方便复制）
    print("\n中心点序号列表:")
    center_indices = [center_details[r]['index'] for r in slowest_regions if r in center_details]
    print(", ".join(map(str, center_indices)))

if __name__ == "__main__":
    main()
