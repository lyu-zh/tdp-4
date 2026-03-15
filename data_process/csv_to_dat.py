#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
将 unique_coordinates_list_filtered_new.csv 转换为 dat 格式
仿照 unique_coordinates_list_filtered.dat 的格式
"""

import csv
import os

def csv_to_dat(csv_file_path, dat_file_path, k=30):
    """
    将 CSV 文件转换为 dat 格式
    
    参数:
        csv_file_path: 输入的 CSV 文件路径
        dat_file_path: 输出的 dat 文件路径
        k: 区域数量（默认30）
    """
    # 读取 CSV 文件
    points = []
    # 尝试多种编码方式
    encodings = ['utf-8', 'utf-8-sig', 'gbk', 'gb2312', 'gb18030']
    file_read = False
    
    for encoding in encodings:
        try:
            with open(csv_file_path, 'r', encoding=encoding) as f:
                lines = f.readlines()
                if len(lines) < 2:
                    continue
                
                # 跳过表头（第一行）
                for line_idx, line in enumerate(lines[1:], start=2):
                    line = line.strip()
                    if not line:
                        continue
                    
                    # 按逗号分割
                    parts = line.split(',')
                    if len(parts) < 3:
                        print(f"警告: 第 {line_idx} 行数据不完整，跳过")
                        continue
                    
                    try:
                        # 第一列是序号（从1开始），转换为从0开始的id
                        point_id = int(parts[0]) - 1
                        longitude = float(parts[1])
                        latitude = float(parts[2])
                        points.append((point_id, longitude, latitude))
                    except (ValueError, IndexError) as e:
                        print(f"警告: 第 {line_idx} 行数据格式错误，跳过: {e}")
                        continue
                
                file_read = True
                break  # 成功读取，退出循环
        except UnicodeDecodeError as e:
            print(f"尝试编码 {encoding} 失败: {e}")
            continue
    
    if not file_read or not points:
        raise ValueError(f"无法读取 CSV 文件: {csv_file_path}，或文件中没有有效数据")
    
    n = len(points)  # 节点个数
    
    # 生成完全图的边（每个节点都与其他所有节点相连）
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            edges.append((i, j))
    
    m = len(edges)  # 边数
    
    print(f"节点数: {n}")
    print(f"区域数: {k}")
    print(f"边数: {m}")
    
    # 写入 dat 文件
    with open(dat_file_path, 'w', encoding='utf-8') as f:
        # 第一行：节点个数
        f.write(f"{n}\n")
        
        # 第二行：区域数量
        f.write(f"{k}\n")
        
        # 接下来 n 行：节点数据（id x y）
        for point_id, longitude, latitude in points:
            f.write(f"{point_id} {longitude} {latitude}\n")
        
        # 然后一行：边数
        f.write(f"{m}\n")
        
        # 接下来 m 行：边数据（a b）
        for a, b in edges:
            f.write(f"{a} {b}\n")
    
    print(f"成功生成 dat 文件: {dat_file_path}")

if __name__ == "__main__":
    # 设置文件路径 - 使用相对路径
    # 假设脚本在 data_process 目录下，CSV 文件在 data/test 目录下
    csv_file = "../data/test/selected_low_ratio_points_top20.csv"
    dat_file = "../data/test/selected_low_ratio_points_top20.dat"
    
    # 转换为绝对路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_file = os.path.abspath(os.path.join(script_dir, csv_file))
    dat_file = os.path.abspath(os.path.join(script_dir, dat_file))
    
    print(f"CSV 文件路径: {csv_file}")
    print(f"DAT 文件路径: {dat_file}")
    
    # 检查 CSV 文件是否存在
    if not os.path.exists(csv_file):
        print(f"错误: CSV 文件不存在: {csv_file}")
        exit(1)
    
    # 转换文件
    # 区域数量 k 可以根据需要调整，这里使用与原文件相同的值 30
    csv_to_dat(csv_file, dat_file, k=30)
    
    print("\n转换完成！")
