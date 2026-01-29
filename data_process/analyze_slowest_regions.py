#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
分析 travel_dist_dual_values_by_date 文件夹中所有 log 文件，
统计最慢的 5 个区域的出现次数
"""

import os
import re
from collections import Counter
from pathlib import Path

def extract_slowest_regions(log_file_path):
    """
    从 log 文件中提取最慢的 5 个区域
    
    返回: 区域编号列表，例如 [134, 107, 20, 163, 69]
    """
    try:
        # 尝试多种编码方式
        encodings = ['utf-8', 'gbk', 'gb2312', 'latin-1']
        content = None
        for encoding in encodings:
            try:
                with open(log_file_path, 'r', encoding=encoding) as f:
                    content = f.read()
                break
            except UnicodeDecodeError:
                continue
        
        if content is None:
            return []
        
        # 查找"最慢的5个区域"部分
        pattern = r'最慢的5个区域:\s*\n((?:\s*区域\s+\d+:\s+[\d.]+s\s*\n){5})'
        match = re.search(pattern, content)
        
        if not match:
            return []
        
        # 提取区域编号
        region_pattern = r'区域\s+(\d+):'
        regions = re.findall(region_pattern, match.group(1))
        
        return [int(r) for r in regions]
    except Exception as e:
        print(f"处理文件 {log_file_path} 时出错: {e}")
        return []

def main():
    # 设置 log 文件目录
    log_dir = Path("output/travel_dist_dual_values_by_date")
    
    if not log_dir.exists():
        print(f"错误: 目录 {log_dir} 不存在")
        return
    
    # 统计所有区域的出现次数
    region_counter = Counter()
    processed_files = 0
    failed_files = 0
    
    # 遍历所有 .txt 文件
    log_files = list(log_dir.glob("*.txt"))
    print(f"找到 {len(log_files)} 个 log 文件")
    
    for log_file in log_files:
        regions = extract_slowest_regions(log_file)
        if regions:
            region_counter.update(regions)
            processed_files += 1
        else:
            failed_files += 1
            if failed_files <= 5:  # 只显示前 5 个失败的文件
                print(f"警告: 无法从 {log_file.name} 中提取区域信息")
    
    print(f"\n处理完成:")
    print(f"  成功处理: {processed_files} 个文件")
    print(f"  失败: {failed_files} 个文件")
    print(f"  总共统计到 {len(region_counter)} 个不同的区域")
    
    # 输出前 10 个出现次数最多的区域
    print(f"\n{'='*60}")
    print("前 10 个出现次数最多的最慢区域:")
    print(f"{'='*60}")
    print(f"{'排名':<6} {'区域编号':<12} {'出现次数':<12} {'占比':<12}")
    print(f"{'-'*60}")
    
    total_count = sum(region_counter.values())
    for rank, (region, count) in enumerate(region_counter.most_common(10), 1):
        percentage = (count / total_count) * 100
        print(f"{rank:<6} {region:<12} {count:<12} {percentage:.2f}%")
    
    # 保存结果到文件
    output_file = "slowest_regions_statistics.txt"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("最慢区域统计结果\n")
        f.write("="*60 + "\n")
        f.write(f"处理文件数: {processed_files}\n")
        f.write(f"失败文件数: {failed_files}\n")
        f.write(f"不同区域数: {len(region_counter)}\n")
        f.write(f"总记录数: {total_count}\n\n")
        f.write("前 10 个出现次数最多的最慢区域:\n")
        f.write("-"*60 + "\n")
        f.write(f"{'排名':<6} {'区域编号':<12} {'出现次数':<12} {'占比':<12}\n")
        f.write("-"*60 + "\n")
        for rank, (region, count) in enumerate(region_counter.most_common(10), 1):
            percentage = (count / total_count) * 100
            f.write(f"{rank:<6} {region:<12} {count:<12} {percentage:.2f}%\n")
    
    print(f"\n结果已保存到: {output_file}")

if __name__ == "__main__":
    main()
