#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
将 output/travel_dist_dual_values_filtered_by_date_new_filled 目录中
所有 CSV 文件的数据除以 1000
"""

import csv
import os
import glob

def divide_csv_by_1000(csv_file_path):
    """
    将 CSV 文件中的所有数值数据除以 1000
    
    参数:
        csv_file_path: CSV 文件路径
    """
    # 读取 CSV 文件
    rows = []
    encodings = ['utf-8', 'utf-8-sig', 'gbk', 'gb2312', 'gb18030']
    file_read = False
    
    for encoding in encodings:
        try:
            with open(csv_file_path, 'r', encoding=encoding, newline='') as f:
                reader = csv.reader(f)
                # 读取表头
                header = next(reader)
                rows.append(header)
                
                # 读取数据行
                for row in reader:
                    if not row:  # 跳过空行
                        continue
                    
                    # 第一列是 PointID，保持不变
                    # 其他列都是数值，需要除以 1000
                    new_row = [row[0]]  # PointID
                    
                    for i in range(1, len(row)):
                        try:
                            value = float(row[i])
                            new_row.append(str(value / 1000.0))
                        except (ValueError, IndexError):
                            # 如果无法转换为浮点数，保持原值
                            new_row.append(row[i])
                    
                    rows.append(new_row)
                
                file_read = True
                break  # 成功读取，退出循环
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(f"读取文件 {csv_file_path} 时出错: {e}")
            return False
    
    if not file_read:
        print(f"错误: 无法读取文件 {csv_file_path}")
        return False
    
    # 写回文件
    try:
        with open(csv_file_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(rows)
        return True
    except Exception as e:
        print(f"写入文件 {csv_file_path} 时出错: {e}")
        return False

def process_all_csv_files(directory_path):
    """
    处理目录中的所有 CSV 文件
    
    参数:
        directory_path: 目录路径
    """
    # 获取所有 CSV 文件
    csv_pattern = os.path.join(directory_path, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"警告: 在目录 {directory_path} 中未找到 CSV 文件")
        return
    
    print(f"找到 {len(csv_files)} 个 CSV 文件")
    
    success_count = 0
    fail_count = 0
    
    for csv_file in csv_files:
        file_name = os.path.basename(csv_file)
        print(f"正在处理: {file_name}...", end=' ')
        
        if divide_csv_by_1000(csv_file):
            print("成功")
            success_count += 1
        else:
            print("失败")
            fail_count += 1
    
    print(f"\n处理完成！")
    print(f"成功: {success_count} 个文件")
    print(f"失败: {fail_count} 个文件")

if __name__ == "__main__":
    # 设置目录路径 - 使用相对路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.abspath(os.path.join(script_dir, "../output/travel_dist_dual_values_filtered_by_date_new_filled"))
    
    print(f"目标目录: {target_dir}")
    
    # 检查目录是否存在
    if not os.path.exists(target_dir):
        print(f"错误: 目录不存在: {target_dir}")
        exit(1)
    
    if not os.path.isdir(target_dir):
        print(f"错误: 路径不是目录: {target_dir}")
        exit(1)
    
    # 处理所有 CSV 文件
    process_all_csv_files(target_dir)
    
    print("\n所有文件处理完成！")
