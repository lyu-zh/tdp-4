import csv
import sys

def parse_value(val_str):
    """解析CSV中的值，返回数值或None"""
    val_str = val_str.strip()
    if val_str == "NaN" or val_str == "Null" or val_str == "":
        return None
    try:
        return float(val_str)
    except:
        return None

def analyze_null_points():
    """分析Null点是否对应最大的前几个值（按行分析）"""
    
    # 读取两个文件
    file1 = "output/travel_dist_dual_values_4_p5.csv"  # 使用3n/k，有Null标记
    file2 = "output/travel_dist_dual_values_1_p5.csv"  # 使用5n/k，所有点都有值或NaN
    
    data1 = []
    data2 = []
    
    # 读取文件1
    with open(file1, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header1 = next(reader)
        for row in reader:
            data1.append(row)
    
    # 读取文件2
    with open(file2, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header2 = next(reader)
        for row in reader:
            data2.append(row)
    
    # 检查列数是否一致
    if len(header1) != len(header2):
        print("错误：两个文件的列数不一致")
        return
    
    num_regions = len(header1) - 1  # 减去PointID列
    
    print("=" * 80)
    print("分析：travel_dist_dual_values_4_p5.csv中为Null的点")
    print("对应到travel_dist_dual_values_1_p5.csv中的值是否是该行（该点）取值最大的前几个")
    print("=" * 80)
    print()
    
    # 对每一行（每个点）进行分析
    for row_idx, (row1, row2) in enumerate(zip(data1, data2)):
        point_id = row1[0]
        
        # 收集该行在文件1中为Null的区域及其在文件2中的值
        null_regions = []
        row_values = []  # 该行在文件2中的所有值
        
        for region_idx in range(num_regions):
            col_idx = region_idx + 1  # +1因为第一列是PointID
            val1 = row1[col_idx].strip()
            val2 = parse_value(row2[col_idx])
            
            # 如果文件1中是Null，记录这个区域
            if val1 == "Null":
                if val2 is not None:
                    null_regions.append((region_idx, val2))
            
            # 收集文件2中该行的所有有效数值
            if val2 is not None:
                row_values.append((region_idx, val2))
        
        # 如果该行有Null点，进行分析
        if len(null_regions) > 0:
            # 按值从大到小排序
            row_values.sort(key=lambda x: x[1], reverse=True)
            null_values = [x[1] for x in null_regions]
            null_values.sort(reverse=True)
            
            # 找出文件1中该行非Null且非NaN的值（即在3n/k范围内被计算的值）
            computed_regions = []
            for region_idx in range(num_regions):
                col_idx = region_idx + 1
                val1 = row1[col_idx].strip()
                val2 = parse_value(row2[col_idx])
                if val1 != "Null" and val1 != "NaN" and val2 is not None:
                    computed_regions.append((region_idx, val2))
            
            print(f"点 {point_id} (行 {row_idx + 2}):")
            print(f"  - Null区域数量: {len(null_regions)}")
            print(f"  - Null区域的值（从大到小）: {[f'{v:.2f}' for v in null_values]}")
            print(f"  - 该行所有值的前{len(null_values)}个最大值: ", end="")
            top_values = [x[1] for x in row_values[:len(null_values)]]
            print([f'{v:.2f}' for v in top_values])
            
            # 检查是否匹配
            if len(null_values) == len(top_values):
                match = all(abs(null_values[i] - top_values[i]) < 1e-6 for i in range(len(null_values)))
                if match:
                    print(f"  [MATCH] Null点的值确实是该行最大的前{len(null_values)}个")
                else:
                    print(f"  [NO MATCH] 不匹配！")
                    print(f"     在3n/k范围内被计算的值: {[f'{v:.2f}' for _, v in computed_regions]}")
                    print(f"     差异详情：")
                    for i in range(min(len(null_values), len(top_values))):
                        if abs(null_values[i] - top_values[i]) >= 1e-6:
                            print(f"       位置{i}: Null值={null_values[i]:.2f}, 最大值={top_values[i]:.2f}")
            else:
                print(f"  [UNKNOWN] 数量不一致（Null={len(null_values)}, 最大值={len(top_values)}）")
            
            print()

if __name__ == "__main__":
    analyze_null_points()
