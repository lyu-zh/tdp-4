"""
筛选 merged_orders_by_date_location.csv 中的订单数据
筛选掉经纬度不在指定范围内的数据（纬度23-23.2，经度113-113.5）
保留的数据按照原来的相对顺序排列
"""

import pandas as pd
import os

def filter_orders_by_coordinates(input_file, output_file=None, 
                                  lat_min=23.0, lat_max=23.2,
                                  lon_min=113.0, lon_max=113.5):
    """
    筛选订单数据，保留指定经纬度范围内的数据
    
    参数:
        input_file: 输入CSV文件路径
        output_file: 输出CSV文件路径，如果为None则自动生成
        lat_min: 最小纬度，默认23.0
        lat_max: 最大纬度，默认23.2
        lon_min: 最小经度，默认113.0
        lon_max: 最大经度，默认113.5
    """
    # 读取CSV文件
    print(f"正在读取文件: {input_file}")
    df = pd.read_csv(input_file, encoding='utf-8-sig')
    
    print(f"原始数据行数: {len(df)}")
    print(f"原始数据列: {df.columns.tolist()}")
    
    # 检查必要的列是否存在
    required_cols = ['经度', '纬度']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"CSV文件中缺少必要的列: {missing_cols}")
    
    # 转换经纬度为数值类型
    print("\n正在转换数据类型...")
    df['经度'] = pd.to_numeric(df['经度'], errors='coerce')
    df['纬度'] = pd.to_numeric(df['纬度'], errors='coerce')
    
    # 删除经纬度为NaN的行
    initial_count = len(df)
    df = df.dropna(subset=['经度', '纬度'])
    valid_count = len(df)
    removed_nan = initial_count - valid_count
    if removed_nan > 0:
        print(f"删除了 {removed_nan} 行包含无效经纬度的数据")
    
    print(f"\n筛选条件:")
    print(f"  纬度范围: {lat_min} - {lat_max}")
    print(f"  经度范围: {lon_min} - {lon_max}")
    
    # 筛选数据（保留在指定范围内的数据）
    filtered_df = df[
        (df['纬度'] >= lat_min) & 
        (df['纬度'] <= lat_max) & 
        (df['经度'] >= lon_min) & 
        (df['经度'] <= lon_max)
    ].copy()
    
    # 保持原来的相对顺序（pandas的筛选操作已经保持了原始顺序）
    # 重置索引以保持连续性，但不改变行的顺序
    filtered_df = filtered_df.reset_index(drop=True)
    
    print(f"\n筛选结果:")
    print(f"  筛选后数据行数: {len(filtered_df)}")
    print(f"  筛选掉的数据行数: {len(df) - len(filtered_df)}")
    print(f"  保留比例: {len(filtered_df)/len(df)*100:.2f}%")
    
    # 如果没有指定输出文件，自动生成
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_dir = os.path.dirname(input_file)
        output_file = os.path.join(output_dir, f"{base_name}_filtered.csv")
    
    # 保存结果
    print(f"\n正在保存结果到: {output_file}")
    filtered_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    
    print(f"\n筛选完成！")
    print(f"前5行数据预览:")
    print(filtered_df.head())
    
    return filtered_df

def main():
    """主函数"""
    # 获取脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, 'data', 'test')
    
    # 文件路径
    input_file = os.path.join(data_dir, 'merged_orders_by_date_location.csv')
    output_file = os.path.join(data_dir, 'merged_orders_by_date_location_filtered_new.csv')
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"输入文件不存在: {input_file}")
    
    # 执行筛选
    filter_orders_by_coordinates(
        input_file=input_file,
        output_file=output_file,
        lat_min=23.0,
        lat_max=23.06,
        lon_min=113.4,
        lon_max=113.48
    )

if __name__ == '__main__':
    main()
