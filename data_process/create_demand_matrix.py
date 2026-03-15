"""
根据merged_orders_by_date_location.csv创建需求矩阵
行：独特的经纬度坐标组合
列：发车日期
值：0或1（1表示该日期在该经纬度有需求，0表示没有）
"""

import pandas as pd
import numpy as np
import os

def create_demand_matrix(input_file, output_file=None, index_file=None):
    """
    创建需求矩阵CSV文件
    
    参数:
        input_file: 输入CSV文件路径
        output_file: 输出CSV文件路径，如果为None则自动生成
        index_file: 包含序号信息的CSV文件路径，如果为None则自动查找
    """
    # 读取CSV文件
    print(f"正在读取文件: {input_file}")
    df = pd.read_csv(input_file, encoding='utf-8-sig')
    
    print(f"原始数据行数: {len(df)}")
    print(f"列名: {df.columns.tolist()}")
    
    # 检查必要的列是否存在
    required_cols = ['发车日期', '经度', '纬度']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"CSV文件中缺少必要的列: {missing_cols}")
    
    # 数据预处理
    print("\n步骤1: 数据预处理...")
    print("-" * 80)
    
    # 转换发车日期为日期类型
    df['发车日期'] = pd.to_datetime(df['发车日期'], errors='coerce')
    
    # 转换经纬度为数值类型
    df['经度'] = pd.to_numeric(df['经度'], errors='coerce')
    df['纬度'] = pd.to_numeric(df['纬度'], errors='coerce')
    
    # 删除无效数据
    initial_count = len(df)
    df = df.dropna(subset=['发车日期', '经度', '纬度'])
    valid_count = len(df)
    removed_count = initial_count - valid_count
    
    print(f"原始数据行数: {initial_count:,}")
    print(f"有效数据行数: {valid_count:,} (删除了 {removed_count:,} 行无效数据)")
    
    # 提取独特的经纬度组合
    print("\n步骤2: 提取独特的经纬度组合...")
    print("-" * 80)
    unique_coords = df[['经度', '纬度']].drop_duplicates().reset_index(drop=True)
    print(f"独特经纬度组合数: {len(unique_coords):,}")
    
    # 尝试加载序号信息
    print("\n步骤2.5: 加载序号信息...")
    print("-" * 80)
    coord_index_map = None
    
    # 如果未指定序号文件，使用默认路径
    if index_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(script_dir)
        # 默认查找 data/test/unique_coordinates_list_filtered.csv，如果不存在则查找 data/unique_coordinates_list.csv
        default_index_file1 = os.path.join(parent_dir, "data", "test", "selected_low_ratio_points_top20.csv")
        default_index_file2 = os.path.join(parent_dir, "data", "unique_coordinates_list.csv")
        
        if os.path.exists(default_index_file1):
            index_file = default_index_file1
        elif os.path.exists(default_index_file2):
            index_file = default_index_file2
        else:
            index_file = None
    
    # 尝试读取序号文件
    if index_file and os.path.exists(index_file):
        try:
            print(f"正在读取序号文件: {index_file}")
            index_df = pd.read_csv(index_file, encoding='utf-8-sig')
            
            # 检查是否有序号列和经纬度列
            if '序号' in index_df.columns and '经度' in index_df.columns and '纬度' in index_df.columns:
                # 创建经纬度到序号的映射
                # 使用字符串匹配来避免浮点数精度问题
                index_df['经度'] = pd.to_numeric(index_df['经度'], errors='coerce')
                index_df['纬度'] = pd.to_numeric(index_df['纬度'], errors='coerce')
                
                # 创建映射字典，使用四舍五入到6位小数来匹配
                coord_index_map = {}
                for _, row in index_df.iterrows():
                    lon = round(row['经度'], 6)
                    lat = round(row['纬度'], 6)
                    coord_index_map[(lon, lat)] = int(row['序号'])
                
                print(f"成功加载 {len(coord_index_map):,} 个坐标的序号信息")
            else:
                print(f"警告: 序号文件缺少必要的列（需要：序号、经度、纬度）")
                coord_index_map = None
        except Exception as e:
            print(f"警告: 读取序号文件失败: {e}")
            coord_index_map = None
    else:
        if index_file:
            print(f"序号文件不存在: {index_file}")
        print("将使用自动生成的序号（从1开始）")
        coord_index_map = None
    
    # 提取独特的发车日期
    print("\n步骤3: 提取独特的发车日期...")
    print("-" * 80)
    unique_dates = sorted(df['发车日期'].unique())
    print(f"独特发车日期数: {len(unique_dates):,}")
    print(f"日期范围: {unique_dates[0]} 到 {unique_dates[-1]}")
    
    # 创建需求矩阵
    print("\n步骤4: 创建需求矩阵...")
    print("-" * 80)
    
    # 创建坐标标识列（用于分组）
    df['坐标标识'] = df.apply(lambda row: f"{row['经度']}_{row['纬度']}", axis=1)
    
    # 使用pivot_table创建矩阵（更高效的方法）
    print("正在创建矩阵...")
    # 先创建一个辅助列，表示有需求（值为1）
    df['需求'] = 1
    
    # 使用pivot_table创建矩阵
    matrix_df = df.pivot_table(
        index=['经度', '纬度'],
        columns='发车日期',
        values='需求',
        aggfunc='max',  # 如果有多个值，取最大值（实际上都是1）
        fill_value=0  # 没有需求的填0
    )
    
    # 重置索引，使经纬度成为普通列
    matrix_df = matrix_df.reset_index()
    
    # 【重要修复】确保包含所有基准坐标点
    # 如果加载了序号文件，需要确保所有坐标点都在矩阵中
    if coord_index_map:
        print("正在检查并补充缺失的坐标点...")
        # 从序号文件中获取所有基准坐标点
        base_coords = set()
        for (lon, lat) in coord_index_map.keys():
            base_coords.add((round(lon, 6), round(lat, 6)))
        
        # 获取矩阵中已有的坐标点
        matrix_coords = set()
        for _, row in matrix_df.iterrows():
            lon = round(row['经度'], 6)
            lat = round(row['纬度'], 6)
            matrix_coords.add((lon, lat))
        
        # 找出缺失的坐标点
        missing_coords = base_coords - matrix_coords
        if missing_coords:
            print(f"发现 {len(missing_coords)} 个缺失的坐标点，正在补充...")
            # 获取所有日期列
            date_cols = [col for col in matrix_df.columns if col not in ['经度', '纬度']]
            
            # 为缺失的坐标点创建行（所有日期列都是0）
            missing_rows = []
            for lon, lat in missing_coords:
                row_data = {'经度': lon, '纬度': lat}
                for date_col in date_cols:
                    row_data[date_col] = 0
                missing_rows.append(row_data)
            
            # 将缺失的行添加到矩阵
            if missing_rows:
                missing_df = pd.DataFrame(missing_rows)
                # 确保列的顺序一致
                missing_df = missing_df[matrix_df.columns]
                # 确保日期列的数据类型是整数
                for date_col in date_cols:
                    if date_col in missing_df.columns:
                        missing_df[date_col] = missing_df[date_col].astype(int)
                matrix_df = pd.concat([matrix_df, missing_df], ignore_index=True)
                print(f"已补充 {len(missing_rows)} 个坐标点")
    
    # 将日期列名转换为字符串格式
    new_columns = {}
    for col in matrix_df.columns:
        if col not in ['经度', '纬度']:
            # 如果是datetime对象，转换为字符串
            if isinstance(col, pd.Timestamp):
                new_columns[col] = col.strftime('%Y-%m-%d')
            else:
                new_columns[col] = str(col)
        else:
            new_columns[col] = col
    
    matrix_df = matrix_df.rename(columns=new_columns)
    
    # 确保所有日期列都存在（如果某些日期没有数据，需要添加）
    date_columns = [date.strftime('%Y-%m-%d') for date in unique_dates]
    for date_col in date_columns:
        if date_col not in matrix_df.columns:
            matrix_df[date_col] = 0
    
    # 重新排序列（经度、纬度在前，然后是日期列）
    cols_order = ['经度', '纬度'] + sorted([col for col in matrix_df.columns if col not in ['经度', '纬度']])
    matrix_df = matrix_df[cols_order]
    
    # 确保数据类型为整数
    for col in date_columns:
        if col in matrix_df.columns:
            matrix_df[col] = matrix_df[col].astype(int)
    
    # 创建输出DataFrame
    print("\n步骤5: 准备输出数据...")
    print("-" * 80)
    
    result_df = matrix_df.copy()
    
    # 添加序号列作为第一列
    print("正在添加序号列...")
    if coord_index_map:
        # 使用从文件加载的序号
        def get_index(row):
            lon = round(row['经度'], 6)
            lat = round(row['纬度'], 6)
            return coord_index_map.get((lon, lat), None)
        
        result_df['序号'] = result_df.apply(get_index, axis=1)
        
        # 检查是否有未匹配的坐标
        unmatched = result_df['序号'].isna().sum()
        if unmatched > 0:
            print(f"警告: 有 {unmatched} 个坐标未在序号文件中找到，将使用自动生成的序号")
            # 为未匹配的坐标生成序号（从最大序号+1开始）
            max_index = result_df['序号'].max() if result_df['序号'].notna().any() else 0
            result_df.loc[result_df['序号'].isna(), '序号'] = range(max_index + 1, max_index + 1 + unmatched)
        
        result_df['序号'] = result_df['序号'].astype(int)
    else:
        # 使用自动生成的序号（从1开始）
        result_df['序号'] = range(1, len(result_df) + 1)
    
    # 重新排序列，将序号放在第一列
    cols_order = ['序号', '经度', '纬度'] + sorted([col for col in result_df.columns 
                                                     if col not in ['序号', '经度', '纬度']])
    result_df = result_df[cols_order]
    
    # 按照序号排序
    print("正在按照序号排序...")
    result_df = result_df.sort_values('序号').reset_index(drop=True)
    
    # 统计信息
    date_cols = [col for col in result_df.columns if col not in ['序号', '经度', '纬度']]
    matrix_values = result_df[date_cols].values
    print(f"输出矩阵形状: {result_df.shape}")
    print(f"总需求数: {matrix_values.sum():,}")
    print(f"平均每个坐标点的需求天数: {matrix_values.sum(axis=1).mean():.2f}")
    print(f"平均每天的需求坐标点数: {matrix_values.sum(axis=0).mean():.2f}")
    
    # 确定输出文件路径
    if output_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(script_dir)
        output_file = os.path.join(parent_dir, "data", "demand_matrix_filtered_new.csv")
    
    # 保存结果
    print(f"\n步骤6: 保存结果...")
    print("-" * 80)
    result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"结果已保存到: {output_file}")
    
    return result_df

if __name__ == "__main__":
    # 设置文件路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    
    # 输入文件：可以根据需要修改
    input_file = os.path.join(parent_dir, "data", "test", "merged_orders_by_date_location_low_ratio.csv")
    
    # 可选：指定序号文件（如果为None，将自动查找）
    index_file = None  # 例如: os.path.join(parent_dir, "data", "test", "unique_coordinates_list_filtered.csv")
    
    # 可选：指定输出文件（如果为None，将使用默认路径）
    output_file = os.path.join(parent_dir, "data", "test", "demand_matrix_low_ratio.csv")  # 例如: os.path.join(parent_dir, "data", "demand_matrix_filtered.csv")
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误: 找不到文件 {input_file}")
        print(f"请检查文件路径是否正确")
    else:
        # 创建需求矩阵
        create_demand_matrix(input_file, output_file=output_file, index_file=index_file)
