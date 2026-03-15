"""
筛选 merged_orders_by_date_location.csv 中的订单数据
只保留经纬度落在 data/test/selected_low_ratio_points_top20.csv 中存在的点上的订单。
保留的数据按照原来的相对顺序排列。
"""

import pandas as pd
import os

# 坐标匹配时保留的小数位数，避免浮点误差导致漏匹配
COORD_DECIMALS = 12


def load_allowed_coordinates(points_file):
    """
    从 selected_low_ratio_points_top20.csv 读取允许的 (经度, 纬度) 集合。
    使用四舍五入到 COORD_DECIMALS 位作为键，避免浮点误差。
    """
    df = pd.read_csv(points_file, encoding='utf-8-sig')
    for col in ['经度', '纬度']:
        if col not in df.columns:
            raise ValueError(f"点文件缺少列: {col}")
    df['经度'] = pd.to_numeric(df['经度'], errors='coerce')
    df['纬度'] = pd.to_numeric(df['纬度'], errors='coerce')
    df = df.dropna(subset=['经度', '纬度'])
    allowed = set(
        (round(row['经度'], COORD_DECIMALS), round(row['纬度'], COORD_DECIMALS))
        for _, row in df.iterrows()
    )
    return allowed


def filter_orders_by_coordinates(
    input_file,
    output_file=None,
    points_file=None,
):
    """
    筛选订单数据，只保留经纬度在 points_file 所列点上的订单。

    参数:
        input_file: 输入订单 CSV 文件路径
        output_file: 输出 CSV 文件路径，若为 None 则自动生成
        points_file: 允许的坐标点 CSV 路径（须含 经度、纬度 列），若为 None 则用默认路径
    """
    if points_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        points_file = os.path.join(project_root, 'data', 'test', 'selected_low_ratio_points_top20.csv')

    if not os.path.exists(points_file):
        raise FileNotFoundError(f"坐标点文件不存在: {points_file}")

    print(f"正在读取允许的坐标点: {points_file}")
    allowed_coords = load_allowed_coordinates(points_file)
    print(f"允许的坐标点数量: {len(allowed_coords)}")

    print(f"\n正在读取订单文件: {input_file}")
    df = pd.read_csv(input_file, encoding='utf-8-sig')

    print(f"原始数据行数: {len(df)}")
    print(f"原始数据列: {df.columns.tolist()}")

    required_cols = ['经度', '纬度']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"订单文件中缺少必要的列: {missing_cols}")

    df['经度'] = pd.to_numeric(df['经度'], errors='coerce')
    df['纬度'] = pd.to_numeric(df['纬度'], errors='coerce')

    initial_count = len(df)
    df = df.dropna(subset=['经度', '纬度'])
    if len(df) < initial_count:
        print(f"删除了 {initial_count - len(df)} 行包含无效经纬度的数据")

    # 用四舍五入后的坐标判断是否在允许集合中（向量化，避免逐行 apply）
    df['_lon_key'] = df['经度'].round(COORD_DECIMALS)
    df['_lat_key'] = df['纬度'].round(COORD_DECIMALS)
    df['_key'] = list(zip(df['_lon_key'], df['_lat_key']))
    filtered_df = df[df['_key'].isin(allowed_coords)].copy()
    # 统计点表中有多少个坐标在订单中被找到（至少有一条订单）
    coords_found = filtered_df[['_lon_key', '_lat_key']].drop_duplicates()
    n_coords_found = len(coords_found)
    filtered_df = filtered_df.drop(columns=['_key', '_lon_key', '_lat_key'])
    filtered_df = filtered_df.reset_index(drop=True)

    print(f"\n筛选结果:")
    print(f"  点表中共 {len(allowed_coords)} 个坐标，其中在订单中被找到的坐标数: {n_coords_found}")
    print(f"  筛选后数据行数: {len(filtered_df)}")
    print(f"  筛选掉的数据行数: {len(df) - len(filtered_df)}")
    if len(df) > 0:
        print(f"  保留比例: {len(filtered_df)/len(df)*100:.2f}%")

    if output_file is None:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_dir = os.path.dirname(input_file)
        output_file = os.path.join(output_dir, f"{base_name}_filtered.csv")

    print(f"\n正在保存结果到: {output_file}")
    filtered_df.to_csv(output_file, index=False, encoding='utf-8-sig')

    print("\n筛选完成！")
    print("前5行数据预览:")
    print(filtered_df.head())

    return filtered_df


def main():
    """主函数"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, 'data', 'test')

    input_file = os.path.join(data_dir, 'merged_orders_by_date_location.csv')
    output_file = os.path.join(data_dir, 'merged_orders_by_date_location_low_ratio.csv')
    points_file = os.path.join(data_dir, 'selected_low_ratio_points_top20.csv')

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    filter_orders_by_coordinates(
        input_file=input_file,
        output_file=output_file,
        points_file=points_file,
    )


if __name__ == '__main__':
    main()
