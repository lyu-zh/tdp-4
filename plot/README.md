# 可视化脚本说明

## 功能描述

`visualize_centers.py` 脚本用于可视化中心点和区域划分，主要功能包括：

1. **读取数据**：
   - `unique_coordinates_list.csv`: 所有点的坐标
   - `selected_centers_p200.csv`: 中心点坐标
   - `demand_matrix.csv`: 需求矩阵（包含日期列）

2. **可视化特性**：
   - 中心点用红色星形标记，带有黑色边框
   - 普通点根据需求值（0或1）用不同颜色显示：
     - 需求为0的点：浅蓝色
     - 需求为1的点：橙色
   - 为每个中心点绘制一个透明圆形区域，包含离它最近的 n/k 个点
     - n = 总点数
     - k = 中心点数
     - 每个中心点包含 n/k 个最近的点

3. **输出文件**：
   - PNG格式图片（高分辨率，300 DPI）
   - PDF格式图片
   - 保存在 `pic/` 文件夹中

## 使用方法

```bash
cd plot
python visualize_centers.py
```

## 配置参数

在脚本中可以修改以下参数：

- `target_date`: 目标日期（默认: "2022-07-17"）
- `data_dir`: 数据文件目录（默认: "../data"）
- `output_dir`: 输出图片目录（默认: "pic"）

## 输出说明

生成的图片文件名格式：
- `centers_visualization_YYYY_MM_DD.png`
- `centers_visualization_YYYY_MM_DD.pdf`

其中 `YYYY_MM_DD` 是目标日期。

## 依赖库

- pandas
- numpy
- matplotlib

安装依赖：
```bash
pip install pandas numpy matplotlib
```
