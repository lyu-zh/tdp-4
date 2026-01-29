#!/bin/bash

# TSP对偶值提取器运行脚本
# 使用方法: ./run_tsp_dual_extractor.sh

# 获取脚本所在目录（支持Windows Git Bash和Linux）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 设置颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== TSP对偶值提取器运行脚本 ===${NC}"

# 检查Java是否安装
if ! command -v javac &> /dev/null; then
    echo -e "${RED}错误: 未找到 javac 命令，请确保已安装 Java JDK${NC}"
    exit 1
fi

if ! command -v java &> /dev/null; then
    echo -e "${RED}错误: 未找到 java 命令，请确保已安装 Java JDK${NC}"
    exit 1
fi

# 显示Java版本
echo -e "${YELLOW}Java 版本:${NC}"
java -version
echo ""

# 设置类路径
# For Windows Git Bash, use Windows path format; for Linux, adjust path accordingly
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    # Windows Git Bash or Cygwin
    GUROBI_JAR="D:/gurobi1300/win64/lib/gurobi.jar"
else
    # Linux - adjust path as needed
    GUROBI_JAR="/opt/gurobi/lib/gurobi.jar"
fi
JAMA_JAR="$SCRIPT_DIR/Jama-1.0.3.jar"
SRC_DIR="$SCRIPT_DIR/src"
OUTPUT_DIR="$SCRIPT_DIR/output"
DATA_DIR="$SCRIPT_DIR/data"

# 检查必要的文件是否存在
if [ ! -f "$GUROBI_JAR" ]; then
    echo -e "${RED}错误: 未找到 gurobi.jar，期望位置: $GUROBI_JAR${NC}"
    if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
        echo -e "${YELLOW}提示: 如果是 Windows Git Bash，请确保 Gurobi 安装在 D:\\gurobi1300\\win64\\lib\\${NC}"
    else
        echo -e "${YELLOW}提示: 如果是 Linux，请修改脚本中的 GUROBI_JAR 路径${NC}"
    fi
    exit 1
fi

if [ ! -f "$JAMA_JAR" ]; then
    echo -e "${YELLOW}警告: 未找到 Jama-1.0.3.jar，期望位置: $JAMA_JAR${NC}"
    echo -e "${YELLOW}继续运行，但可能缺少某些依赖...${NC}"
fi

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 设置编译和运行的类路径
CLASSPATH="$GUROBI_JAR:$JAMA_JAR:$SRC_DIR"

# 编译Java文件
echo -e "${GREEN}正在编译 Java 文件...${NC}"
javac -encoding UTF-8 -cp "$CLASSPATH" -d "$SRC_DIR" "$SRC_DIR/TSPDualExtractor.java"

if [ $? -ne 0 ]; then
    echo -e "${RED}编译失败！请检查错误信息。${NC}"
    exit 1
fi

echo -e "${GREEN}编译成功！${NC}"
echo ""

# 检查数据文件是否存在
POINTS_FILE="$DATA_DIR/unique_coordinates_list.csv"
CENTERS_FILE="$DATA_DIR/selected_centers_p200.csv"

if [ ! -f "$POINTS_FILE" ]; then
    echo -e "${YELLOW}警告: 未找到数据文件: $POINTS_FILE${NC}"
    echo -e "${YELLOW}程序将使用默认路径，如果文件不存在会报错${NC}"
fi

if [ ! -f "$CENTERS_FILE" ]; then
    echo -e "${YELLOW}警告: 未找到中心点文件: $CENTERS_FILE${NC}"
    echo -e "${YELLOW}程序将使用默认路径，如果文件不存在会报错${NC}"
fi

# 运行程序
echo -e "${GREEN}开始运行 TSP对偶值提取器...${NC}"
echo -e "${YELLOW}注意: 这可能需要较长时间，请耐心等待...${NC}"
echo ""

# 设置最大堆内存（可根据需要调整）
# -Xmx4g 表示最大堆内存为4GB，可根据机器配置调整
java -Xmx4g -cp "$CLASSPATH" TSPDualExtractor

# 检查运行结果
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}程序运行完成！${NC}"
    echo -e "${GREEN}输出文件位于: $OUTPUT_DIR${NC}"
else
    echo ""
    echo -e "${RED}程序运行失败！请检查错误信息。${NC}"
    exit 1
fi
