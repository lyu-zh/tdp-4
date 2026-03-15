#!/usr/bin/env bash
# 直接运行 DRTest1.java，不依赖 IDE
# 用法：在项目根目录执行 bash run_DRTest1.sh 或 ./run_DRTest1.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 检查 Java
if ! command -v javac &>/dev/null; then
    echo "[错误] 未找到 javac，请安装 JDK。"
    exit 1
fi
if ! command -v java &>/dev/null; then
    echo "[错误] 未找到 java，请安装 JDK。"
    exit 1
fi

# 路径：优先环境变量，否则默认
if [[ -n "$GUROBI_JAR" ]]; then
    :
elif [[ -n "$GUROBI_HOME" && -f "$GUROBI_HOME/lib/gurobi.jar" ]]; then
    export GUROBI_JAR="$GUROBI_HOME/lib/gurobi.jar"
else
    # Windows Git Bash 常见路径
    if [[ -f "D:/gurobi1300/win64/lib/gurobi.jar" ]]; then
        export GUROBI_JAR="D:/gurobi1300/win64/lib/gurobi.jar"
    elif [[ -f "/opt/gurobi/lib/gurobi.jar" ]]; then
        export GUROBI_JAR="/opt/gurobi/lib/gurobi.jar"
    else
        echo "[错误] 未找到 gurobi.jar。请设置 GUROBI_JAR 或 GUROBI_HOME。"
        exit 1
    fi
fi

JAMA_JAR="$SCRIPT_DIR/Jama-1.0.3.jar"
SRC="$SCRIPT_DIR/src"
OUT="$SCRIPT_DIR/out/production/TDP"

[[ -f "$JAMA_JAR" ]] || { echo "[错误] 未找到 Jama-1.0.3.jar"; exit 1; }
mkdir -p "$OUT"

# 编译
echo "正在编译..."
CP="$OUT:$GUROBI_JAR:$JAMA_JAR"
javac -encoding UTF-8 -cp "$CP" -sourcepath "$SRC" -d "$OUT" "$SRC/DRTest1.java"
echo "编译成功."
echo

# 运行
echo "运行 DRTest1..."
java -cp "$OUT:$GUROBI_JAR:$JAMA_JAR" DRTest1
