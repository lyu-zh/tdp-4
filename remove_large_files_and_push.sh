#!/bin/bash
# 从 Git 历史中移除超过 GitHub 100MB 限制的大文件，并推送到 origin/master
# 请在项目根目录下用 Git Bash 运行: bash remove_large_files_and_push.sh

set -e
echo "=== 1. 添加 .gitignore 并暂存当前更改 ==="
if [ -f .gitignore ]; then
  git add .gitignore
else
  echo "提示: 未找到 .gitignore，请确保项目根目录已有 .gitignore 再运行"
fi
git stash push -u -m "before remove large files" 2>/dev/null || true

echo "=== 2. 从全部提交历史中移除大文件 ==="
export FILTER_BRANCH_SQUELCH_WARNING=1
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch data/merged_orders_by_date_location.csv data/test/merged_orders_by_date_location.csv data/delivery_2023-11-19_reorganized.csv" \
  --prune-empty -- --all

echo "=== 3. 恢复暂存的更改 ==="
git stash pop || true

echo "=== 4. 推送到 GitHub（会改写远程历史，因远程为空故直接 push）==="
git push -u origin master --force

echo "=== 完成 ==="
