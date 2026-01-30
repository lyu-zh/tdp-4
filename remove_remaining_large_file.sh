#!/bin/bash
# 从历史中移除遗漏的大文件 data/merged_orders_by_date_location.csv 并再次推送
# 在项目根目录用 Git Bash 运行: bash remove_remaining_large_file.sh

set -e
echo "=== 1. 删除上次 filter-branch 的备份（否则无法再次运行）==="
git for-each-ref --format="%(refname)" refs/original/ 2>/dev/null | xargs -n 1 git update-ref -d 2>/dev/null || true

echo "=== 2. 暂存当前更改 ==="
git stash push -u -m "before remove remaining large file" 2>/dev/null || true

echo "=== 3. 从全部历史中移除 data/merged_orders_by_date_location.csv ==="
export FILTER_BRANCH_SQUELCH_WARNING=1
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch data/merged_orders_by_date_location.csv" \
  --prune-empty -- --all

echo "=== 4. 恢复暂存 ==="
git stash pop 2>/dev/null || true

echo "=== 5. 推送到 GitHub ==="
git push -u origin master --force

echo "=== 完成 ==="
