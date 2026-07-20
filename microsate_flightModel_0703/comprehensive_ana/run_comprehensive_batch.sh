#!/bin/bash
# 批量跑 comprehensive_ana：遍历 data_microsat，result_microsat 里已有同名大文件夹则跳过
#
# 用法: bash run_comprehensive_batch.sh
# 可选: bash run_comprehensive_batch.sh .pkg   # 只处理指定后缀

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ONE="$SCRIPT_DIR/run_comprehensive.sh"
DATA_DIR="$SCRIPT_DIR/data_secondFSE"
RESULT_DIR="$SCRIPT_DIR/result_secondFSE"

if [[ -n "$1" ]]; then
	FILE_EXTS=("$1")
else
	FILE_EXTS=(".pkg" ".dat" ".bin")
fi

if [[ ! -f "$RUN_ONE" ]]; then
	echo "找不到 $RUN_ONE"
	exit 1
fi

mkdir -p "$DATA_DIR" "$RESULT_DIR"

BATCH_LOCK="$SCRIPT_DIR/.batch.lock"
exec 9>"$BATCH_LOCK"
if ! flock -n 9; then
	echo "另一个 batch 正在运行，退出"
	exit 0
fi

has_result_dir() {
	local name="$1"
	local ext="$2"
	[[ -d "$RESULT_DIR/$name" ]] && return 0
	# 兼容旧版输出目录名带 .dat 后缀
	if [[ "$ext" == ".dat" && -d "$RESULT_DIR/${name}.dat" ]]; then
		return 0
	fi
	return 1
}

n_skip=0
n_run=0
n_fail=0

echo "数据目录: $DATA_DIR"
echo "结果目录: $RESULT_DIR"
echo "后缀: ${FILE_EXTS[*]}"
echo "=========================================="

for ext in "${FILE_EXTS[@]}"; do
	for binfile in "$DATA_DIR"/*"$ext"; do
		[[ -f "$binfile" ]] || continue

		binbasename="$(basename "$binfile" "$ext")"
		outroot="$RESULT_DIR/$binbasename"

		if has_result_dir "$binbasename" "$ext"; then
			echo "[skip] $binbasename  (已存在 $outroot)"
			((n_skip++)) || true
			continue
		fi

		echo "=========================================="
		echo "[run]  $binfile"
		echo "输出:   $outroot"
		echo "=========================================="

		if bash "$RUN_ONE" "$binfile"; then
			((n_run++)) || true
		else
			echo "[warn] $binbasename 处理失败"
			((n_fail++)) || true
		fi
	done
done

echo "=========================================="
echo "批量完成: 跳过=$n_skip  新处理=$n_run  失败=$n_fail"
echo "结果目录: $RESULT_DIR"
