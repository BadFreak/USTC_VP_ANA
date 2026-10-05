#!/bin/bash
# 批量跑 comprehensive_ana：遍历 data_*，result_* 里已有同名大文件夹则跳过
#
# 用法: bash run_comprehensive_batch.sh
#       bash run_comprehensive_batch.sh .pkg
#       bash run_comprehensive_batch.sh -f --decode-only --modes cosmic
#       bash run_comprehensive_batch.sh --data-dir DIR --result-dir DIR
# 输入/输出目录由本脚本决定（或命令行覆盖），调用 run_comprehensive.sh 时显式传入。
# 选项:
#   -f, --force       已有结果也强制重跑（传给 run_comprehensive.sh）
#   --decode-only     只跑 vdecode（传给 run_comprehensive.sh）
#   --modes LIST      只跑指定模式，如 cosmic 或 cosmic,ped
#   --data-dir DIR    输入数据目录
#   --result-dir DIR  输出结果目录

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ONE="$SCRIPT_DIR/run_comprehensive.sh"
DATA_DIR="$SCRIPT_DIR/data_thirdFSE"
RESULT_DIR="$SCRIPT_DIR/result_thirdFSE"

FORCE_RUN=false
RUN_ARGS=()
FILE_EXT=""

while [[ $# -gt 0 ]]; do
	case "$1" in
		-f|--force)
			FORCE_RUN=true
			RUN_ARGS+=("-f")
			shift
			;;
		--decode-only)
			RUN_ARGS+=("--decode-only")
			shift
			;;
		--modes)
			shift
			[[ $# -gt 0 ]] || { echo "--modes 需要参数"; exit 1; }
			RUN_ARGS+=("--modes" "$1")
			shift
			;;
		--data-dir)
			shift
			[[ $# -gt 0 ]] || { echo "--data-dir 需要目录"; exit 1; }
			DATA_DIR="$1"
			shift
			;;
		--result-dir)
			shift
			[[ $# -gt 0 ]] || { echo "--result-dir 需要目录"; exit 1; }
			RESULT_DIR="$1"
			shift
			;;
		-h|--help)
			echo "用法: bash run_comprehensive_batch.sh [-f] [--decode-only] [--modes cosmic] [--data-dir DIR] [--result-dir DIR] [.dat|.pkg|.bin]"
			exit 0
			;;
		-*)
			echo "未知选项: $1"
			exit 1
			;;
		*)
			FILE_EXT="$1"
			shift
			;;
	esac
done

if [[ -n "$FILE_EXT" ]]; then
	FILE_EXTS=("$FILE_EXT")
else
	FILE_EXTS=(".pkg" ".dat" ".bin")
fi

if [[ ! -f "$RUN_ONE" ]]; then
	echo "找不到 $RUN_ONE"
	exit 1
fi

mkdir -p "$DATA_DIR" "$RESULT_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"
RESULT_DIR="$(cd "$RESULT_DIR" && pwd)"

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
echo "force=$FORCE_RUN args: ${RUN_ARGS[*]}"
echo "=========================================="

for ext in "${FILE_EXTS[@]}"; do
	for binfile in "$DATA_DIR"/*"$ext"; do
		[[ -f "$binfile" ]] || continue

		binbasename="$(basename "$binfile" "$ext")"
		outroot="$RESULT_DIR/$binbasename"

		if ! $FORCE_RUN && has_result_dir "$binbasename" "$ext"; then
			echo "[skip] $binbasename  (已存在 $outroot)"
			((n_skip++)) || true
			continue
		fi

		echo "=========================================="
		echo "[run]  $binfile"
		echo "输出:   $outroot"
		echo "=========================================="

		if bash "$RUN_ONE" "${RUN_ARGS[@]}" \
			--data-dir "$DATA_DIR" --result-dir "$RESULT_DIR" "$binfile"; then
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
