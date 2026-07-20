#!/bin/bash
# CSI chn1 (CellID==100000) 对比：从 input/csi 填谱 -> 归一化叠加 + MPV 折线
# 无 onlyCalo。只在本目录写输出。
#
# 用法:
#   bash run_compare_csi_chn1.sh
#   bash run_compare_csi_chn1.sh -f
#   bash run_compare_csi_chn1.sh --plot-only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOT_PY="$SCRIPT_DIR/plot_csi_chn1_compare.py"
INPUT_DIR="$SCRIPT_DIR/input/csi"
OUTPUT_DIR="$SCRIPT_DIR/output/csi"
COMPARE_DIR="$OUTPUT_DIR/compare"

FORCE=false
PLOT_ONLY=false

while [[ $# -gt 0 ]]; do
	case "$1" in
		-f|--force) FORCE=true; shift ;;
		--plot-only) PLOT_ONLY=true; shift ;;
		-h|--help)
			echo "用法: bash run_compare_csi_chn1.sh [-f|--force] [--plot-only]"
			exit 0
			;;
		*) echo "未知选项: $1"; exit 1 ;;
	esac
done

# 无 onlyCalo
TESTS=(firstFSE secondFSE specialItem vibration)

echo "INPUT : $INPUT_DIR"
echo "OUTPUT: $OUTPUT_DIR"
echo "CSI channel: chn1 (CellID=100000)"

mkdir -p "$COMPARE_DIR" "$OUTPUT_DIR"

extra=()
$FORCE && extra+=(--force)
$PLOT_ONLY && extra+=(--plot-only)

python3 "$PLOT_PY" \
	--tests "${TESTS[@]}" \
	--input-dir "$INPUT_DIR" \
	--output-dir "$OUTPUT_DIR" \
	--compare-dir "$COMPARE_DIR" \
	"${extra[@]}" | tee "$COMPARE_DIR/run.log"

echo "完成: $COMPARE_DIR"
ls -la "$COMPARE_DIR"
