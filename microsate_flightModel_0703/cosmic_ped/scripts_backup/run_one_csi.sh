#!/bin/bash
# 单文件跑 decode + draw；请用: bash run_one.sh

FILE_TYPE=".bin"
VPLOT_ARGS="--no-sig-fit"
anaCsI=true
anaHK=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data_microsat"
RESULT_DIR="$SCRIPT_DIR/result_microsat/csi"
mkdir -p "$RESULT_DIR"

# 改这里：data_test 下的文件名，或任意绝对路径
binfile="$DATA_DIR/CHAN_A(2026年04月26日20时40分46秒)$FILE_TYPE"
echo "binfile: $binfile"
# binfile="$DATA_DIR/calo_1$FILE_TYPE"

binbasename=$(basename "$binfile" "$FILE_TYPE")
resultdir="$RESULT_DIR/$binbasename"
mkdir -p "$resultdir"
cd "$resultdir" || exit 1

echo "[Step 1] decode..."
time (/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/vdecode/build/vdecode combine "$binfile" > log_decode 2>&1)
resultfile=$(ls result_*.root 2>/dev/null | head -1)


if $anaCsI; then
	echo "[Step 2b] CsI..."
	time (root -l -b -q "../../../draw_csi.cxx(\"$resultfile\")" > log_draw_csi 2>&1)
fi
if $anaHK; then
	echo "[Step 2c] HK..."
	time (root -l -b -q "../../../draw_hk.cxx(\"$resultfile\")" > log_draw_hk 2>&1)
fi

cd "$SCRIPT_DIR"
echo "Done: $resultdir"
