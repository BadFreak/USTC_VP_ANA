#!/bin/bash
# 单文件跑 decode + draw；请用: bash run_one.sh

FILE_TYPE=".dat"
VPLOT_ARGS1=" " # --no-ped-fit 
VPLOT_ARGS2=" " # --no-sig-fit
# VPLOT_ARGS1="--no-ped-fit"
# VPLOT_ARGS2="--no-sig-fit"
anaCalo=true
anaCsI=true
anaHK=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data_microsat"
RESULT_DIR="$SCRIPT_DIR/result_microsat"
mkdir -p "$RESULT_DIR"

# 改这里：data_test 下的文件名，或任意绝对路径
# binfile="$DATA_DIR/CsI+量能器_CHAN_A(2026年04月27日19时39分39秒)$FILE_TYPE"
binfile="$DATA_DIR/VLAST_data_260528_05262200$FILE_TYPE"
configfile="$SCRIPT_DIR/../config/triggerIDToBeDrawn.yaml"
echo "binfile: $binfile"
# binfile="$DATA_DIR/calo_1$FILE_TYPE"

binbasename=$(basename "$binfile" "$FILE_TYPE")
resultdir="$RESULT_DIR/$binbasename"
mkdir -p "$resultdir"
cd "$resultdir" || exit 1

echo "[Step 1] decode..."
time (/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/vdecode/build/vdecode combine "$binfile" "$configfile" > log_decode 2>&1)
resultfile=$(ls result_*.root 2>/dev/null | head -1)
cp $resultfile result.root

PLOT_DECODE="$SCRIPT_DIR/../scripts/plot_calo_csi_after_decode.py"
if [[ -f "$PLOT_DECODE" && -n "$resultfile" ]]; then
	echo "[Step 1b] caloTree check plots..."
	python3 "$PLOT_DECODE" "$resultfile" "." > log_plot_decode_check 2>&1 \
		|| echo "plot_calo_csi_after_decode 失败"
fi

if $anaCalo; then
	echo "[Step 2a] Calo..."
	time (root -l -b -q "../../scripts/pack_calo_by_trigger.cxx(\"$resultfile\",\"pack_$resultfile\")" > log_pack_calo 2>&1)
	cp pack_$resultfile pack.root 
	time (root -l -b -q "../../scripts/draw_calo.cxx(1,\"pack_$resultfile\",1,200,false)" > log_draw_calo 2>&1)
	# time (python ../../scripts/VPlot_calo.py $VPLOT_ARGS > log_plot_calo 2>&1)
	# time (root -l -b -q "../../scripts/draw_calo.cxx(1,\"$resultfile\",1,200,true)" > log_draw_calo 2>&1)
	echo "resultfile: $resultfile"
	time (python ../../scripts/VPlot_calo.py $VPLOT_ARGS1 $VPLOT_ARGS2 > log_plot_calo 2>&1)
	# time (root -l -b -q "../../scripts/pack_calo_by_trigger.cxx(\"$resultfile\",\"pack_$resultfile\")" > log_pack_calo 2>&1)
	# time (root -l -b -q "../../scripts/draw_new.cxx(1,\"pack_$resultfile\",1,300,true)" > log_draw_calo 2>&1)
    # time (python ../../scripts/VPlot.py $VPLOT_ARGS > log_plot 2>&1)
fi
if $anaCsI; then
	echo "[Step 2b] CsI..."
	time (root -l -b -q "../../scripts/draw_csi.cxx(\"$resultfile\")" > log_draw_csi 2>&1)
fi
if $anaHK; then
	echo "[Step 2c] HK..."
	time (root -l -b -q "../../scripts/draw_hk.cxx(\"$resultfile\")" > log_draw_hk 2>&1)
fi

cd "$SCRIPT_DIR"
echo "Done: $resultdir"
