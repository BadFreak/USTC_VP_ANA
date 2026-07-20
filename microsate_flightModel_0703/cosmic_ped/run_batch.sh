#!/bin/bash
# Process files in the data directory
# Change FILE_TYPE to ".bin" or ".dat" to select file type

# ===== Configuration: Change this to select file type =====
# FILE_TYPE=".pkg"  # Options: ".bin" or ".dat"
# VPlot.py 额外参数：关闭 hmh/hbh/hml/hbl 上 Landau 拟合与 MPV 时设为 --no-sig-fit
VPLOT_ARGS1=""
VPLOT_ARGS2=""
FILE_TYPE=".dat"  # Options: ".bin" or ".dat"
anaCalo=true
anaCsI=true
anaHK=false
# ===========================================================

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data_microsat"
RESULT_DIR="$SCRIPT_DIR/result_microsat"
# DATA_DIR="$SCRIPT_DIR/data_test"
# RESULT_DIR="$SCRIPT_DIR/result_test"
mkdir -p $RESULT_DIR

# Find files of the specified type in data directory
configfile="$SCRIPT_DIR/../config/triggerIDToBeDrawn.yaml"
binfiles=$(find "$DATA_DIR" -maxdepth 1 -name "*${FILE_TYPE}" -type f | sort)

# Process each file
for binfile in $binfiles; do
    # Extract filename without path and extension for result directory name
    binbasename=$(basename "$binfile" "$FILE_TYPE")
    resultdir="$RESULT_DIR/$binbasename"
    
    # 若结果目录已存在且已有 result_*.root，则跳过 decode，只执行 draw 和 plot
    skip_decode=0
	if false; then
    # if [ -d "$resultdir" ]; then
        existing_result=$(find "$resultdir" -maxdepth 1 -name "result_*.root" -type f 2>/dev/null | head -1)
        if [ -n "$existing_result" ]; then
            skip_decode=1
        fi
    fi

    echo "=========================================="
    echo "Processing: $binbasename"
    echo "Output directory: $resultdir"
    [ "$skip_decode" -eq 1 ] && echo "(skipping decode, re-running draw & plot only)"
    echo "=========================================="

    mkdir -p "$resultdir"
    cd "$resultdir" || exit 1

    if [ "$skip_decode" -eq 0 ]; then
        # Step 1: decode（仅当未跳过时执行）
        echo "[Step 1] decode: running..."
        time (/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/vdecode/build/vdecode combine "$binfile" "$configfile" > log_decode 2>&1)
        resultfile=$(ls result_*.root 2>/dev/null | head -1)
        cp $resultfile result.root
        if [ -z "$resultfile" ]; then
            echo "Warning: No result_*.root file generated for $binbasename, skipping..."
            cd "$SCRIPT_DIR"
            continue
        fi
        echo "Decoded file: $resultfile"
    else
        resultfile=$(ls result_*.root 2>/dev/null | head -1)
        if [ -z "$resultfile" ]; then
            echo "Warning: skip_decode=1 but no result_*.root in $resultdir, skipping..."
            cd "$SCRIPT_DIR"
            continue
        fi
        echo "Using existing file: $resultfile"
    fi
    # echo $resultfile

    # Step 2: draw（每次执行）
   
    if $anaCalo; then
        echo "[Step 2a] Calo..."
        time (root -l -b -q "../../scripts/pack_calo_by_trigger.cxx(\"$resultfile\",\"pack_$resultfile\")" > log_pack_calo 2>&1)
        cp pack_$resultfile pack.root 
        time (root -l -b -q "../../scripts/draw_calo.cxx(1,\"pack_$resultfile\",1,200,true)" > log_draw_calo 2>&1)
        # time (python ../../VPlot_calo.py $VPLOT_ARGS > log_plot_calo 2>&1)
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


    # # Step 3: plot（每次执行）
    # echo "[Step 3] plot: running..."
    # time (python ../../VPlot.py $VPLOT_ARGS > log_plot 2>&1)
    
    echo "Completed processing: $binbasename"
    
    # Return to script directory for next iteration
    cd "$SCRIPT_DIR" || exit 1
done

echo "All files processed successfully!"

