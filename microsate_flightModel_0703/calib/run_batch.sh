#!/bin/bash
FILE_TYPE=".dat"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data_microsat"
RESULT_DIR="$SCRIPT_DIR/result_microsat"
mkdir -p $RESULT_DIR
binfiles=`ls $DATA_DIR/*${FILE_TYPE}`
configfile="$SCRIPT_DIR/../config/triggerIDToBeDrawn.yaml"
logfileDir="$RESULT_DIR/$binfileDir"
mkdir -p logfileDir
for binfile in $binfiles; do
    echo "Processing: $binfile"
    binfileDir=$(basename "$binfile" $FILE_TYPE)
    mkdir -p  $RESULT_DIR/$binfileDir && cd $RESULT_DIR/$binfileDir
    /home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/vdecode/build/vdecode combine "$binfile" "$configfile" > log_decode 2>&1
    resultfile=`ls result_*.root`
    echo "resultfile: $resultfile"
    mkdir -p log
    root -l -b -q $SCRIPT_DIR/scripts/draw_calo.cxx\(1,\"${resultfile}\"\) > log/log_draw_calo 2>&1
    root -l -b -q $SCRIPT_DIR/scripts/calib_calo.cxx\(\) > log/log_calib_calo 2>&1
    python $SCRIPT_DIR/scripts/VPlot_calo.py > log/log_plot_calo 2>&1
    root -l -b -q $SCRIPT_DIR/scripts/draw_csi.cxx\(1,\"${resultfile}\"\) > log/log_draw_csi 2>&1
    root -l -b -q $SCRIPT_DIR/scripts/calib_csi.cxx\(\) > log/log_calib_csi 2>&1
    python $SCRIPT_DIR/scripts/VPlot_csi.py > log/log_plot_csi 2>&1
    cd ../../
done
