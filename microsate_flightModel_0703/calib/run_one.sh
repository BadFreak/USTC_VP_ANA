#!/bin/bash
FILE_TYPE=".dat"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data_microsat"
RESULT_DIR="$SCRIPT_DIR/result_microsat"
mkdir -p $RESULT_DIR
# binfile="$DATA_DIR/CsI+量能器_CHAN_A(2026年04月27日15时46分52秒)$FILE_TYPE"
binfile="$DATA_DIR/VLAST_data_260528_1430$FILE_TYPE"
echo "Processing: $binfile"
binfileDir=$(basename "$binfile" $FILE_TYPE)
configfile="$SCRIPT_DIR/../config/triggerIDToBeDrawn.yaml"
mkdir -p  $RESULT_DIR/$binfileDir && cd $RESULT_DIR/$binfileDir
logfileDir="$RESULT_DIR/$binfileDir"
mkdir -p log
time (/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/vdecode/build/vdecode combine "$binfile" "$configfile" > log/log_decode 2>&1)
resultfile=`ls result_*.root`
echo "resultfile: $resultfile"
PLOT_DECODE="$SCRIPT_DIR/../scripts/plot_calo_csi_after_decode.py"
if [[ -f "$PLOT_DECODE" && -n "$resultfile" ]]; then
	python3 "$PLOT_DECODE" "$resultfile" "." > log/log_plot_decode_check 2>&1 \
		|| echo "plot_calo_csi_after_decode 失败"
fi
# root -l -b -q "../../pack_calo_by_trigger.cxx(\"$resultfile\",\"pack_$resultfile\")" > log_pack_calo 2>&1    
root -l -b -q $SCRIPT_DIR/scripts/draw_calo.cxx\(1,\"${resultfile}\"\) > log/log_draw_calo 2>&1
root -l -b -q $SCRIPT_DIR/scripts/calib_calo.cxx\(\) > log/log_calib_calo 2>&1
python $SCRIPT_DIR/scripts/VPlot_calo.py > log/log_plot_calo 2>&1
root -l -b -q $SCRIPT_DIR/scripts/draw_csi.cxx\(1,\"${resultfile}\"\) > log/log_draw_csi 2>&1
root -l -b -q $SCRIPT_DIR/scripts/calib_csi.cxx\(\) > log/log_calib_csi 2>&1
python $SCRIPT_DIR/scripts/VPlot_csi.py > log/log_plot_csi 2>&1
cd ../../
