#!/bin/bash
binfiles=`ls data/*.pkg`
for binfile in $binfiles; do
    binfile=$(realpath "$binfile")
    binfileDir=$(basename "$binfile" .pkg)
    mkdir -p  result/$binfileDir && cd result/$binfileDir
    ../../../../../../build/vdecode combine $binfile > log_decode 2>&1
    resultfile=`ls result_*.root`
    echo "resultfile: $resultfile"
    root -l -b -q ../../draw.cxx\(1,\"${resultfile}\"\) > log_draw 2>&1
    # peaks.root（供 VPlot plot_calib）由 hist.root 经 calib.cxx 生成
    root -l -b -q ../../calib.cxx\(\) > log_calib 2>&1
    python ../../VPlot.py > log_plot 2>&1
    cd ../../
done
