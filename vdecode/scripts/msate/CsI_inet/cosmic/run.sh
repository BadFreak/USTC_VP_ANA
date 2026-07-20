#!/bin/bash
binfile=`ls *.bin`
vdecode amp $binfile > log_decode 2>&1
resultfile=`ls result_*.root`
root -l -b -q draw.cxx\(1,\"$resultfile\"\) > log_draw 2>&1
python VPlot.py > log_plot 2>&1
