#!/bin/bash
# 示例:
#   bash run.sh
#   bash run.sh /path/to/hist_calo.root main_hlr_12 out.png

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HISTFILE="${1:-/home/test/wangjiaxuan/VLAST-P/USTC_VP_ANA/microsate_flightModel_260525/comprehensive_ana/result_microsat/CsI_Calo_2606_all/cosmic/hist_calo.root}"
HISTNAME="${2:-main_hlg_12}"
OUTPNG="${3:-$SCRIPT_DIR/$HISTNAME.png}"

cd "$SCRIPT_DIR"
root -l -b -q "draw_one_ratio.cxx+(\"$HISTFILE\",\"$HISTNAME\",\"$OUTPNG\")"
