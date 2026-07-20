#!/bin/bash
# 对一个 .dat：在 result/<文件名>/ 下生成 cosmic ped amp wave 四个子目录
# package_mode_id: 0=cosmic, 4=ped, 8=amp, 12=wave（与 cosmic_ped / calib 的 yaml 一致）
#
# 用法: bash run_comprehensive.sh [data_microsat/xxx.dat]
#       bash run_comprehensive.sh -f [data_microsat/xxx.dat]   # 强制重跑，不 skip
# 改下面开关，只跑需要的模式
FILE_NAME=VLAST_data_260706_1848
RUN_COSMIC=true
RUN_PED=true
RUN_AMP=true
RUN_WAVE=true
FORCE_RUN=true

FILE_TYPE=".dat"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$MODEL_DIR/.." && pwd)"

VDECODE="$REPO_ROOT/vdecode/build/vdecode"
COSMIC_PED="$MODEL_DIR/cosmic_ped"
CALIB="$MODEL_DIR/calib"
DATA_DIR="$SCRIPT_DIR/data_zhangjiang"
RESULT_DIR="$SCRIPT_DIR/result_zhangjiang"

CFG_DIR="$SCRIPT_DIR/config"
CFG_COSMIC="$CFG_DIR/cosmic.yaml"
CFG_PED="$CFG_DIR/ped.yaml"
CFG_AMP="$CFG_DIR/amp.yaml"
CFG_WAVE="$CFG_DIR/wave.yaml"
PLOT_AFTER_DECODE="$SCRIPT_DIR/scripts/plot_calo_csi_after_decode.py"
FIND_LACK_TRIGGER="$SCRIPT_DIR/scripts/find_lack_trigger.cxx"
DECODE_SUMMARY_REPORT="$SCRIPT_DIR/scripts/write_decode_summary_report.py"

write_decode_summary_report() {
	if [[ ! -f "$DECODE_SUMMARY_REPORT" ]]; then
		echo "未找到 $DECODE_SUMMARY_REPORT"
		return 0
	fi
	python3 "$DECODE_SUMMARY_REPORT" "$outroot" "$binbasename" \
		|| echo "[warn] 生成解析报告失败"
}

find_lack_trigger_check() {
	local name="$1"
	local resultfile="$2"
	local packfile="${3:-}"
	if [[ ! -f "$FIND_LACK_TRIGGER" ]]; then
		echo "[$name] 未找到 $FIND_LACK_TRIGGER"
		return 0
	fi
	if [[ -n "$packfile" ]]; then
		time root -l -b -q "$FIND_LACK_TRIGGER(\"$resultfile\",\"$packfile\")" \
			>log/log_find_lack_trigger 2>&1 \
			|| echo "[$name] find_lack_trigger 失败"
	else
		time root -l -b -q "$FIND_LACK_TRIGGER(\"$resultfile\")" \
			>log/log_find_lack_trigger 2>&1 \
			|| echo "[$name] find_lack_trigger 失败"
	fi
}

if [[ ! -x "$VDECODE" ]]; then
	echo "找不到 vdecode: $VDECODE （请先编译 vdecode）"
	exit 1
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
		-f|--force)
			FORCE_RUN=true
			shift
			;;
		-h|--help)
			echo "用法: bash run_comprehensive.sh [-f|--force] [data.pkg]"
			echo "  -f, --force  忽略已完成检查，强制重跑"
			exit 0
			;;
		-*)
			echo "未知选项: $1"
			exit 1
			;;
		*)
			break
			;;
	esac
done

binfile="${1:-$DATA_DIR/$FILE_NAME$FILE_TYPE}"
if [[ ! -f "$binfile" ]]; then
	echo "文件不存在: $binfile"
	exit 1
fi

if [[ "$binfile" == *.* ]]; then
	binext=".${binfile##*.}"
else
	binext="$FILE_TYPE"
fi
binbasename="$(basename "$binfile" "$binext")"
outroot="$RESULT_DIR/$binbasename"
mkdir -p "$outroot"

MODES=(cosmic ped amp wave)

mode_done() {
	local name="$1"
	local work="$outroot/$name"
	[[ -d "$work" ]] && compgen -G "$work/result_*.root" >/dev/null
}

run_cosmic_ped() {
	local name="$1"
	local cfg="$2"
	local work="$outroot/$name"
	mkdir -p "$work/log"
	cd "$work"

	echo "======== $name ========"
	cp -f "$cfg" config.yaml

	time "$VDECODE" combine "$binfile" config.yaml >log/log_decode 2>&1
	local resultfile
	resultfile="$(ls result_*.root 2>/dev/null | head -1)"
	if [[ -z "$resultfile" ]]; then
		echo "[$name] decode 未生成 result_*.root"
		return 1
	fi
	cp -f "$resultfile" result.root

	if [[ -f "$PLOT_AFTER_DECODE" ]]; then
		python3 "$PLOT_AFTER_DECODE" "$resultfile" "$work" >log/log_plot_decode_check 2>&1 \
			|| echo "[$name] plot_calo_csi_after_decode 失败"
	else
		echo "[$name] 未找到 $PLOT_AFTER_DECODE"
	fi

	time root -l -b -q "$COSMIC_PED/scripts/pack_calo_by_trigger.cxx(\"$resultfile\",\"pack_$resultfile\")" \
		>log/log_pack_calo 2>&1
	cp -f "pack_$resultfile" pack.root

	find_lack_trigger_check "$name" "$resultfile" "pack.root"

	# cosmic: 恰好 1 个晶体击中 (ADC-ped>200)；ped: 不筛选
	local n_crystal_hit=-1
	local hit_threshold=200
	if [[ "$name" == "cosmic" ]]; then
		n_crystal_hit=1
	fi
	time root -l -b -q "$COSMIC_PED/scripts/draw_calo.cxx(1,\"pack_$resultfile\",${n_crystal_hit},${hit_threshold},false)" \
		>log/log_draw_calo 2>&1
	time python "$COSMIC_PED/scripts/VPlot_calo.py" >log/log_plot_calo 2>&1

	time root -l -b -q "$COSMIC_PED/scripts/draw_csi.cxx(\"$resultfile\")" >log/log_draw_csi 2>&1
}

run_calib() {
	local name="$1"
	local cfg="$2"
	local work="$outroot/$name"
	mkdir -p "$work/log"
	cd "$work"

	echo "======== $name ========"
	cp -f "$cfg" config.yaml

	time "$VDECODE" combine "$binfile" config.yaml >log/log_decode 2>&1
	local resultfile
	resultfile="$(ls result_*.root 2>/dev/null | head -1)"
	if [[ -z "$resultfile" ]]; then
		echo "[$name] decode 未生成 result_*.root"
		return 1
	fi

	if [[ -f "$PLOT_AFTER_DECODE" ]]; then
		python3 "$PLOT_AFTER_DECODE" "$resultfile" "$work" >log/log_plot_decode_check 2>&1 \
			|| echo "[$name] plot_calo_csi_after_decode 失败"
	else
		echo "[$name] 未找到 $PLOT_AFTER_DECODE"
	fi

	find_lack_trigger_check "$name" "$resultfile"

	time root -l -b -q "$CALIB/scripts/draw_calo.cxx(1,\"${resultfile}\")" >log/log_draw_calo 2>&1
	time root -l -b -q "$CALIB/scripts/calib_calo.cxx()" >log/log_calib_calo 2>&1
	time python "$CALIB/scripts/VPlot_calo.py" >log/log_plot_calo 2>&1

	time root -l -b -q "$CALIB/scripts/draw_csi.cxx(1,\"${resultfile}\")" >log/log_draw_csi 2>&1
	time root -l -b -q "$CALIB/scripts/calib_csi.cxx()" >log/log_calib_csi 2>&1
	time python "$CALIB/scripts/VPlot_csi.py" >log/log_plot_csi 2>&1
}

echo "输入: $binfile"
echo "输出根目录: $outroot"
echo "开关: cosmic=$RUN_COSMIC ped=$RUN_PED amp=$RUN_AMP wave=$RUN_WAVE force=$FORCE_RUN"

if $RUN_COSMIC; then
	if ! $FORCE_RUN && mode_done cosmic; then
		echo "[skip] cosmic 已完成"
	else
		$FORCE_RUN && mode_done cosmic && echo "[force] 重新运行 cosmic"
		run_cosmic_ped cosmic "$CFG_COSMIC" || echo "[warn] cosmic 失败"
	fi
fi
if $RUN_PED; then
	if ! $FORCE_RUN && mode_done ped; then
		echo "[skip] ped 已完成"
	else
		$FORCE_RUN && mode_done ped && echo "[force] 重新运行 ped"
		run_cosmic_ped ped "$CFG_PED" || echo "[warn] ped 失败"
	fi
fi
if $RUN_AMP; then
	if ! $FORCE_RUN && mode_done amp; then
		echo "[skip] amp 已完成"
	else
		$FORCE_RUN && mode_done amp && echo "[force] 重新运行 amp"
		run_calib amp "$CFG_AMP" || echo "[warn] amp 失败"
	fi
fi
if $RUN_WAVE; then
	if ! $FORCE_RUN && mode_done wave; then
		echo "[skip] wave 已完成"
	else
		$FORCE_RUN && mode_done wave && echo "[force] 重新运行 wave"
		run_calib wave "$CFG_WAVE" || echo "[warn] wave 失败"
	fi
fi

write_decode_summary_report

echo "完成: $outroot"
