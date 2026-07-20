#!/bin/bash
# 对 input/calo 下各 result_*.root：按 comprehensive cosmic 流程
#   pack_calo_by_trigger -> draw_calo(1-hit, thr=200) -> 中间晶体 HG main/backup MIP 对比图
# 只在本目录写输出，不改外部脚本。
#
# 用法:
#   bash run_compare_center_mip.sh
#   bash run_compare_center_mip.sh -f          # 强制重跑 pack / hist
#   bash run_compare_center_mip.sh --plot-only # 仅画图（需已有 hist_calo.root）

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # microsate_flightModel_0703
COSMIC_PED="$MODEL_DIR/cosmic_ped"
PACK_CXX="$COSMIC_PED/scripts/pack_calo_by_trigger.cxx"
DRAW_CXX="$COSMIC_PED/scripts/draw_calo.cxx"
PLOT_PY="$SCRIPT_DIR/plot_center_mip_compare.py"

INPUT_DIR="$SCRIPT_DIR/input/calo"
OUTPUT_DIR="$SCRIPT_DIR/output/calo"
COMPARE_DIR="$OUTPUT_DIR/compare"
CENTER_CRYSTAL=12
N_CRYSTAL_HIT=1
HIT_THRESHOLD=200

FORCE=false
PLOT_ONLY=false

while [[ $# -gt 0 ]]; do
	case "$1" in
		-f|--force) FORCE=true; shift ;;
		--plot-only) PLOT_ONLY=true; shift ;;
		-h|--help)
			echo "用法: bash run_compare_center_mip.sh [-f|--force] [--plot-only]"
			exit 0
			;;
		*) echo "未知选项: $1"; exit 1 ;;
	esac
done

# 测试顺序固定，legend / pad 标号与此一致
TESTS=(firstFSE secondFSE specialItem vibration onlyCalo)

pack_entries() {
	local packfile="$1"
	python3 -c "import ROOT; f=ROOT.TFile.Open('$packfile'); t=f.Get('caloTree'); print(t.GetEntries() if t else 0); f.Close()" 2>/dev/null || echo 0
}

is_already_packed_result() {
	# 启发式：前若干 entry 的 TriggerID 几乎不重复且 Cell 数接近 100
	local resultfile="$1"
	python3 - "$resultfile" <<'PY'
import sys
import ROOT
path = sys.argv[1]
f = ROOT.TFile.Open(path)
if not f or f.IsZombie():
    sys.exit(1)
t = f.Get("caloTree")
if not t:
    sys.exit(1)
n = min(40, int(t.GetEntries()))
if n < 4:
    sys.exit(1)
tids = []
nc = []
for i in range(n):
    t.GetEntry(i)
    tids.append(int(t.TriggerID))
    nc.append(len(t.CellID))
f.Close()
unique = len(set(tids))
avg_nc = sum(nc) / float(len(nc))
# 已 pack：Trigger 几乎唯一，且 cell 数接近 4*25=100
ok = unique >= n - 2 and avg_nc >= 80
sys.exit(0 if ok else 1)
PY
}

process_one() {
	local name="$1"
	local resultfile="$INPUT_DIR/result_${name}.root"
	local work="$OUTPUT_DIR/$name"
	mkdir -p "$work/log"

	if [[ ! -f "$resultfile" ]]; then
		echo "[$name] 缺少输入: $resultfile"
		return 1
	fi

	# 损坏 / 截断的 ROOT 提前报错
	if ! python3 - "$resultfile" <<'PY'
import sys
import ROOT
f = ROOT.TFile.Open(sys.argv[1])
ok = f and (not f.IsZombie()) and f.Get("caloTree") is not None
if f:
    f.Close()
sys.exit(0 if ok else 1)
PY
	then
		echo "[$name] 输入 ROOT 无法打开或无 caloTree（可能截断）: $resultfile"
		return 1
	fi

	local packfile="$work/pack.root"
	local histfile="$work/hist_calo.root"

	if ! $PLOT_ONLY; then
		local npack=0
		if [[ -f "$packfile" ]]; then
			npack="$(pack_entries "$packfile")"
		fi
		# 空 pack（例如对已 packed 的 vibration/onlyCalo 误跑 pack）视为无效，强制重建
		if $FORCE || [[ ! -f "$packfile" ]] || [[ "${npack:-0}" -le 0 ]]; then
			if is_already_packed_result "$resultfile"; then
				echo "[$name] 输入已是 packed 形态，复制为 pack.root"
				cp -f "$resultfile" "$packfile"
			else
				echo "[$name] pack_calo_by_trigger ..."
				(
					cd "$work"
					root -l -b -q "$PACK_CXX(\"$resultfile\",\"pack.root\")" \
						>log/log_pack_calo 2>&1
				)
			fi
		else
			echo "[$name] 复用已有 pack.root (entries=$npack)"
		fi

		if [[ ! -f "$packfile" ]]; then
			echo "[$name] 未生成 pack.root"
			return 1
		fi

		npack="$(pack_entries "$packfile")"
		if [[ "${npack:-0}" -le 0 ]]; then
			echo "[$name] pack.root 无事例 (entries=$npack)"
			return 1
		fi
		echo "[$name] pack entries=$npack"

		if $FORCE || [[ ! -f "$histfile" ]]; then
			echo "[$name] draw_calo 1-hit thr=$HIT_THRESHOLD ..."
			(
				cd "$work"
				# 与 run_comprehensive.sh cosmic 一致: draw_calo(1,"pack_...",1,200,false)
				root -l -b -q "$DRAW_CXX(1,\"pack.root\",${N_CRYSTAL_HIT},${HIT_THRESHOLD},false)" \
					>log/log_draw_calo 2>&1
			)
		else
			echo "[$name] 复用已有 hist_calo.root"
		fi
	fi

	if [[ ! -f "$histfile" ]]; then
		echo "[$name] 缺少 hist_calo.root"
		return 1
	fi
	echo "[$name] OK  hist=$histfile"
}

echo "INPUT : $INPUT_DIR"
echo "OUTPUT: $OUTPUT_DIR"
echo "center crystal index: $CENTER_CRYSTAL (hmh_${CENTER_CRYSTAL} / hbh_${CENTER_CRYSTAL})"
echo "filter: nCrystalHit==$N_CRYSTAL_HIT, threshold=$HIT_THRESHOLD"

ok_tests=()
for name in "${TESTS[@]}"; do
	if process_one "$name"; then
		ok_tests+=("$name")
	else
		echo "[warn] $name 跳过"
	fi
done

if [[ ${#ok_tests[@]} -eq 0 ]]; then
	echo "没有任何可用测试，退出"
	exit 1
fi

mkdir -p "$COMPARE_DIR"
echo "======== 画对比图 (${#ok_tests[@]} tests) ========"
python3 "$PLOT_PY" \
	--tests "${ok_tests[@]}" \
	--output-dir "$OUTPUT_DIR" \
	--compare-dir "$COMPARE_DIR" \
	--crystal "$CENTER_CRYSTAL"

echo "完成: $COMPARE_DIR"
ls -la "$COMPARE_DIR"
