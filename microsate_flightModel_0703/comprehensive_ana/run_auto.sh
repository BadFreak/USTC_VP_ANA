#!/bin/bash
# 每 10 秒检查 data_microsat 是否有新增/变更；有则后台跑 run_comprehensive_batch.sh
# batch 运行期间继续探测；若发现新文件则等当前 batch 结束后立即启动下一轮
#
# 用法:
#   bash run_auto.sh              # 前台循环
#   bash run_auto.sh .pkg         # 只监视 .pkg
#   nohup bash run_auto.sh &      # 后台运行

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data_microsat"
BATCH="$SCRIPT_DIR/run_comprehensive_batch.sh"
STATE_FILE="$SCRIPT_DIR/.data_microsat_watch_state"
LOCK_FILE="$SCRIPT_DIR/.run_auto.lock"
RESULT_DIR="$SCRIPT_DIR/result_microsat"
LOG_DIR="$SCRIPT_DIR/logs"
LOG_FILE="$LOG_DIR/run_auto.log"

INTERVAL_SEC=10
BATCH_PID=""
BATCH_JUST_FINISHED=false
FOLLOWUP_AFTER_BATCH=false

ensure_log_dir() {
	if [[ -e "$LOG_DIR" && ! -d "$LOG_DIR" ]]; then
		echo "错误: $LOG_DIR 已存在且不是目录，请改名后重试" >&2
		exit 1
	fi
	mkdir -p "$LOG_DIR"
}

if [[ -n "$1" ]]; then
	FILE_EXTS=("$1")
	BATCH_ARGS=("$1")
else
	FILE_EXTS=(".pkg" ".dat" ".bin")
	BATCH_ARGS=()
fi

mkdir -p "$DATA_DIR" "$RESULT_DIR"
ensure_log_dir

exec 8>"$LOCK_FILE"
if ! flock -n 8; then
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] 已有 run_auto 在运行，退出" | tee -a "$LOG_FILE"
	exit 1
fi

log() {
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

snapshot_data_dir() {
	local f ext
	for ext in "${FILE_EXTS[@]}"; do
		for f in "$DATA_DIR"/*"$ext"; do
			[[ -f "$f" ]] || continue
			stat -c '%n %s %Y' "$f" 2>/dev/null || stat -f '%N %z %m' "$f"
		done
	done | sort
}

batch_running() {
	[[ -n "$BATCH_PID" ]] && kill -0 "$BATCH_PID" 2>/dev/null
}

has_result_dir() {
	local name="$1"
	local ext="$2"
	[[ -d "$RESULT_DIR/$name" ]] && return 0
	if [[ "$ext" == ".dat" && -d "$RESULT_DIR/${name}.dat" ]]; then
		return 0
	fi
	return 1
}

has_pending_files() {
	local f ext b
	for ext in "${FILE_EXTS[@]}"; do
		for f in "$DATA_DIR"/*"$ext"; do
			[[ -f "$f" ]] || continue
			b="$(basename "$f" "$ext")"
			has_result_dir "$b" "$ext" || return 0
		done
	done
	return 1
}

reap_batch() {
	BATCH_JUST_FINISHED=false
	if [[ -z "$BATCH_PID" ]]; then
		return 0
	fi
	if kill -0 "$BATCH_PID" 2>/dev/null; then
		return 0
	fi
	if wait "$BATCH_PID"; then
		log "batch 结束 pid=$BATCH_PID (成功)"
	else
		log "batch 结束 pid=$BATCH_PID (失败 exit=$?)"
	fi
	BATCH_PID=""
	BATCH_JUST_FINISHED=true
}

start_batch() {
	if batch_running; then
		return 1
	fi
	bash "$BATCH" "${BATCH_ARGS[@]}" >>"$LOG_FILE" 2>&1 &
	BATCH_PID=$!
	log "已启动 batch pid=$BATCH_PID"
	return 0
}

try_start_batch() {
	local cur="$1"
	if start_batch; then
		printf '%s\n' "$cur" >"$STATE_FILE"
		FOLLOWUP_AFTER_BATCH=false
		return 0
	fi
	return 1
}

if [[ ! -f "$BATCH" ]]; then
	log "错误: 找不到 $BATCH"
	exit 1
fi

log "启动 run_auto: 间隔=${INTERVAL_SEC}s 监视=${DATA_DIR} 后缀=${FILE_EXTS[*]}"

while true; do
	reap_batch
	cur="$(snapshot_data_dir)"
	prev=""
	[[ -f "$STATE_FILE" ]] && prev="$(cat "$STATE_FILE")"

	has_snapshot_change=false
	[[ "$cur" != "$prev" ]] && has_snapshot_change=true

	needs_work=false
	if has_snapshot_change || has_pending_files; then
		needs_work=true
	fi

	if batch_running; then
		if $needs_work; then
			FOLLOWUP_AFTER_BATCH=true
			log "检测到新文件/待处理，batch 运行中 (pid=$BATCH_PID)，待结束后立即处理"
		fi
	elif $BATCH_JUST_FINISHED && $FOLLOWUP_AFTER_BATCH; then
		log "batch 已结束，立即处理积压文件"
		try_start_batch "$cur"
	elif [[ ! -f "$STATE_FILE" ]]; then
		log "首次检查"
		if [[ -z "$cur" ]]; then
			log "data_microsat 暂无数据文件"
			: >"$STATE_FILE"
		elif $needs_work; then
			try_start_batch "$cur"
		fi
	elif $needs_work; then
		if $has_snapshot_change; then
			log "检测到 data_microsat 有更新"
		else
			log "仍有未处理文件，启动 batch"
		fi
		try_start_batch "$cur"
	else
		log "无更新，跳过"
		FOLLOWUP_AFTER_BATCH=false
	fi

	# batch 刚结束且已立即启动下一轮时，不再多等一个周期
	if $BATCH_JUST_FINISHED && batch_running; then
		continue
	fi

	sleep "$INTERVAL_SEC"
done
