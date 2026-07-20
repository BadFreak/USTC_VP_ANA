#!/usr/bin/env python3
"""汇总 cosmic / ped / amp / wave 的 time_range 与 find_lack_trigger 结果。"""
import os
import re
import sys
from datetime import datetime

MODES = ("cosmic", "ped", "amp", "wave")

# 各模式期望解包事例数（准则）
EXPECTED_EVENTS = {
    "ped": 17000,
    "amp": 32000,
    "wave": 32000,
}

_RE_CALO = re.compile(r"\[ComReader\] \[info\] CALO:\s*(\d+)")
_RE_CSI = re.compile(r"\[ComReader\] \[info\] CsITK:\s*(\d+)")
_RE_TOTAL_HEAD = re.compile(r"\[ComReader\] \[info\] Total head:\s*(\d+)")
_RE_ENTRIES = re.compile(r"^# entries=(\d+)\s*$")


def _read_text(path):
    if not path or not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as f:
        return f.read().rstrip()


def _find_time_range_files(mode_dir):
    calo_path = csi_path = combined_path = None
    if not os.path.isdir(mode_dir):
        return calo_path, csi_path, combined_path
    for name in sorted(os.listdir(mode_dir)):
        if not name.endswith("_time_range.txt"):
            continue
        path = os.path.join(mode_dir, name)
        if name.endswith("_calo_time_range.txt"):
            calo_path = path
        elif name.endswith("_csi_time_range.txt"):
            csi_path = path
        else:
            combined_path = path
    return calo_path, csi_path, combined_path


def _parse_kv_lines(text):
    out = {}
    if not text:
        return out
    for line in text.splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1]
    return out


def _time_range_summary(tr):
    if not tr:
        return None
    return f"[{tr.get('beijing_min', '?')} .. {tr.get('beijing_max', '?')}]"


def _parse_gap_count(text):
    if not text:
        return None
    for line in reversed(text.splitlines()):
        if line.startswith("# gap_count="):
            return line.split("=", 1)[1].strip()
    return None


def _parse_entries(text):
    if not text:
        return None
    for line in text.splitlines():
        m = _RE_ENTRIES.match(line)
        if m:
            return int(m.group(1))
    return None


def _parse_first_trigger_status(text):
    if not text:
        return None
    first_tid = None
    first_ok = None
    for line in text.splitlines():
        if line.startswith("first TriggerID:"):
            first_tid = line.split(":", 1)[1].strip()
        elif line.startswith("OKAY:") or line.startswith("ERROR:"):
            first_ok = line
    if first_tid is None:
        return None
    return f"TriggerID={first_tid}" + (f", {first_ok}" if first_ok else "")


def _parse_decode_log(mode_dir):
    """从 log/log_decode 提取解包事例数与校验报错数。"""
    log_path = os.path.join(mode_dir, "log", "log_decode")
    text = _read_text(log_path)
    if text is None:
        return None
    calo = csi = total_head = None
    crc_err = 0
    sum_err = 0
    for line in text.splitlines():
        if "CRC check failed" in line:
            crc_err += 1
        elif "Accumulation summation check failed" in line:
            sum_err += 1
        m = _RE_CALO.search(line)
        if m:
            calo = int(m.group(1))
        m = _RE_CSI.search(line)
        if m:
            csi = int(m.group(1))
        m = _RE_TOTAL_HEAD.search(line)
        if m:
            total_head = int(m.group(1))
    return {
        "calo": calo,
        "csi": csi,
        "total_head": total_head,
        "crc_err": crc_err,
        "sum_err": sum_err,
        "log_path": log_path,
    }


def _fmt_count(v):
    return "?" if v is None else str(v)


def _expected_note(mode, n):
    exp = EXPECTED_EVENTS.get(mode)
    if exp is None or n is None:
        return ""
    if n == exp:
        return f" (期望={exp}, OK)"
    return f" (期望={exp}, 不符)"


def _append_time_range_section(lines, label, path):
    lines.append(f"### time_range / {label}")
    text = _read_text(path)
    if text is None:
        lines.append(f"(无 {label} time_range 文件)")
    else:
        lines.append(f"# file: {os.path.basename(path)}")
        lines.extend(text.splitlines())
    lines.append("")


def write_report(outroot, basename):
    lines = [
        f"# 解析报告: {basename}",
        f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"# 目录: {outroot}",
        "",
        "## 准则",
        "# ped 期望约 17000 事例；amp / wave 期望约 32000 事例",
        "# cosmic 无固定期望事例数",
        "",
        "## 解包与校验概要",
    ]

    for mode in MODES:
        mode_dir = os.path.join(outroot, mode)
        if not os.path.isdir(mode_dir):
            lines.append(f"- {mode}: (未运行)")
            continue

        dec = _parse_decode_log(mode_dir)
        if dec is None:
            # 回退：从 lack_trigger 的 # entries= 读取
            csi_n = _parse_entries(_read_text(os.path.join(mode_dir, "csi_lack_trigger.txt")))
            calo_n = _parse_entries(_read_text(os.path.join(mode_dir, "calo_lack_trigger.txt")))
            n_ref = calo_n if calo_n is not None else csi_n
            lines.append(
                f"- {mode}: 解包事例 calo={_fmt_count(calo_n)} csi={_fmt_count(csi_n)}"
                f"{_expected_note(mode, n_ref)}; CRC报错=(无 log_decode); "
                f"累加和校验报错=(无 log_decode)"
            )
            continue

        n_ref = dec["calo"] if dec["calo"] is not None else dec["csi"]
        lines.append(
            f"- {mode}: 解包事例 CALO={_fmt_count(dec['calo'])} "
            f"CsITK={_fmt_count(dec['csi'])} TotalHead={_fmt_count(dec['total_head'])}"
            f"{_expected_note(mode, n_ref)}; "
            f"CRC报错={dec['crc_err']}; 累加和校验报错={dec['sum_err']}"
        )

    lines.append("")
    lines.append("## 概要")

    for mode in MODES:
        mode_dir = os.path.join(outroot, mode)
        if not os.path.isdir(mode_dir):
            lines.append(f"- {mode}: (未运行)")
            continue

        calo_tr_path, csi_tr_path, combined_tr_path = _find_time_range_files(mode_dir)
        calo_tr = _parse_kv_lines(_read_text(calo_tr_path))
        csi_tr = _parse_kv_lines(_read_text(csi_tr_path))
        combined_tr = _parse_kv_lines(_read_text(combined_tr_path))
        csi_text = _read_text(os.path.join(mode_dir, "csi_lack_trigger.txt"))
        calo_lack_text = _read_text(os.path.join(mode_dir, "calo_lack_trigger.txt"))

        time_parts = []
        calo_sum = _time_range_summary(calo_tr)
        csi_sum = _time_range_summary(csi_tr)
        if calo_sum:
            time_parts.append(f"calo {calo_sum}")
        if csi_sum:
            time_parts.append(f"csi {csi_sum}")
        if not time_parts:
            all_sum = _time_range_summary(combined_tr)
            time_parts.append(f"time {all_sum}" if all_sum else "time (无 time_range 文件)")

        parts = [", ".join(time_parts)]

        csi_gap = _parse_gap_count(csi_text)
        if csi_text is None:
            parts.append("csi (无 lack_trigger)")
        elif "# empty" in csi_text:
            parts.append("csi empty")
        else:
            parts.append(f"csi gap_count={csi_gap if csi_gap is not None else '?'}")

        if calo_lack_text is None:
            parts.append("calo (无 lack_trigger)")
        elif "# empty" in calo_lack_text:
            parts.append("calo empty")
        else:
            calo_gap = _parse_gap_count(calo_lack_text)
            parts.append(f"calo gap_count={calo_gap if calo_gap is not None else '?'}")

        lines.append(f"- {mode}: " + ", ".join(parts))

    lines.append("")
    for mode in MODES:
        mode_dir = os.path.join(outroot, mode)
        lines.extend(["=" * 72, f"## {mode}", ""])

        if not os.path.isdir(mode_dir):
            lines.append("(未运行或目录不存在)")
            lines.append("")
            continue

        calo_tr_path, csi_tr_path, combined_tr_path = _find_time_range_files(mode_dir)
        _append_time_range_section(lines, "calo", calo_tr_path)
        _append_time_range_section(lines, "csi", csi_tr_path)
        if combined_tr_path:
            _append_time_range_section(lines, "combined", combined_tr_path)

        for tag in ("csi", "calo"):
            lack_path = os.path.join(mode_dir, f"{tag}_lack_trigger.txt")
            lines.append(f"### find_lack_trigger / {tag}")
            lack_text = _read_text(lack_path)
            if lack_text is None:
                lines.append(f"(无 {tag}_lack_trigger.txt)")
            else:
                status = _parse_first_trigger_status(lack_text)
                gap = _parse_gap_count(lack_text)
                if status:
                    lines.append(f"# summary: {status}, gap_count={gap}")
                lines.extend(lack_text.splitlines())
            lines.append("")

    report_path = os.path.join(outroot, f"{basename}_decode_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines).rstrip() + "\n")
    print(f"saved {report_path}")
    return report_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"用法: {sys.argv[0]} <outroot> [basename]", file=sys.stderr)
        sys.exit(1)
    root = sys.argv[1]
    name = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(root.rstrip("/"))
    write_report(root, name)
