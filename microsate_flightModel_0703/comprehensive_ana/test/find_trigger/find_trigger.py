#!/usr/bin/env python3
"""从 Gamma* 列表文件读取 (trigger_id, time_code)，输出低 12 位 trigger 与 time_code。"""

import argparse
import glob
import os
import sys


def extract_low_12bits(trigger_id):
    """从 16 位 trigger_id 中提取低 12 位（0 ~ 4095）。"""
    return trigger_id & 0xFFF


def load_gamma_file(path):
    data = []
    with open(path, encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                trigger_id = int(parts[0])
                time_code = int(parts[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: cannot parse {line!r}") from exc
            data.append((trigger_id, time_code))
    return data


def resolve_gamma_path(pattern_or_path):
    if os.path.isfile(pattern_or_path):
        return pattern_or_path
    matches = sorted(glob.glob(pattern_or_path))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"No Gamma file found for: {pattern_or_path}")
    msg = "\n".join(f"  {m}" for m in matches)
    raise FileNotFoundError(
        f"Multiple Gamma files match {pattern_or_path!r}:\n{msg}\nPlease specify one explicitly."
    )


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="Read Gamma* list and output low-12-bit trigger_id / time_code arrays."
    )
    parser.add_argument(
        "gamma_file",
        nargs="?",
        default=os.path.join(script_dir, "Gamma*"),
        help="Gamma list file path or glob (default: Gamma* beside this script)",
    )
    args = parser.parse_args()

    try:
        path = resolve_gamma_path(args.gamma_file)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1
    data = load_gamma_file(path)
    if not data:
        print(f"Warning: no entries in {path}", file=sys.stderr)
        return 1

    print(f"# source: {path} ({len(data)} entries)", file=sys.stderr)
    print("原始triggerid -> 新的triggerid (低12位)")
    print("2^12 = 4096，取值范围：0 ~ 4095")
    print("-" * 50)

    print("[", end="")
    for original_id, _ in data:
        print(f"{extract_low_12bits(original_id)}", end=", ")
    print("]")

    print("[", end="")
    for _, time_code in data:
        print(f"{time_code}", end=", ")
    print("]")

    for original_id, time_code in data:
        print(f"{extract_low_12bits(original_id):6d}, {time_code}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
