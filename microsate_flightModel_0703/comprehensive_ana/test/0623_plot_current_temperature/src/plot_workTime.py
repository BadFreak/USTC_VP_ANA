# -*- coding: utf-8 -*-
"""工作时段模式：workTime —— 保留断点，关机段留空
图存到 result/ ，文件名带 _workTime 后缀。"""

import telemetry_lib as t

t.GAP_MODE = "workTime"
t.run_all(show=True)
print("workTime 模式完成。")
