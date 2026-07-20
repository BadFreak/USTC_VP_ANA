# -*- coding: utf-8 -*-
"""断点模式：dots —— 散点不连线，无跨断点误导
图存到 result/ ，文件名带 _dots 后缀。"""

import telemetry_lib as t

t.GAP_MODE = "dots" 
t.run_all(show=True)
print("dots 模式完成。")
