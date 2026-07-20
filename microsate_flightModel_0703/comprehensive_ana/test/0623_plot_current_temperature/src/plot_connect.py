# -*- coding: utf-8 -*-
"""断点模式：connect —— 扔掉无效点画连续线
图存到 result/ ，文件名带 _connect 后缀。"""

import telemetry_lib as t

t.GAP_MODE = "connect"
t.run_all(show=True)
print("connect 模式完成。")
