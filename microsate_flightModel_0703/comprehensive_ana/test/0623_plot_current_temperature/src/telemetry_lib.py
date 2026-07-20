# -*- coding: utf-8 -*-
"""
遥测画图 公共库（不直接运行，被 plot_workTime/plot_dots/plot_connect.py 调用）
放在 0623/src/ 下。依赖: pip install pandas matplotlib openpyxl
"""

from pathlib import Path
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------- 路径 ----------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RESULT_DIR = PROJECT_DIR / "result"
RESULT_DIR.mkdir(exist_ok=True)

# ---------- 中文字体（自动注册本机可用字体） ----------
def _setup_cjk_font():
    candidates = [
        "/usr/share/fonts/google-noto-cjk/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/google-noto-cjk/NotoSansCJK-DemiLight.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/google-droid-sans-fonts/DroidSansFallbackFull.ttf",
        "/System/Library/Fonts/PingFang.ttc",
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/simhei.ttf",
    ]
    for path in candidates:
        if not Path(path).exists():
            continue
        try:
            fm.fontManager.addfont(path)
            name = fm.FontProperties(fname=path).get_name()
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            return name
        except (OSError, ValueError):
            continue
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    return None


_setup_cjk_font()

TIME_COL = "地面时间"

# ==================== 可调参数 ====================
CSI_CURRENT_VALID = (50, 440)    # 径迹 CsI 电流有效量程 / 纵轴
CALO_CURRENT_VALID = (50, 1400)   # 量能器电流有效量程 / 纵轴
TEMP_VALID = (15, 45)             # 温度有效量程(℃) / 纵轴
TEMP_YLIM = TEMP_VALID
GAP_MODE = "workTime"         # 断点方式，由入口文件设置: workTime / dots / connect
# 量能器电流：((cur[11:0]*2500/4096 - 245.1)/20)/divisor，单位 mA
FEE_CUR0_DIVISOR = 0.01    # 监测点1 (cur0)，FEE1~4
FEE_CUR1_DIVISOR = 0.125   # 监测点2 (cur1)，FEE1~4
# FEE1 温度：cur[11:0] 
# t = NUM / (BETA + sqrt(BETA^2 + GAMMA*(DELTA - ln(6408.7*tempx/(2500-0.64087*tempx))))) - 273.15
FEE1_TEMP_NUM = -172843.44828
FEE1_TEMP_BETA = -4622.53337
FEE1_TEMP_GAMMA = 345686.89656
FEE1_TEMP_DELTA = -6.001188
FEE1_TEMP_LN_NUM = 6408.7
FEE1_TEMP_LN_DEN_OFFSET = 2500.0
FEE1_TEMP_LN_DEN_SCALE = 0.64087
FEE1_TEMP_KELVIN_OFFSET = 273.15
# =================================================


def load_data(filename):
    path = DATA_DIR / filename
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, encoding="utf-8-sig")
    if TIME_COL in df.columns:
        df[TIME_COL] = pd.to_datetime(df[TIME_COL], errors="coerce")
    return df


def pick(df, *candidates):
    """从多个候选列名里挑第一个存在的（兼容不同包列名写法不一致）。"""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def clean_signal(series, valid=None, low12=False):
    """滤无效值：溢出占位(>=1e6)一律滤；若给了 valid=(lo,hi) 则范围外也滤成 NaN。
    low12=True 时先取 16 位值的低 12 位 [11:0]（量能器 FEE1 电流/温度）。"""
    s = pd.to_numeric(series, errors="coerce")
    if low12:
        s = (s.astype("Int64") & 0x0FFF).astype(float)
    s = s.mask(s >= 1e6)
    if valid is not None:
        lo, hi = valid
        s = s.mask((s < lo) | (s > hi))
    return s


def fee_cur12(series, low12=False):
    """量能器 raw → cur[11:0]。"""
    s = pd.to_numeric(series, errors="coerce")
    s = s.mask(s >= 1e6)
    if low12:
        s = (s.astype("Int64") & 0x0FFF).astype(float)
    return s


def fee_current_ma(series, point_idx, low12=False):
    """监测点1(cur0)/2(cur1) → mA。point_idx: 0 或 1。"""
    cur = fee_cur12(series, low12=low12)
    divisor = FEE_CUR0_DIVISOR if point_idx == 0 else FEE_CUR1_DIVISOR
    ma = ((cur * 2500.0 / 4096.0 - 245.1) / 20.0) / divisor
    lo, hi = CALO_CURRENT_VALID
    return ma.mask((ma < lo) | (ma > hi))


def fee1_temp_c(series):
    """FEE1 温度：raw 取低 12 位 [11:0]，tempx=cur*2500/4096，再按文档公式换算为 ℃。"""
    tempx = fee_cur12(series, low12=True)
    denom_inner = FEE1_TEMP_LN_DEN_OFFSET - FEE1_TEMP_LN_DEN_SCALE * tempx
    ln_arg = (FEE1_TEMP_LN_NUM * tempx) / denom_inner
    ok = tempx.notna() & (denom_inner > 0) & (ln_arg > 0)
    inner = FEE1_TEMP_DELTA - np.log(ln_arg.where(ok))
    radicand = FEE1_TEMP_BETA ** 2 + FEE1_TEMP_GAMMA * inner
    ok = ok & (radicand >= 0)
    denom = FEE1_TEMP_BETA + np.sqrt(radicand.where(ok))
    ok = ok & (denom != 0)
    temp_c = (FEE1_TEMP_NUM / denom) - FEE1_TEMP_KELVIN_OFFSET
    temp_c = temp_c.where(ok)
    lo, hi = TEMP_VALID
    return temp_c.mask((temp_c < lo) | (temp_c > hi))


def apply_time_axis(ax):
    """横轴两行：上行日期，下行时:分（含分钟）。"""
    locator = mdates.AutoDateLocator(minticks=5, maxticks=12)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    ax.tick_params(axis="x", labelsize=8, pad=3)
    for label in ax.get_xticklabels():
        label.set_ha("center")


def plot_series_on(ax, df, col_map, title, ylabel, valid=None, legend_ncol=1, low12=False, ylim=None):
    """在子图上画多条曲线。col_map={图例名: 列名}；valid=有效量程；
    legend_ncol=图例列数(通道多时设大一点, 如探测器温度18路)。
    画法由全局 GAP_MODE 决定：
      workTime -> 保留 NaN，线自动断开（关机段留空）
      dots    -> 只画点不连线，无跨断点误导
      connect -> 扔掉无效点，画连续线（会横跨关机段）
    """
    x = df[TIME_COL]
    for label, col in col_map.items():
        if col is None:
            continue
        y = clean_signal(df[col], valid, low12=low12)
        if GAP_MODE == "dots":
            ax.plot(x, y, label=label, linestyle="none", marker=".", markersize=1)
        elif GAP_MODE == "connect":
            m = y.notna()
            ax.plot(x[m], y[m], label=label, linewidth=0.8)
        elif GAP_MODE == "workTime":
            ax.plot(x, y, label=label, linewidth=0.8)
    ax.set_title(title, fontsize=12)
    ax.set_ylabel(ylabel)
    fs = 8 if legend_ncol == 1 else 6
    ax.legend(fontsize=fs, loc="upper right", ncol=legend_ncol)
    ax.grid(True, alpha=0.3)
    if ylim is not None:
        ax.set_ylim(*ylim)


def plot_value_series_on(ax, x, value_map, title, ylabel, legend_ncol=1, ylim=None):
    """value_map={图例名: 已换算 Series}。"""
    for label, y in value_map.items():
        if y is None:
            continue
        if GAP_MODE == "dots":
            ax.plot(x, y, label=label, linestyle="none", marker=".", markersize=1)
        elif GAP_MODE == "connect":
            m = y.notna()
            ax.plot(x[m], y[m], label=label, linewidth=0.8)
        elif GAP_MODE == "workTime":
            ax.plot(x, y, label=label, linewidth=0.8)
    ax.set_title(title, fontsize=12)
    ax.set_ylabel(ylabel)
    fs = 8 if legend_ncol == 1 else 6
    ax.legend(fontsize=fs, loc="upper right", ncol=legend_ncol)
    ax.grid(True, alpha=0.3)
    if ylim is not None:
        ax.set_ylim(*ylim)


# ========== 遥测包3 径迹CsI ==========
def plot_csi():
    df = load_data("kx12_1_VLAST探路者遥测包3_0302.csv")
    currents = {
        "电流监测点1": pick(df, "径迹CsI电源电流监测点1-[11:0]"),
        "电流监测点2": pick(df, "径迹CsI电源电流监测点2-[11:0]"),
        "电流监测点3": pick(df, "径迹CsI电源电流监测点3-[11:0]"),
    }
    temps = {
        "温度点1": pick(df, "径迹CsI-FEE温度点1-[11:0]"),
        "温度点2": pick(df, "径迹CsI-FEE温度点2-[11:0]"),
        "温度点3": pick(df, "径迹CsI-FEE温度点3-[11:0]"),
        "温度点4": pick(df, "径迹CsI-FEE温度点4-[11:0]"),
    }
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    plot_series_on(ax1, df, currents, "径迹CsI 电源电流", "电流",
                   valid=CSI_CURRENT_VALID, ylim=CSI_CURRENT_VALID)
    plot_series_on(ax2, df, temps, "径迹CsI-FEE 温度", "温度 (℃)", valid=TEMP_VALID, ylim=TEMP_YLIM)
    ax2.set_xlabel(TIME_COL)
    fig.suptitle(f"遥测包3 · 径迹CsI  [{GAP_MODE}]", fontsize=15, y=0.98)
    apply_time_axis(ax2)
    fig.subplots_adjust(top=0.92, hspace=0.28, bottom=0.12)
    out = RESULT_DIR / f"遥测包3_径迹CsI_{GAP_MODE}.png"
    fig.savefig(out, dpi=150)
    print("已保存:", out.name)


# ========== 遥测包4 量能器 FEE1~4 ==========
def plot_calorimeter():
    df = load_data("kx12_1_VLAST探路者遥测包4_0304.xlsx")
    fig, axes = plt.subplots(4, 2, figsize=(15, 19), sharex=True)
    for i, n in enumerate([1, 2, 3, 4]):
        if n == 1:
            currents = {
                "电流点1": pick(df, "量能器FEE1-电源电流监测点1"),
                "电流点2": pick(df, "量能器FEE1-电源电流监测点2"),
            }
        else:
            currents = {
                "电流点1": pick(df, f"量能器FEE{n}-电源电流监测点1-data[11:0]", f"量能器FEE{n}-电源电流监测点1"),
                "电流点2": pick(df, f"量能器FEE{n}-电源电流监测点2-data[11:0]", f"量能器FEE{n}-电源电流监测点2"),
            }
        if n == 1:
            temps = {f"温度点{k}": pick(df, f"量能器FEE1-温度监测点{k}") for k in [1, 2, 3, 4]}
        else:
            temps = {
                f"温度点{k}": pick(df, f"量能器FEE{n}-温度监测点{k}-data[11:0]", f"量能器FEE{n}-温度监测点{k}")
                for k in [1, 2, 3, 4]
            }
        if n == 1:
            x = df[TIME_COL]
            cur_ma = {}
            if currents["电流点1"]:
                cur_ma["电流点1"] = fee_current_ma(df[currents["电流点1"]], 0, low12=True)
            if currents["电流点2"]:
                cur_ma["电流点2"] = fee_current_ma(df[currents["电流点2"]], 1, low12=True)
            plot_value_series_on(axes[i, 0], x, cur_ma, f"FEE{n} 电源电流", "电流 (mA)", ylim=CALO_CURRENT_VALID)
            temp_c = {}
            for label, col in temps.items():
                if col:
                    temp_c[label] = fee1_temp_c(df[col])
            plot_value_series_on(axes[i, 1], x, temp_c, f"FEE{n} 温度", "温度 (℃)", ylim=TEMP_YLIM)
        else:
            plot_series_on(axes[i, 1], df, temps, f"FEE{n} 温度", "温度 (℃)", valid=TEMP_VALID, ylim=TEMP_YLIM)
            plot_series_on(axes[i, 0], df, currents, f"FEE{n} 电源电流", "电流",
                           valid=CALO_CURRENT_VALID, ylim=CALO_CURRENT_VALID)
    axes[-1, 0].set_xlabel(TIME_COL); axes[-1, 1].set_xlabel(TIME_COL)
    fig.suptitle(f"遥测包4 · 量能器 FEE1~4  [{GAP_MODE}]", fontsize=16, y=0.995)
    apply_time_axis(axes[-1, 0])
    fig.subplots_adjust(top=0.95, hspace=0.4, bottom=0.10)
    out = RESULT_DIR / f"遥测包4_量能器FEE_{GAP_MODE}.png"
    fig.savefig(out, dpi=150)
    print("已保存:", out.name)


def run_all(show=False):
    plot_csi()
    plot_calorimeter()
    if show:
        plt.show()
