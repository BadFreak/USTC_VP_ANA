#!/usr/bin/env python3
"""中间晶体高增益 main/backup：1-hit MIP 谱 Landau×Gauss 拟合对比图。

每个通道一张 canvas、单个 pad，叠加各测试归一化谱；
legend 为测试名 + 事例数 N + MPV。
拟合参数与 cosmic_ped/scripts/VPlot_calo.py::plot_sig(hmh/hbh) 一致。
"""
from __future__ import annotations

import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# 与 VPlot_calo.plot_sig(hmh/hbh) 一致
FIT_XLO, FIT_XHI = 5500.0, 30000.0
SV_WIDTH, SV_SIG = 800.0, 400.0
PLLO_MP, PLHI_MP = 5500.0, 20000.0
PLLO_W, PLHI_W = 100.0, 3000.0
PLLO_S, PLHI_S = 50.0, 2000.0

PLOT_XMIN, PLOT_XMAX = 5000.0, 30000.0
# 归一化分母：直方图中 ADC > NORM_ADC_MIN 的事例数（bin 积分）
NORM_ADC_MIN = 5000.0

# 各测试线色
COLOURS = [
    ROOT.kCyan + 2,
    ROOT.kBlue + 1,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kGreen + 2,
    ROOT.kBlack,
    ROOT.kRed + 1,
]


def _langaufun(x, par):
    if par[3] <= 0.0:
        return 0.0
    invsq2pi = 0.3989422804014
    mpshift = -0.22278298
    np = 100.0
    sc = 5.0
    mpc = par[1] - mpshift * par[0]
    xlow = x[0] - sc * par[3]
    xupp = x[0] + sc * par[3]
    step = (xupp - xlow) / np
    sumv = 0.0
    for i in range(1, int(np / 2) + 1):
        xx = xlow + (i - 0.5) * step
        fland = ROOT.TMath.Landau(xx, mpc, par[0]) / par[0]
        sumv += fland * ROOT.TMath.Gaus(x[0], xx, par[3])
        xx = xupp - (i - 0.5) * step
        fland = ROOT.TMath.Landau(xx, mpc, par[0]) / par[0]
        sumv += fland * ROOT.TMath.Gaus(x[0], xx, par[3])
    return par[2] * step * sumv * invsq2pi / par[3]


def fit_landau_gauss(hist, tag):
    """返回 (ffit, mpv) 或 None。"""
    if hist is None or hist.Integral() <= 0:
        return None
    nbins = hist.GetNbinsX()
    blo = max(1, min(hist.GetXaxis().FindBin(FIT_XLO), nbins))
    bhi = max(1, min(hist.GetXaxis().FindBin(FIT_XHI), nbins))
    if blo > bhi:
        blo, bhi = bhi, blo
    integral_fit = hist.Integral(blo, bhi)
    if integral_fit <= 0:
        return None

    imax = blo
    for b in range(blo, bhi + 1):
        if hist.GetBinContent(b) > hist.GetBinContent(imax):
            imax = b
    mp_init = hist.GetBinCenter(imax)
    mp_init = max(PLLO_MP, min(PLHI_MP, mp_init))

    ffit = ROOT.TF1("langau_%s" % tag, _langaufun, FIT_XLO, FIT_XHI, 4)
    ffit.SetParameters(SV_WIDTH, mp_init, integral_fit, SV_SIG)
    ffit.SetParNames("Width", "MP", "Area", "GSigma")
    ffit.SetParLimits(0, PLLO_W, PLHI_W)
    ffit.SetParLimits(1, PLLO_MP, PLHI_MP)
    ffit.SetParLimits(2, integral_fit * 0.01, 1e10)
    ffit.SetParLimits(3, PLLO_S, PLHI_S)
    res = hist.Fit(ffit, "RB0QS")
    mpv = ffit.GetParameter(1)
    ok = PLLO_MP <= mpv <= PLHI_MP
    if res is not None and hasattr(res, "Status"):
        st = res.Status()
        if st not in (0, 1, 3):
            ok = False
    if not ok:
        return None
    return ffit, mpv


def style_hist(h, colour):
    # HIST 阶梯线更适合 ADC 谱叠加对比（比散点 P 更清晰）
    h.SetMarkerSize(0)
    h.SetMarkerStyle(1)
    h.SetLineColor(colour)
    h.SetLineWidth(2)
    h.SetLineStyle(1)
    h.SetFillColor(0)
    h.GetXaxis().SetTitle("ADC")
    h.GetYaxis().SetTitle("Entries / N(ADC>%.0f)" % NORM_ADC_MIN)
    h.GetYaxis().SetTitleOffset(1.25)
    h.GetXaxis().CenterTitle(True)
    h.GetYaxis().CenterTitle(True)
    h.GetXaxis().SetTitleSize(0.045)
    h.GetYaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.040)
    h.GetYaxis().SetLabelSize(0.040)


def count_above_adc(hist, adc_min=NORM_ADC_MIN):
    """统计 ADC > adc_min 的事例数（从第一个中心 > adc_min 的 bin 积到末尾）。"""
    if hist is None:
        return 0
    axis = hist.GetXaxis()
    # FindBin 返回包含 adc_min 的 bin；要严格 > adc_min，从下一个 bin 起算
    b0 = axis.FindBin(adc_min)
    if axis.GetBinCenter(b0) <= adc_min:
        b0 += 1
    b0 = max(1, min(b0, hist.GetNbinsX()))
    return int(hist.Integral(b0, hist.GetNbinsX()))


def load_center_hist(hist_path, hname):
    fin = ROOT.TFile.Open(hist_path, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError("cannot open %s" % hist_path)
    h = fin.Get(hname)
    if h is None:
        fin.Close()
        raise RuntimeError("%s not found in %s" % (hname, hist_path))
    h = h.Clone("%s_%s" % (hname, os.path.basename(os.path.dirname(hist_path))))
    h.SetDirectory(0)
    fin.Close()
    return h


def plot_channel_canvas(tests, histos, channel_tag, title, outpng, crystal=12):
    """单 pad 叠加所有测试谱。"""
    c = ROOT.TCanvas("c_%s" % channel_tag, title, 1024, 768)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.10)
    ROOT.gPad.SetGridx(True)
    ROOT.gPad.SetGridy(True)
    ROOT.gPad.SetTickx(1)
    ROOT.gPad.SetTicky(1)
    ROOT.gPad.SetFrameLineWidth(2)

    keep = []
    drawn = []  # (test, hist, n_norm, mpv, fdraw)
    ymax = 0.0

    for i, (test, hist) in enumerate(zip(tests, histos)):
        colour = COLOURS[i % len(COLOURS)]
        # 归一分母：该测试自身 ADC>5000 的事例数
        n_norm = count_above_adc(hist, NORM_ADC_MIN)
        fit_res = fit_landau_gauss(hist, "%s_%s" % (channel_tag, test))
        mpv = None
        fdraw = None
        if fit_res is not None:
            ffit, mpv = fit_res
            fdraw = ffit.Clone("langau_draw_%s_%s" % (channel_tag, test))
            if n_norm > 0:
                fdraw.SetParameter(2, ffit.GetParameter(2) / float(n_norm))
            fdraw.SetLineColor(colour)
            fdraw.SetLineStyle(2)
            fdraw.SetLineWidth(2)
            keep.append(fdraw)
            print(
                "%s %s MPV=%.1f N(>%.0f)=%d"
                % (channel_tag, test, mpv, NORM_ADC_MIN, n_norm)
            )
        else:
            print(
                "%s %s fit failed N(>%.0f)=%d"
                % (channel_tag, test, NORM_ADC_MIN, n_norm),
                file=sys.stderr,
            )

        if n_norm > 0:
            hist.Scale(1.0 / n_norm)

        style_hist(hist, colour)
        hist.GetXaxis().SetRangeUser(PLOT_XMIN, PLOT_XMAX)
        hist.SetTitle("")
        ymax = max(ymax, hist.GetMaximum())
        drawn.append((test, hist, n_norm, mpv, fdraw))
        keep.append(hist)

    ymax = max(ymax * 1.25, 1e-4)

    for i, (test, hist, n_norm, mpv, fdraw) in enumerate(drawn):
        hist.SetMinimum(0.0)
        hist.SetMaximum(ymax)
        if i == 0:
            hist.Draw("HIST")
        else:
            hist.Draw("HIST SAME")
        if fdraw is not None:
            fdraw.Draw("LSAME")

    # Crystal / 主备高增益 / 1-hit 筛选
    if channel_tag == "hmh":
        ch_txt = "Main HG"
    elif channel_tag == "hbh":
        ch_txt = "Backup HG"
    else:
        ch_txt = channel_tag
    cry_lab = ROOT.TLatex()
    cry_lab.SetNDC(True)
    cry_lab.SetTextFont(62)
    cry_lab.SetTextSize(0.040)
    cry_lab.DrawLatex(0.14, 0.93, "Crystal %d  %s  (1-hit)" % (crystal, ch_txt))
    keep.append(cry_lab)

    # legend：测试名 + N(>5000) + MPV（略向左移）
    n = len(drawn)
    leg_y1 = max(0.45, 0.88 - 0.055 * n)
    leg = ROOT.TLegend(0.38, leg_y1, 0.78, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.030)
    for test, hist, n_norm, mpv, _ in drawn:
        if mpv is not None:
            label = "%s  N(>%.0f)=%d  MPV=%.0f" % (test, NORM_ADC_MIN, n_norm, mpv)
        else:
            label = "%s  N(>%.0f)=%d" % (test, NORM_ADC_MIN, n_norm)
        leg.AddEntry(hist, label, "l")
    leg.Draw()
    keep.append(leg)

    c.Update()
    c.SaveAs(outpng)
    print("saved", outpng)
    c._keep = keep
    # 返回各测试 MPV（拟合失败为 None）
    return {t: mpv for t, _, _, mpv, _ in drawn}


# MPV 折线图横轴顺序（按用户指定）
MPV_TREND_ORDER = ["onlyCalo", "vibration", "firstFSE", "specialItem", "secondFSE"]
# 与 MPV_TREND_ORDER 一一对应的温度（摄氏度）
MPV_TREND_TEMP_C = [34, 28, 26, 26, 34]


def plot_mpv_trend(mpv_hmh, mpv_hbh, outpng, crystal=12):
    """上 pad：MPV；下 pad：相对第一次的波动率 (MPV/MPV0-1)*100%。"""
    labels = [t for t in MPV_TREND_ORDER if t in mpv_hmh or t in mpv_hbh]
    if not labels:
        print("no tests for MPV trend", file=sys.stderr)
        return None

    # 温度与最终 labels 对齐（按 MPV_TREND_ORDER 顺序）
    temps = []
    for t in labels:
        idx = MPV_TREND_ORDER.index(t)
        temps.append(MPV_TREND_TEMP_C[idx])

    c = ROOT.TCanvas("c_mpv_trend", "MPV trend", 1024, 900)
    c.Divide(1, 2, 0.001, 0.01)

    # ---------- 上：MPV ----------
    p1 = c.cd(1)
    p1.SetLeftMargin(0.12)
    p1.SetRightMargin(0.05)
    p1.SetTopMargin(0.12)
    p1.SetBottomMargin(0.02)
    p1.SetGridx(True)
    p1.SetGridy(True)

    n = len(labels)
    frame = ROOT.TH1F("h_mpv_frame", "", n, 0, n)
    frame.SetDirectory(0)
    for i in range(n):
        frame.GetXaxis().SetBinLabel(i + 1, "")
    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTickLength(0.03)
    frame.GetYaxis().SetTitle("MPV (ADC)")
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetTitleOffset(1.05)
    frame.SetStats(0)

    def _make_graph(mpv_map, name, colour, marker):
        xs, ys = [], []
        for i, lab in enumerate(labels):
            v = mpv_map.get(lab)
            if v is None:
                continue
            xs.append(i + 0.5)
            ys.append(float(v))
        if not xs:
            return None
        gr = ROOT.TGraph(len(xs))
        gr.SetName(name)
        for i, (x, y) in enumerate(zip(xs, ys)):
            gr.SetPoint(i, x, y)
        gr.SetLineColor(colour)
        gr.SetMarkerColor(colour)
        gr.SetLineWidth(2)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerSize(1.4)
        return gr

    g_main = _make_graph(mpv_hmh, "g_hmh", ROOT.kCyan + 2, 20)
    g_back = _make_graph(mpv_hbh, "g_hbh", ROOT.kMagenta + 1, 21)
    if g_main is None and g_back is None:
        print("no MPV points", file=sys.stderr)
        return None

    ys_all = []
    for g in (g_main, g_back):
        if g is None:
            continue
        for i in range(g.GetN()):
            ys_all.append(g.GetPointY(i))
    ymin = min(ys_all) - 500.0
    ymax = max(ys_all) + 500.0
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.Draw("AXIS")

    cry_lab = ROOT.TLatex()
    cry_lab.SetNDC(True)
    cry_lab.SetTextFont(62)
    cry_lab.SetTextSize(0.055)
    cry_lab.DrawLatex(0.14, 0.88, "Crystal %d  Main/Backup HG  (1-hit)" % crystal)

    keep = [frame, cry_lab, p1]
    if g_main is not None:
        g_main.Draw("PL SAME")
        keep.append(g_main)
    if g_back is not None:
        g_back.Draw("PL SAME")
        keep.append(g_back)

    pts = []
    for g, mpv_map, dy in (
        (g_main, mpv_hmh, 120.0),
        (g_back, mpv_hbh, -180.0),
    ):
        if g is None:
            continue
        for i, lab in enumerate(labels):
            v = mpv_map.get(lab)
            if v is None:
                continue
            pt = ROOT.TLatex()
            pt.SetTextFont(42)
            pt.SetTextSize(0.040)
            pt.SetTextAlign(22)
            pt.SetTextColor(g.GetMarkerColor())
            pt.DrawLatex(i + 0.5, v + dy, "%.0f" % v)
            pts.append(pt)
    keep.append(pts)

    # Main/Backup legend 仍在 pad1
    leg = ROOT.TLegend(0.65, 0.72, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    if g_main is not None:
        leg.AddEntry(g_main, "Main HG (hmh)", "lp")
    if g_back is not None:
        leg.AddEntry(g_back, "Back HG (hbh)", "lp")
    leg.Draw()
    keep.append(leg)

    # ---------- 下：相对第一次的波动率 ----------
    p2 = c.cd(2)
    p2.SetLeftMargin(0.12)
    p2.SetRightMargin(0.05)
    p2.SetTopMargin(0.02)
    p2.SetBottomMargin(0.34)
    p2.SetGridx(True)
    p2.SetGridy(True)

    def _make_rel_graph(mpv_map, name, colour, marker):
        ref = None
        for lab in labels:
            if mpv_map.get(lab) is not None:
                ref = float(mpv_map[lab])
                break
        if ref is None or ref == 0:
            return None
        xs, ys = [], []
        for i, lab in enumerate(labels):
            v = mpv_map.get(lab)
            if v is None:
                continue
            xs.append(i + 0.5)
            ys.append((float(v) / ref - 1.0) * 100.0)
        if not xs:
            return None
        gr = ROOT.TGraph(len(xs))
        gr.SetName(name)
        for i, (x, y) in enumerate(zip(xs, ys)):
            gr.SetPoint(i, x, y)
        gr.SetLineColor(colour)
        gr.SetMarkerColor(colour)
        gr.SetLineWidth(2)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerSize(1.4)
        return gr, ref

    g_main_r = _make_rel_graph(mpv_hmh, "g_hmh_rel", ROOT.kCyan + 2, 20)
    g_back_r = _make_rel_graph(mpv_hbh, "g_hbh_rel", ROOT.kMagenta + 1, 21)

    frame2 = ROOT.TH1F("h_mpv_rel_frame", "", n, 0, n)
    frame2.SetDirectory(0)
    for i, lab in enumerate(labels):
        frame2.GetXaxis().SetBinLabel(
            i + 1, "#splitline{%s}{~%d^{#circ}C}" % (lab, temps[i])
        )
    frame2.SetTitle("")
    frame2.GetXaxis().SetLabelSize(0.055)
    frame2.GetXaxis().LabelsOption("h")
    frame2.GetXaxis().SetTitle("")  # 温度说明用 TLatex，紧挨横轴标签下方
    frame2.GetYaxis().SetTitle("#DeltaMPV / MPV_{0} [%]")
    frame2.GetYaxis().CenterTitle(True)
    frame2.GetYaxis().SetTitleSize(0.055)
    frame2.GetYaxis().SetLabelSize(0.050)
    frame2.GetYaxis().SetTitleOffset(1.05)
    frame2.SetStats(0)

    ys_rel = []
    for item in (g_main_r, g_back_r):
        if item is None:
            continue
        g, _ = item
        for i in range(g.GetN()):
            ys_rel.append(g.GetPointY(i))
    if ys_rel:
        span = max(abs(min(ys_rel)), abs(max(ys_rel)), 1.0)
        frame2.SetMinimum(-span * 1.45)
        frame2.SetMaximum(span * 1.45)
    else:
        frame2.SetMinimum(-10)
        frame2.SetMaximum(10)
    frame2.Draw("AXIS")

    # 零线
    zero = ROOT.TLine(0, 0, n, 0)
    zero.SetLineStyle(2)
    zero.SetLineColor(ROOT.kGray + 2)
    zero.SetLineWidth(2)
    zero.Draw()

    note = ROOT.TLatex()
    note.SetNDC(True)
    note.SetTextFont(42)
    note.SetTextSize(0.045)
    note.DrawLatex(0.14, 0.88, "relative to first (%s)" % labels[0])

    keep.extend([p2, frame2, zero, note])
    if g_main_r is not None:
        g_main_r[0].Draw("PL SAME")
        keep.append(g_main_r[0])
        print("MPV trend Main ref(%s)=%.1f" % (labels[0], g_main_r[1]))
    if g_back_r is not None:
        g_back_r[0].Draw("PL SAME")
        keep.append(g_back_r[0])
        print("MPV trend Back ref(%s)=%.1f" % (labels[0], g_back_r[1]))

    for item in (g_main_r, g_back_r):
        if item is None:
            continue
        g, _ = item
        for i in range(g.GetN()):
            x = g.GetPointX(i)
            y = g.GetPointY(i)
            pt = ROOT.TLatex()
            pt.SetTextFont(42)
            pt.SetTextSize(0.040)
            pt.SetTextAlign(22)
            pt.SetTextColor(g.GetMarkerColor())
            pt.DrawLatex(x, y + (1.5 if y >= 0 else -2.2), "%.1f%%" % y)
            keep.append(pt)

    # 紧挨 pad2 横轴分类标签正下方
    p2.cd()
    temp_note = ROOT.TLatex()
    temp_note.SetNDC(True)
    temp_note.SetTextFont(42)
    temp_note.SetTextSize(0.045)
    temp_note.SetTextAlign(22)
    temp_note.DrawLatex(
        0.50, 0.08, "average temperature, deviation within 2^{#circ}C"
    )
    keep.append(temp_note)

    c.Update()
    c.SaveAs(outpng)
    print("saved", outpng)
    for lab, tc in zip(labels, temps):
        print("  %s -> %d C" % (lab, tc))
    c._keep = keep
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tests", nargs="+", required=True)
    ap.add_argument("--output-dir", required=True, help="各测试子目录根 (含 <test>/hist_calo.root)")
    ap.add_argument("--compare-dir", required=True)
    ap.add_argument("--crystal", type=int, default=12)
    args = ap.parse_args()

    os.makedirs(args.compare_dir, exist_ok=True)
    cry = args.crystal
    hmh_name = "hmh_%d" % cry
    hbh_name = "hbh_%d" % cry

    tests, hmh_list, hbh_list = [], [], []
    for name in args.tests:
        hist_path = os.path.join(args.output_dir, name, "hist_calo.root")
        if not os.path.isfile(hist_path):
            print("[skip] no hist:", hist_path, file=sys.stderr)
            continue
        try:
            hmh = load_center_hist(hist_path, hmh_name)
            hbh = load_center_hist(hist_path, hbh_name)
        except Exception as e:
            print("[skip] %s: %s" % (name, e), file=sys.stderr)
            continue
        tests.append(name)
        hmh_list.append(hmh)
        hbh_list.append(hbh)

    if not tests:
        print("no valid hist_calo.root", file=sys.stderr)
        sys.exit(1)

    out_main = os.path.join(
        args.compare_dir, "center_cry%d_hmh_mip_compare.png" % cry
    )
    out_back = os.path.join(
        args.compare_dir, "center_cry%d_hbh_mip_compare.png" % cry
    )
    mpv_hmh = plot_channel_canvas(
        tests,
        hmh_list,
        "hmh",
        "Crystal %d Main HG MIP (1-hit)" % cry,
        out_main,
        crystal=cry,
    )
    mpv_hbh = plot_channel_canvas(
        tests,
        hbh_list,
        "hbh",
        "Crystal %d Back HG MIP (1-hit)" % cry,
        out_back,
        crystal=cry,
    )

    out_trend = os.path.join(
        args.compare_dir, "center_cry%d_mpv_trend.png" % cry
    )
    plot_mpv_trend(mpv_hmh, mpv_hbh, out_trend, crystal=cry)


if __name__ == "__main__":
    main()
