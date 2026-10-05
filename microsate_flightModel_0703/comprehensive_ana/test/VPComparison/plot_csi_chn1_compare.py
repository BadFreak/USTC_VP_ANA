#!/usr/bin/env python3
"""CSI chn1 (CellID==100000) MIP 对比：各测试归一化谱叠加 + MPV 折线。

谱：CellADC - CellPLAT；拟合参数与 cosmic_ped/scripts/draw_csi.cxx 一致。
无 onlyCalo。只在 VPComparison 内写输出。
"""
from __future__ import annotations

import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# 与 draw_csi.cxx FitLandauGauss 一致
FIT_XLO, FIT_XHI = 1000.0, 14000.0
SV_WIDTH, SV_SIG = 300.0, 150.0
PLLO_MP, PLHI_MP = 1500.0, 12000.0
PLLO_W, PLHI_W = 50.0, 1500.0
PLLO_S, PLHI_S = 20.0, 500.0

# 与 draw_csi CellADC_4x2 画图范围一致
PLOT_XMIN, PLOT_XMAX = 1000.0, 14000.0
NORM_ADC_MIN = 1000.0

CSI_CHN1_CELLID = 100000
HIST_NAME = "h_chn1"
NBINS = 200
HIST_XMIN, HIST_XMAX = 1000.0, 14000.0

COLOURS = [
    ROOT.kCyan + 2,
    ROOT.kBlue + 1,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kGreen + 2,
    ROOT.kBlack,
]

# MPV 折线横轴顺序（无 onlyCalo）
MPV_TREND_ORDER = ["vibration", "firstFSE", "specialItem", "secondFSE"]
# 与 MPV_TREND_ORDER 一一对应的温度（摄氏度）
MPV_TREND_TEMP_C = [24, 24, 24, 31]


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


def count_above_adc(hist, adc_min=NORM_ADC_MIN):
    if hist is None:
        return 0
    axis = hist.GetXaxis()
    b0 = axis.FindBin(adc_min)
    if axis.GetBinCenter(b0) <= adc_min:
        b0 += 1
    b0 = max(1, min(b0, hist.GetNbinsX()))
    return int(hist.Integral(b0, hist.GetNbinsX()))


def fill_chn1_hist(resultfile, outfile, cell_id=CSI_CHN1_CELLID):
    """从 csiTree 填 CellADC-CellPLAT (CellID==chn1)，写入 hist_csi.root。"""
    fin = ROOT.TFile.Open(resultfile, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError("cannot open %s" % resultfile)
    tree = fin.Get("csiTree")
    if not tree:
        fin.Close()
        raise RuntimeError("no csiTree in %s" % resultfile)

    cut = "CellID==%d" % cell_id
    n = tree.Draw(
        "CellADC-CellPLAT>>htmp_csi(%d,%g,%g)" % (NBINS, HIST_XMIN, HIST_XMAX),
        cut,
        "goff",
    )
    htmp = tree.GetHistogram()
    if htmp is None:
        fin.Close()
        raise RuntimeError("Draw failed for %s" % resultfile)
    h = htmp.Clone(HIST_NAME)
    h.SetTitle("CSI chn1 CellADC-CellPLAT")
    h.SetDirectory(0)
    fin.Close()

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    fout = ROOT.TFile.Open(outfile, "RECREATE")
    h.Write()
    fout.Close()
    print("saved %s  entries=%d  drawn=%d" % (outfile, int(h.GetEntries()), int(n)))
    return h


def load_hist(hist_path, hname=HIST_NAME):
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


def style_hist(h, colour):
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


def plot_overlay(tests, histos, outpng, channel_label="chn1"):
    c = ROOT.TCanvas("c_csi_%s" % channel_label, "", 1024, 768)
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
    drawn = []
    ymax = 0.0

    for i, (test, hist) in enumerate(zip(tests, histos)):
        colour = COLOURS[i % len(COLOURS)]
        n_norm = count_above_adc(hist, NORM_ADC_MIN)
        fit_res = fit_landau_gauss(hist, "csi_%s_%s" % (channel_label, test))
        mpv = None
        fdraw = None
        if fit_res is not None:
            ffit, mpv = fit_res
            fdraw = ffit.Clone("langau_csi_%s_%s" % (channel_label, test))
            if n_norm > 0:
                fdraw.SetParameter(2, ffit.GetParameter(2) / float(n_norm))
            fdraw.SetLineColor(colour)
            fdraw.SetLineStyle(2)
            fdraw.SetLineWidth(2)
            keep.append(fdraw)
            print("csi %s %s MPV=%.1f N(>%.0f)=%d" % (channel_label, test, mpv, NORM_ADC_MIN, n_norm))
        else:
            print(
                "csi %s %s fit failed N(>%.0f)=%d" % (channel_label, test, NORM_ADC_MIN, n_norm),
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

    cry_lab = ROOT.TLatex()
    cry_lab.SetNDC(True)
    cry_lab.SetTextFont(62)
    cry_lab.SetTextSize(0.045)
    cry_lab.DrawLatex(0.14, 0.93, "CSI %s (CellID=%d)" % (channel_label, CSI_CHN1_CELLID))
    keep.append(cry_lab)

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
    return {t: mpv for t, _, _, mpv, _ in drawn}


def plot_mpv_trend(mpv_map, outpng, channel_label="chn1"):
    """上 pad：MPV；下 pad：相对第一次的波动率；横轴带温度。"""
    labels = [t for t in MPV_TREND_ORDER if t in mpv_map and mpv_map[t] is not None]
    if not labels:
        print("no tests for CSI MPV trend", file=sys.stderr)
        return None

    temps = []
    for t in labels:
        idx = MPV_TREND_ORDER.index(t)
        temps.append(MPV_TREND_TEMP_C[idx])

    c = ROOT.TCanvas("c_csi_mpv_trend", "CSI MPV trend", 1024, 900)
    c.Divide(1, 2, 0.001, 0.001)

    # ---------- 上：MPV ----------
    p1 = c.cd(1)
    p1.SetLeftMargin(0.12)
    p1.SetRightMargin(0.05)
    p1.SetTopMargin(0.12)
    p1.SetBottomMargin(0.02)
    p1.SetGridx(True)
    p1.SetGridy(True)

    n = len(labels)
    frame = ROOT.TH1F("h_csi_mpv_frame", "", n, 0, n)
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

    gr = ROOT.TGraph(len(labels))
    gr.SetName("g_csi_chn1")
    ys = []
    for i, lab in enumerate(labels):
        v = float(mpv_map[lab])
        gr.SetPoint(i, i + 0.5, v)
        ys.append(v)
    gr.SetLineColor(ROOT.kCyan + 2)
    gr.SetMarkerColor(ROOT.kCyan + 2)
    gr.SetLineWidth(2)
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.4)

    ymin = min(ys) - 400.0
    ymax = max(ys) + 400.0
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.Draw("AXIS")

    cry_lab = ROOT.TLatex()
    cry_lab.SetNDC(True)
    cry_lab.SetTextFont(62)
    cry_lab.SetTextSize(0.055)
    cry_lab.DrawLatex(0.14, 0.88, "CSI %s (CellID=%d)" % (channel_label, CSI_CHN1_CELLID))

    gr.Draw("PL SAME")

    pts = []
    for i, lab in enumerate(labels):
        v = float(mpv_map[lab])
        pt = ROOT.TLatex()
        pt.SetTextFont(42)
        pt.SetTextSize(0.040)
        pt.SetTextAlign(22)
        pt.SetTextColor(ROOT.kCyan + 2)
        pt.DrawLatex(i + 0.5, v + 80.0, "%.0f" % v)
        pts.append(pt)

    leg = ROOT.TLegend(0.65, 0.72, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry(gr, "CSI %s" % channel_label, "lp")
    leg.Draw()

    keep = [frame, cry_lab, gr, pts, leg, p1]

    # ---------- 下：相对第一次的波动率 ----------
    p2 = c.cd(2)
    p2.SetLeftMargin(0.12)
    p2.SetRightMargin(0.05)
    p2.SetTopMargin(0.02)
    p2.SetBottomMargin(0.32)
    p2.SetGridx(True)
    p2.SetGridy(True)

    ref = float(mpv_map[labels[0]])
    gr_r = ROOT.TGraph(len(labels))
    gr_r.SetName("g_csi_chn1_rel")
    ys_rel = []
    for i, lab in enumerate(labels):
        y = (float(mpv_map[lab]) / ref - 1.0) * 100.0
        gr_r.SetPoint(i, i + 0.5, y)
        ys_rel.append(y)
    gr_r.SetLineColor(ROOT.kCyan + 2)
    gr_r.SetMarkerColor(ROOT.kCyan + 2)
    gr_r.SetLineWidth(2)
    gr_r.SetMarkerStyle(20)
    gr_r.SetMarkerSize(1.4)

    frame2 = ROOT.TH1F("h_csi_mpv_rel_frame", "", n, 0, n)
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

    span = max(abs(min(ys_rel)), abs(max(ys_rel)), 1.0)
    frame2.SetMinimum(-span * 1.45)
    frame2.SetMaximum(span * 1.45)
    frame2.Draw("AXIS")

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

    gr_r.Draw("PL SAME")

    rpts = []
    for i in range(gr_r.GetN()):
        x = gr_r.GetPointX(i)
        y = gr_r.GetPointY(i)
        pt = ROOT.TLatex()
        pt.SetTextFont(42)
        pt.SetTextSize(0.040)
        pt.SetTextAlign(22)
        pt.SetTextColor(ROOT.kCyan + 2)
        pt.DrawLatex(x, y + (1.5 if y >= 0 else -2.2), "%.1f%%" % y)
        rpts.append(pt)

    # 紧挨 pad2 横轴分类标签正下方（不用 SetTitle，避免 offset 过大/裁切）
    p2.cd()
    temp_note = ROOT.TLatex()
    temp_note.SetNDC(True)
    temp_note.SetTextFont(42)
    temp_note.SetTextSize(0.045)
    temp_note.SetTextAlign(22)
    temp_note.DrawLatex(
        0.50, 0.12, "average temperature, deviation within 2^{#circ}C"
    )

    keep.extend([p2, frame2, zero, note, gr_r, rpts, temp_note])
    c.Update()
    c.SaveAs(outpng)
    print("saved", outpng)
    print("MPV trend CSI ref(%s)=%.1f" % (labels[0], ref))
    for lab, tc in zip(labels, temps):
        print("  %s -> %d C" % (lab, tc))
    c._keep = keep
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tests", nargs="+", required=True)
    ap.add_argument("--input-dir", required=True, help="input/csi 目录")
    ap.add_argument("--output-dir", required=True, help="output/csi 根目录")
    ap.add_argument("--compare-dir", required=True)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.compare_dir, exist_ok=True)

    tests = []
    histos = []
    for name in args.tests:
        resultfile = os.path.join(args.input_dir, "result_%s.root" % name)
        work = os.path.join(args.output_dir, name)
        hist_path = os.path.join(work, "hist_csi.root")
        os.makedirs(work, exist_ok=True)

        if not args.plot_only and (args.force or not os.path.isfile(hist_path)):
            if not os.path.isfile(resultfile):
                print("[skip] missing", resultfile, file=sys.stderr)
                continue
            try:
                fill_chn1_hist(resultfile, hist_path)
            except Exception as e:
                print("[skip] fill %s: %s" % (name, e), file=sys.stderr)
                continue
        if not os.path.isfile(hist_path):
            print("[skip] no hist:", hist_path, file=sys.stderr)
            continue
        try:
            h = load_hist(hist_path)
        except Exception as e:
            print("[skip] load %s: %s" % (name, e), file=sys.stderr)
            continue
        tests.append(name)
        histos.append(h)

    if not tests:
        print("no valid CSI hist", file=sys.stderr)
        sys.exit(1)

    out_overlay = os.path.join(args.compare_dir, "csi_chn1_mip_compare.png")
    out_trend = os.path.join(args.compare_dir, "csi_chn1_mpv_trend.png")
    mpv_map = plot_overlay(tests, histos, out_overlay, channel_label="chn1")
    plot_mpv_trend(mpv_map, out_trend, channel_label="chn1")


if __name__ == "__main__":
    main()
