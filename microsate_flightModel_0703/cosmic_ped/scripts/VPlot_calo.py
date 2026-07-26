# An example python script for making a tidy plot with several lines on it
# I'm using TProfiles because that's what is in the TIDA output, but
# the standard histogram is a "TH1F" and the exact same commands work for it.
# If you want to use graphs instead of histograms you can check out
# the other example script in here.

import ROOT
import math
import os
import sys

# Imports ATLAS style for plotting
# You have to have set it up first (see README for instructions)
# You can run it without this but it will have an ugly stats box and so on
# that you'd have to turn off manually.
# import VPStyle
# ROOT.SetVPStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")

# hmh/hbh/hml/hbl 多 pad 图是否做 Landau×Gauss 拟合并画红线与 MPV；也可用 plot_sig(..., do_fit=False) 或命令行 --no-sig-fit
SIG_FIT_ENABLED = True

PLOT_ALL_ROOT = "plot_all.root"

# main_hlr / back_hlr 的 TGraph 线性拟合区间（与 draw_one_ratio 一致）
HL_RATIO_FIT_XLO = 200.
HL_RATIO_FIT_XHI = 16000.


def _write_canvas_to_plot_all(canvas, key_name):
    """把 canvas 写入 plot_all.root（同名则覆盖）。"""
    mode = "UPDATE" if os.path.isfile(PLOT_ALL_ROOT) else "RECREATE"
    fout = ROOT.TFile.Open(PLOT_ALL_ROOT, mode)
    if not fout or fout.IsZombie():
        print("[VPlot] cannot open %s for write" % PLOT_ALL_ROOT)
        return
    fout.cd()
    canvas.Write(key_name, ROOT.TObject.kOverwrite)
    fout.Close()
    print("[VPlot] wrote %s -> %s" % (key_name, PLOT_ALL_ROOT))


def _hl_ratio_graph_name(ratio_name, channel):
    """main_hlr -> main_hlg_i，back_hlr -> back_hlg_i"""
    if ratio_name.endswith("_hlr"):
        return ratio_name[:-4] + "_hlg_" + str(channel)
    return ratio_name + "_" + str(channel)


def _fit_hl_ratio_graph(gr, tag):
    if gr is None or gr.GetN() <= 0:
        return None, None
    ffit = ROOT.TF1("pol1_" + tag, "pol1", HL_RATIO_FIT_XLO, HL_RATIO_FIT_XHI)
    ffit.SetLineColor(ROOT.kRed)
    ffit.SetLineWidth(2)
    gr.Fit(ffit, "RB0QS")
    return ffit, ffit.GetParameter(1)


def _draw_hl_ratio_slope_label(channel, slope):
    pt = ROOT.TPaveText(0.14, 0.72, 0.56, 0.90, "NDC")
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.048)
    pt.SetTextColor(ROOT.kBlack)
    pt.AddText("%d  slope: %.3f" % (channel, slope))
    pt.Draw()
    return pt


def _style_hl_ratio_graph(gr, colour):
    gr.SetMarkerColor(colour)
    gr.SetMarkerStyle(3)
    gr.SetMarkerSize(0.55)
    gr.SetLineColor(colour)
    gr.SetLineWidth(1)
    gr.SetFillColor(0)


def _hk_yrange_from_tgraph(gr, is_temp):
    """由 TGraph 估计 Y 轴范围；优先用物理合理点，否则用全部有限点（避免全 -1 时 y 轴为 0~1 看不见）。"""
    if not gr:
        return 0.0, 1.0
    np = gr.GetN()
    if np <= 0:
        return 0.0, 1.0
    ys = []
    ys_all = []
    for ip in range(np):
        y = gr.GetPointY(ip)
        if not math.isfinite(y):
            continue
        ys_all.append(y)
        if y > -0.5:
            ys.append(y)
    use = ys if ys else ys_all
    if not use:
        return 0.0, 1.0
    mn = min(use)
    mx = max(use)
    span = mx - mn
    if span <= 0:
        span = max(abs(mx), 1.0) * 0.1
    if is_temp:
        return mn - span, mx + span
    return mn - 2.0 * span, mx + 2.0 * span


def _xrange_from_tgraph(gr):
    if not gr or gr.GetN() <= 0:
        return 0.0, 1.0
    xs = [gr.GetPointX(i) for i in range(gr.GetN()) if math.isfinite(gr.GetPointX(i))]
    if not xs:
        return 0.0, 1.0
    return min(xs), max(xs)


def _style_hk_tgraph(gr, colour, lw=2):
    gr.SetLineColor(colour)
    gr.SetLineWidth(lw)
    gr.SetMarkerSize(0)
    gr.GetXaxis().SetTitle("hkTree row index (decode order)")
    gr.GetYaxis().SetTitleOffset(1.2)


def plot_ct():
    """读取 temp_hist.root（draw_temp 写的 TGraph f{TelFEE}c* / f{TelFEE}t*），遥测 FEE 为 2 与 5–8。"""
    if not os.path.isfile("temp_hist.root"):
        return
    fin = ROOT.TFile.Open("temp_hist.root", "READ")
    if not fin or fin.IsZombie():
        return

    cols_t = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]
    cols_c = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2]
    fee_ids = (2, 5, 6, 7, 8)

    for fee in fee_ids:
        cname_t = "fee%d_temp" % fee
        ct = ROOT.TCanvas(cname_t, "FEE %d temperatures" % fee, 900, 900)
        ct.Divide(2, 2)
        for it in range(4):
            ct.cd(it + 1)
            gr = fin.Get("f%dt%d" % (fee, it))
            if not gr:
                continue
            _style_hk_tgraph(gr, cols_t[it])
            gr.GetYaxis().SetTitle("T [#circ C]")
            np = gr.GetN()
            if np <= 0:
                lat = ROOT.TLatex(0.35, 0.5, "no points (TelFEE=%d)" % fee)
                lat.SetNDC(True)
                lat.Draw()
            else:
                lo, hi = _hk_yrange_from_tgraph(gr, True)
                gr.GetYaxis().SetRangeUser(lo, hi)
                xlo, xhi = _xrange_from_tgraph(gr)
                if xhi > xlo:
                    gr.GetXaxis().SetRangeUser(xlo, xhi)
                # batch 下 TGraph 必须用 AL 才能画出坐标轴与折线，仅 L 常为空白
                gr.Draw("AL")
            ROOT.gPad.Modified()
            ROOT.gPad.Update()
        ct.SaveAs("fee%d_temp.png" % fee)

        cname_i = "fee%d_cur" % fee
        cc = ROOT.TCanvas(cname_i, "FEE %d currents" % fee, 1200, 420)
        cc.Divide(3, 1)
        for ic in range(3):
            cc.cd(ic + 1)
            gr = fin.Get("f%dc%d" % (fee, ic))
            if not gr:
                continue
            _style_hk_tgraph(gr, cols_c[ic])
            gr.GetYaxis().SetTitle("I [mA]")
            np = gr.GetN()
            if np <= 0:
                lat = ROOT.TLatex(0.35, 0.5, "no points (TelFEE=%d)" % fee)
                lat.SetNDC(True)
                lat.Draw()
            else:
                lo, hi = _hk_yrange_from_tgraph(gr, False)
                gr.GetYaxis().SetRangeUser(lo, hi)
                xlo, xhi = _xrange_from_tgraph(gr)
                if xhi > xlo:
                    gr.GetXaxis().SetRangeUser(xlo, xhi)
                gr.Draw("AL")
            ROOT.gPad.Modified()
            ROOT.gPad.Update()
        cc.SaveAs("fee%d_current.png" % fee)

    fin.Close()


def plot_sig(name, do_fit=None):
    if do_fit is None:
        do_fit = SIG_FIT_ENABLED
    infile = ROOT.TFile.Open("hist_calo.root", "READ")
    histos = []
    fit_info_list = []  # 每个 pad 的 (ffit, integral_above, mpv)
    if name in ["hmh", "hbh"]:
        # 与 draw_calo 临时增益 hist 临时范围一致（5000–30000）
        fit_xlo, fit_xhi = 5000., 30000.
        sv_mp, sv_width, sv_sig = 9000., 800., 400.
        pllo_mp, plhi_mp = 5000., 20000.
        pllo_w, plhi_w = 100., 3000.
        pllo_s, plhi_s = 50., 2000.
    else:
        fit_xlo, fit_xhi = 600., 3000.
        sv_mp, sv_width, sv_sig = 800., 150., 80.
        pllo_mp, plhi_mp = 600., 2000.
        pllo_w, plhi_w = 30., 500.
        pllo_s, plhi_s = 10., 300.

    for i in range(0, 25):
        histo = infile.Get(name + "_" + str(i))
        if histo is None or not hasattr(histo, "SetDirectory"):
            continue
        histo.SetDirectory(0)
        integral_full = histo.Integral()
        fit_info_i = None
        if do_fit and integral_full > 0:
            hfit = histo.Clone()  # 拟合用归一化后的直方图在拟合区间的形状
            nbins = hfit.GetNbinsX()
            blo = max(1, min(hfit.GetXaxis().FindBin(fit_xlo), nbins))
            bhi = max(1, min(hfit.GetXaxis().FindBin(fit_xhi), nbins))
            if blo > bhi:
                blo, bhi = bhi, blo
            integral_fit = hfit.Integral(blo, bhi)
            if integral_fit > 0:
                imax = blo
                for b in range(blo, bhi + 1):
                    if hfit.GetBinContent(b) > hfit.GetBinContent(imax):
                        imax = b
                mp_init = hfit.GetBinCenter(imax)
                mp_init = max(pllo_mp, min(plhi_mp, mp_init))
                fname = "langau_%s_%d" % (name, i)
                ffit = ROOT.TF1(fname, _langaufun, fit_xlo, fit_xhi, 4)
                ffit.SetParameters(sv_width, mp_init, integral_fit, sv_sig)
                ffit.SetParNames("Width", "MP", "Area", "GSigma")
                ffit.SetParLimits(0, pllo_w, plhi_w)
                ffit.SetParLimits(1, pllo_mp, plhi_mp)
                ffit.SetParLimits(2, integral_fit * 0.01, 1e10)
                ffit.SetParLimits(3, pllo_s, plhi_s)
                res = hfit.Fit(ffit, "RB0QS")
                mpv = ffit.GetParameter(1)
                ok = pllo_mp <= mpv <= plhi_mp
                if res is not None and hasattr(res, "Status"):
                    st = res.Status()
                    if st not in (0, 1, 3):
                        ok = False
                if ok:
                    fit_info_i = (ffit, integral_full, mpv)
        histos.append(histo)
        fit_info_list.append(fit_info_i)

    infile.Close()
    n_ok = sum(1 for x in fit_info_list if x is not None)
    if do_fit and name in ["hmh", "hbh", "hml", "hbl"] and n_ok == 0 and len(histos) > 0:
        print("[VPlot] plot_sig(%s): no pad had successful fit (0/%d)" % (name, len(histos)))

    c = ROOT.TCanvas("canvas", '', 0, 0, 1500, 1500)
    c.Divide(5, 5)
    c.SetGridx(0)
    c.SetGridy(0)
    goodColours = [ROOT.kCyan + 2, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange, ROOT.kBlack]
    pad_fit_refs = []  # (index, ffit, integral_above, mpv) 第二遍画拟合线
    pad_mpv_refs = []  # (index, mpv_val) 第三遍画 MPV 文字

    for index, histo in enumerate(histos):
        colour = goodColours[0]
        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(0)
        histo.SetMarkerStyle(0)
        histo.SetLineColor(colour)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)
        histo.SetFillColor(0)
        histo.GetXaxis().SetTitle("ADC")
        histo.GetYaxis().SetTitle("Entries")
        histo.GetYaxis().SetTitleOffset(1.5)
        # 横轴与 Landau 拟合区间 fit_xlo/fit_xhi 一致（见函数开头）
        x1 = max(fit_xlo, histo.GetXaxis().GetXmin())
        x2 = min(fit_xhi, histo.GetXaxis().GetXmax())
        histo.GetXaxis().SetRangeUser(x1, x2)
        c.cd(index + 1)
        ROOT.gPad.SetLogx(False)
        ROOT.gPad.SetLogy(False)

        if index == 0:
            histo.Draw("H")
        else:
            histo.Draw("H")
        if do_fit and name in ["hmh", "hbh", "hml", "hbl"] and index < len(fit_info_list) and fit_info_list[
            index] is not None:
            pad_fit_refs.append((index,) + fit_info_list[index])

        legend = ROOT.TLegend(0.6, 0.72, 0.92, 0.92)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.SetBorderSize(0)
        legend.SetLineColor(0)
        legend.SetLineStyle(1)
        legend.SetLineWidth(3)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.AddEntry(histo, str(index))
        legend.Draw()

    # 第二遍：逐个 pad 画拟合曲线（不归一化，与计数直方图一致）
    fit_curve_keep = []
    for item in pad_fit_refs:
        index, ffit, _, mpv_val = item[0], item[1], item[2], item[3]
        c.cd(index + 1)
        ffit_draw = ffit.Clone("langau_norm_%s_%d" % (name, index))
        ffit_draw.SetLineColor(ROOT.kRed)
        ffit_draw.SetLineWidth(2)
        ffit_draw.Draw("LSAME")
        fit_curve_keep.append(ffit_draw)
        pad_mpv_refs.append((index, mpv_val))

    # 第三遍：逐个 pad 画 MPV（左上角）
    pt_keep = []
    for index, mpv_val in pad_mpv_refs:
        c.cd(index + 1)
        pt = ROOT.TPaveText(0.58, 0.72, 0.88, 0.90, "NDC")
        pt.SetName("mpv_%s_%d" % (name, index))
        pt.SetBorderSize(0)
        pt.SetFillStyle(0)
        pt.SetTextFont(42)
        pt.SetTextSize(0.07)
        pt.SetTextColor(ROOT.kBlack)
        pt.AddText("MPV = %.0f" % mpv_val)
        # pt.AddText("Entries = %d" % histos[index].GetEntries())
        # 对hmh和hbh，统计大于2000ADC的entries
        if name in ["hmh", "hbh"]:
            hist = histos[index]
            xaxis = hist.GetXaxis()
            bin_2000 = xaxis.FindBin(2000)
            entries_above_2000 = int(hist.Integral(bin_2000, hist.GetNbinsX()))
            pt.AddText("Count>2000 = %d" % entries_above_2000)
        pt.Draw()
        pt_keep.append(pt)

    key = "Sig_" + name
    c.SaveAs(key + ".png")
    _write_canvas_to_plot_all(c, key)


# 廊道卷积高斯 (Landau*Gaussian) 拟合，基于 ROOT langaus.C
def _langaufun(x, par):
    # par[0]=Width(Landau), par[1]=MP, par[2]=Area, par[3]=GSigma
    if par[3] <= 0.:
        return 0.
    invsq2pi = 0.3989422804014
    mpshift = -0.22278298
    np = 100.0
    sc = 5.0
    mpc = par[1] - mpshift * par[0]
    xlow = x[0] - sc * par[3]
    xupp = x[0] + sc * par[3]
    step = (xupp - xlow) / np
    s = 0.0
    for i in range(1, int(np / 2) + 1):
        xx = xlow + (i - 0.5) * step
        fland = ROOT.TMath.Landau(xx, mpc, par[0]) / par[0]
        s += fland * ROOT.TMath.Gaus(x[0], xx, par[3])
        xx = xupp - (i - 0.5) * step
        fland = ROOT.TMath.Landau(xx, mpc, par[0]) / par[0]
        s += fland * ROOT.TMath.Gaus(x[0], xx, par[3])
    return par[2] * step * s * invsq2pi / par[3]


def get_mpv_list(name, hist_path="hist_calo.root"):
    """对 name (hmh/hbh/hml/hbl) 从 hist_path 读直方图，做 Landau*Gauss 拟合，返回 [(crystal_index, MPV), ...]，拟合失败为 (i, None)。"""
    infile = ROOT.TFile.Open(hist_path, "READ")
    if not infile or infile.IsZombie():
        return []
    if name in ["hmh", "hbh"]:
        fit_xlo, fit_xhi = 5000., 30000.
        sv_mp, sv_width, sv_sig = 9000., 800., 400.
        pllo_mp, plhi_mp = 5000., 20000.
        pllo_w, plhi_w = 100., 3000.
        pllo_s, plhi_s = 50., 2000.
    else:
        fit_xlo, fit_xhi = 600., 3000.
        sv_mp, sv_width, sv_sig = 800., 150., 80.
        pllo_mp, plhi_mp = 600., 2000.
        pllo_w, plhi_w = 30., 500.
        pllo_s, plhi_s = 10., 300.
    out = []
    for i in range(25):
        histo = infile.Get(name + "_" + str(i))
        if histo is None or not hasattr(histo, "SetDirectory"):
            out.append((i, None))
            continue
        histo.SetDirectory(0)
        hfit = histo.Clone()
        nbins = hfit.GetNbinsX()
        blo = max(1, min(hfit.GetXaxis().FindBin(fit_xlo), nbins))
        bhi = max(1, min(hfit.GetXaxis().FindBin(fit_xhi), nbins))
        if blo > bhi:
            blo, bhi = bhi, blo
        integral_fit = hfit.Integral(blo, bhi)
        if integral_fit <= 0:
            out.append((i, None))
            continue
        imax = blo
        for b in range(blo, bhi + 1):
            if hfit.GetBinContent(b) > hfit.GetBinContent(imax):
                imax = b
        mp_init = hfit.GetBinCenter(imax)
        mp_init = max(pllo_mp, min(plhi_mp, mp_init))
        fname = "langau_mpv_%s_%d" % (name, i)
        ffit = ROOT.TF1(fname, _langaufun, fit_xlo, fit_xhi, 4)
        ffit.SetParameters(sv_width, mp_init, integral_fit, sv_sig)
        ffit.SetParNames("Width", "MP", "Area", "GSigma")
        ffit.SetParLimits(0, pllo_w, plhi_w)
        ffit.SetParLimits(1, pllo_mp, plhi_mp)
        ffit.SetParLimits(2, integral_fit * 0.01, 1e10)
        ffit.SetParLimits(3, pllo_s, plhi_s)
        res = hfit.Fit(ffit, "RB0QS")
        mpv = ffit.GetParameter(1)
        ok = pllo_mp <= mpv <= plhi_mp
        if res is not None and hasattr(res, "Status"):
            if res.Status() not in (0, 1, 3):
                ok = False
        out.append((i, round(mpv, 2) if ok else None))
    infile.Close()
    return out


def plot_sig_hg_langaus():
    """对 Main_HG(hmh) 和 Back_HG(hbh) 在 5500-30000 ADC 做 Landau*Gaussian 拟合并出图"""
    fit_xlo, fit_xhi = 5500., 30000.
    infile = ROOT.TFile.Open("hist.root", "READ")
    for hname, title, outname in [("hmh", "Main HG", "Sig_hmh_langaus"), ("hbh", "Back HG", "Sig_hbh_langaus")]:
        histo = infile.Get(hname)
        if histo is None or not hasattr(histo, "Fit"):
            continue
        histo = histo.Clone(hname + "_fit")
        histo.SetDirectory(0)
        nbins = histo.GetNbinsX()
        bin_lo = max(1, min(histo.GetXaxis().FindBin(fit_xlo), nbins))
        bin_hi = max(1, min(histo.GetXaxis().FindBin(fit_xhi), nbins))
        if bin_lo > bin_hi:
            bin_lo, bin_hi = bin_hi, bin_lo
        integral = histo.Integral(bin_lo, bin_hi)
        if integral <= 0.:
            continue
        # 初值: Landau Width, MP, Area, GSigma (ADC 尺度)
        sv = [800., 9000., integral * 0.1, 400.]
        pllo = [100., 5500., 1., 50.]
        plhi = [3000., 20000., 1e10, 2000.]
        fname = "langaufcn_" + hname
        ffit = ROOT.TF1(fname, _langaufun, fit_xlo, fit_xhi, 4)
        ffit.SetParameters(sv[0], sv[1], sv[2], sv[3])
        ffit.SetParNames("Width", "MP", "Area", "GSigma")
        for i in range(4):
            ffit.SetParLimits(i, pllo[i], plhi[i])
        histo.Fit(ffit, "RB0Q")
        c = ROOT.TCanvas("c_" + hname, "", 800, 600)
        histo.GetXaxis().SetRangeUser(fit_xlo, fit_xhi)
        histo.SetTitle(title + " (5500-30000 ADC);ADC;Entries")
        histo.Draw("H")
        ffit.SetTitle("Landau#otimesGauss")
        ffit.SetLineColor(ROOT.kRed)
        ffit.Draw("LSAME")
        c.BuildLegend()
        c.SaveAs(outname + ".png")
    infile.Close()


def plot_ratio(name):
    infile = ROOT.TFile.Open("hist_calo.root", "READ")
    is_hl = "hlr" in name
    histos = []
    graphs = []

    for i in range(0, 25):
        histo = infile.Get(name + "_" + str(i))
        if histo is not None:
            histo.SetDirectory(0)
            if histo.Integral() > 0.:
                histo.Scale(1. / histo.Integral())
        histos.append(histo)

        gr = None
        if is_hl:
            gr = infile.Get(_hl_ratio_graph_name(name, i))
            if gr is not None:
                gr = gr.Clone()
        graphs.append(gr)

    infile.Close()

    c = ROOT.TCanvas("canvas", '', 0, 0, 1500, 1500)
    c.Divide(5, 5)
    c.SetGridx(0)
    c.SetGridy(0)

    goodColours = [ROOT.kCyan + 2, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange, ROOT.kBlack]
    colour = goodColours[0]
    keep = []  # 防止 PyROOT 回收对象导致前面 pad 被清空

    for index in range(25):
        c.cd(index + 1)
        ROOT.gPad.SetLogx(False)
        ROOT.gPad.SetLogy(False)
        ROOT.gPad.SetGridx(True)
        ROOT.gPad.SetGridy(True)

        if is_hl and graphs[index] is not None and graphs[index].GetN() > 0:
            gr = graphs[index]
            _style_hl_ratio_graph(gr, colour)
            frame = ROOT.TH2F(
                "hframe_hlg_%s_%d" % (name, index), "", 100, 0., 2e4, 100, 0., 2e3)
            frame.GetXaxis().SetTitle("HighADC")
            frame.GetYaxis().SetTitle("LowADC")
            frame.GetYaxis().SetTitleOffset(2.1)
            frame.Draw()
            gr.Draw("P SAME")
            ffit, slope = _fit_hl_ratio_graph(gr, "%s_%d" % (name, index))
            gr.Draw("P SAME")
            keep.extend([frame, gr])
            if ffit is not None:
                ffit.Draw("LSAME")
                keep.append(ffit)
                if slope is not None:
                    keep.append(_draw_hl_ratio_slope_label(index, slope))
            continue

        histo = histos[index]
        if histo is None:
            continue

        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(0)
        histo.SetMarkerStyle(0)
        histo.SetLineColor(colour)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)
        histo.SetFillColor(0)
        histo.GetYaxis().SetTitleOffset(1.5)

        xax, yax = histo.GetXaxis(), histo.GetYaxis()
        if name.find("hl") != -1:
            x1, x2 = max(1, xax.GetXmin()), min(2e4, xax.GetXmax())
            y1, y2 = max(1, yax.GetXmin()), min(2e3, yax.GetXmax())
            histo.GetXaxis().SetRangeUser(x1, x2)
            histo.GetYaxis().SetRangeUser(y1, y2)
            histo.GetXaxis().SetTitle("HighADC")
            histo.GetYaxis().SetTitle("LowADC")
        elif name.find("high") != -1:
            x1, x2 = max(1, xax.GetXmin()), min(2e4, xax.GetXmax())
            y1, y2 = max(1, yax.GetXmin()), min(2e4, yax.GetXmax())
            histo.GetXaxis().SetRangeUser(x1, x2)
            histo.GetYaxis().SetRangeUser(y1, y2)
            histo.GetXaxis().SetTitle("Main_HighADC")
            histo.GetYaxis().SetTitle("Back_HighADc")
        else:
            x1, x2 = max(1, xax.GetXmin()), min(2e3, xax.GetXmax())
            y1, y2 = max(1, yax.GetXmin()), min(2e3, yax.GetXmax())
            histo.GetXaxis().SetRangeUser(x1, x2)
            histo.GetYaxis().SetRangeUser(y1, y2)
            histo.GetXaxis().SetTitle("Main_LowADC")
            histo.GetYaxis().SetTitle("Back_LowADC")

        ROOT.gPad.SetRightMargin(0.12)
        histo.Draw("COLZ")
        legend = ROOT.TLegend(0.6, 0.72, 0.92, 0.92)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.SetBorderSize(0)
        legend.SetLineColor(0)
        legend.SetLineStyle(1)
        legend.SetLineWidth(3)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.AddEntry(histo, str(index))
        legend.Draw()
        keep.extend([histo, legend])

    c.Update()
    c.SaveAs("Ratio_" + name + ".png")


def plot_ped(name, do_fit=None):
    if do_fit is None:
        do_fit = '--no-ped-fit' not in sys.argv
    infile = ROOT.TFile.Open("hist_calo.root", "READ")
    histos = []
    fit_info_list = []  # 每个pad的(mean, sigma)
    for i in range(0, 25):
        histo = infile.Get(name + "_" + str(i))
        histo.SetDirectory(0)
        mean, sigma = None, None
        if do_fit:
            fit_result = histo.Fit("gaus", "SQ", "", 500, 1800)
            fit = histo.GetFunction("gaus")
            if fit:
                mean = fit.GetParameter(1)
                sigma = fit.GetParameter(2)
        histos.append(histo)
        fit_info_list.append((mean, sigma))
    infile.Close()

    c = ROOT.TCanvas("canvas", '', 0, 0, 1500, 1500)
    c.Divide(5, 5)
    c.SetLogx(False)
    c.SetLogy(True)
    goodColours = [ROOT.kCyan + 2, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange, ROOT.kBlack]

    pt_keep = []
    for index, histo in enumerate(histos):
        colour = goodColours[0]
        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(0)
        histo.SetMarkerStyle(0)
        histo.SetLineColor(colour)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)
        histo.SetFillColor(0)
        histo.GetXaxis().SetTitle("ADC")
        histo.GetYaxis().SetTitle("Counts")
        histo.GetYaxis().SetTitleOffset(1.5)
        c.cd(index + 1)
        ROOT.gPad.SetLogy(True)
        max_y = histo.GetMaximum()
        min_y = 1.0
        if max_y > 0:
            histo.GetYaxis().SetRangeUser(min_y, max_y * 1.3)
        else:
            histo.GetYaxis().SetRangeUser(0.5, 2)

        if index == 0:
            histo.Draw("H")
        else:
            histo.Draw("H")

        mean, sigma = fit_info_list[index]
        if do_fit:
            pt = ROOT.TPaveText(0.58, 0.72, 0.88, 0.90, "NDC")
            pt.SetName(f"ped_fit_{name}_{index}")
            pt.SetBorderSize(0)
            pt.SetFillStyle(0)
            pt.SetTextFont(42)
            pt.SetTextSize(0.07)
            pt.SetTextColor(ROOT.kBlack)
            pt.AddText(f"Mean = {mean:.0f}" if mean is not None else "Mean = N/A")
            pt.AddText(f"Sigma = {sigma:.2f}" if sigma is not None else "Sigma = N/A")
            pt.Draw()
            pt_keep.append(pt)
    key = "Plat_" + name
    c.SaveAs(key + ".png")
    _write_canvas_to_plot_all(c, key)


def plot_2dsingle(name, leg, txt, r1, r2):
    # Load some histos from the example file
    # (Gaussian limits from 2016 TLA conf)
    infile = ROOT.TFile.Open("hist_calo.root", "READ")
    histos = []
    legendLines = []
    if "14" in name:
        r1 = 3e2
        r2 = 1e5

    mini = 9999.
    maxi = 0.

    histo = infile.Get(name)
    histo.SetDirectory(0)
    histos.append(histo)
    legendLines.append(leg)
    # Close the input file
    inf/ile.Close()

    # Make a canvas to put the plot on.
    # We don't want log axes for this plot,
    # but if you do you can control them here.
    c = ROOT.TCanvas("canvas", '', 0, 0, 1024, 768)
    c.SetLogx(False)
    c.SetLogy(False)
    c.SetGridx(0)
    c.SetGridy(0)

    # Decide what x and y range to use in the display.
    xRange = [-3.2, 3.2]
    yRange = [0.9, 1.1]

    # Decide what colours to use.
    # These ones look decent, but obviously use
    # whatever you like best.
    goodColours = [ROOT.kCyan + 2, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange, ROOT.kBlack]

    # Make a legend.
    # These are the locations of the left side, bottom side, right
    # side, and top, as fractions of the canvas.
    legend = ROOT.TLegend(0.6, 0.72, 0.92, 0.92)
    # Make the text a nice fond, and big enough
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    # A few more formatting things .....
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(3)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    # Draw each histogram.
    # You really shouldn't put two histograms with different
    # x axes on the same plot - I'm only doing it here
    # to show you how to draw multiple plots on the same
    # canvas.
    for histo, line in zip(histos, legendLines):

        index = histos.index(histo)
        colour = goodColours[index]

        # Set up marker to look nice
        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(0)
        histo.SetMarkerStyle(20 + index)

        # Set up line to look nice
        histo.SetLineColor(colour)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)

        # Make sure we don't get a fill
        histo.SetFillColor(0)

        # Label my axes!!
        histo.GetXaxis().SetTitle("ADC")
        histo.GetYaxis().SetTitle("Fraction")
        # Move the label around if you want
        histo.GetYaxis().SetTitleOffset(1.5)

        # Set the limit for the axes
        # r1=histo.GetBinCenter(histo.GetMinimumBin())
        # r2=histo.GetBinCenter(histo.GetMaximumBin())*50.
        # histo.GetXaxis().SetRangeUser(r1,r2)
        # histo.GetYaxis().SetRangeUser(r1,r2)

        if index == 0:
            histo.Draw("COLZ")  # Draw data points (you'll get error bars by default)
        else:
            histo.Draw("H SAME")  # SAME means don't get rid of the previous stuff on the canvas

        # Fill entry into legend
        # "PL" means both the line and point style
        # will show up in the legend.
        legend.AddEntry(histo, line, "PL")

    # Actually draw the legend
    legend.Draw()

    # This is one way to draw text on the plot
    myLatex = ROOT.TLatex()
    myLatex.SetTextColor(ROOT.kBlack)
    myLatex.SetNDC()

    # Put an VLAST-P Internal label
    # I think it has to be Helvetica
    myLatex.SetTextSize(0.05)
    myLatex.SetTextFont(72)
    # These are the x and y coordinates of the bottom left corner of the text
    # as fractions of the canvas
    myLatex.DrawLatex(0.58, 0.68, "VLAST-P")
    # Now we switch back to normal font for the "Internal"
    myLatex.SetTextFont(42)
    myLatex.DrawLatex(0.75, 0.68, txt)

    # Update the canvas
    c.Update()

    # Save the output as a .eps, a .C, and a .root
    c.SaveAs(name + leg + txt + ".png")


def plot_single(name, leg, txt, r1, r2):
    # Load some histos from the example file
    # (Gaussian limits from 2016 TLA conf)
    infile = ROOT.TFile.Open("hist_calo.root", "READ")
    histos = []
    legendLines = []
    if "14" in name:
        r1 = 3e2
        r2 = 1e5

    mini = 9999.
    maxi = 0.

    histo = infile.Get(name)
    histo.SetDirectory(0)
    if (histo.Integral() > 0.):
        histo.Scale(1. / histo.Integral())
    histos.append(histo)
    legendLines.append(leg)
    # Close the input file
    infile.Close()

    # Make a canvas to put the plot on.
    # We don't want log axes for this plot,
    # but if you do you can control them here.
    c = ROOT.TCanvas("canvas", '', 0, 0, 1024, 768)
    c.SetLogx(True)
    c.SetLogy(False)
    c.SetGridx(0)
    c.SetGridy(0)

    # Decide what x and y range to use in the display.
    xRange = [-3.2, 3.2]
    yRange = [0.9, 1.1]

    # Decide what colours to use.
    # These ones look decent, but obviously use
    # whatever you like best.
    goodColours = [ROOT.kCyan + 2, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange, ROOT.kBlack]

    # Make a legend.
    # These are the locations of the left side, bottom side, right
    # side, and top, as fractions of the canvas.
    legend = ROOT.TLegend(0.6, 0.72, 0.92, 0.92)
    # Make the text a nice fond, and big enough
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    # A few more formatting things .....
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(3)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    # Draw each histogram.
    # You really shouldn't put two histograms with different
    # x axes on the same plot - I'm only doing it here
    # to show you how to draw multiple plots on the same
    # canvas.
    for histo, line in zip(histos, legendLines):

        index = histos.index(histo)
        colour = goodColours[index]

        # Set up marker to look nice
        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(0)
        histo.SetMarkerStyle(20 + index)

        # Set up line to look nice
        histo.SetLineColor(colour)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)

        # Make sure we don't get a fill
        histo.SetFillColor(0)

        # Label my axes!!
        histo.GetXaxis().SetTitle("ADC")
        histo.GetYaxis().SetTitle("Fraction")
        # Move the label around if you want
        histo.GetYaxis().SetTitleOffset(1.5)

        # Set the limit for the axes
        # r1=histo.GetBinCenter(histo.GetMinimumBin())
        # r2=histo.GetBinCenter(histo.GetMaximumBin())*50.
        histo.GetXaxis().SetRangeUser(r1, r2)
        histo.GetYaxis().SetRangeUser(0.0, histo.GetMaximum() * 1.4)

        if index == 0:
            histo.Draw("H")  # Draw data points (you'll get error bars by default)
        else:
            histo.Draw("H SAME")  # SAME means don't get rid of the previous stuff on the canvas

        # Fill entry into legend
        # "PL" means both the line and point style
        # will show up in the legend.
        legend.AddEntry(histo, line, "PL")

    # Actually draw the legend
    legend.Draw()

    # This is one way to draw text on the plot
    myLatex = ROOT.TLatex()
    myLatex.SetTextColor(ROOT.kBlack)
    myLatex.SetNDC()

    # Put an VLAST-P Internal label
    # I think it has to be Helvetica
    myLatex.SetTextSize(0.05)
    myLatex.SetTextFont(72)
    # These are the x and y coordinates of the bottom left corner of the text
    # as fractions of the canvas
    myLatex.DrawLatex(0.58, 0.68, "VLAST-P")
    # Now we switch back to normal font for the "Internal"
    myLatex.SetTextFont(42)
    myLatex.DrawLatex(0.75, 0.68, txt)

    # Update the canvas
    c.Update()

    # Save the output as a .eps, a .C, and a .root
    c.SaveAs(name + leg + txt + ".png")


def plot_all():
    # Load some histos from the example file
    # (Gaussian limits from 2016 TLA conf)
    infile = ROOT.TFile.Open("hist_calo.root", "READ")
    histos = []
    legendLines = ["Main HG", "Main LG", "Back HG", "Back LG"]
    hname = ["hpmh", "hpml", "hpbh", "hpbl"]

    mini = 9999.
    maxi = 0.

    for i in hname:
        histo = infile.Get(i)
        histo.SetDirectory(0)
        if (histo.Integral() > 0.):
            histo.Scale(1. / histo.Integral())
        histos.append(histo)
    # Close the input file
    infile.Close()

    # Make a canvas to put the plot on.
    # We don't want log axes for this plot,
    # but if you do you can control them here.
    c = ROOT.TCanvas("canvas", '', 0, 0, 1024, 768)
    c.SetLogx(True)
    c.SetLogy(True)
    c.SetGridx(0)
    c.SetGridy(0)

    # Decide what x and y range to use in the display.
    xRange = [-3.2, 3.2]
    yRange = [0.9, 1.1]

    # Decide what colours to use.
    # These ones look decent, but obviously use
    # whatever you like best.
    goodColours = [ROOT.kCyan + 2, ROOT.kBlue + 1, ROOT.kMagenta + 1, ROOT.kOrange, ROOT.kBlack]

    # Make a legend.
    # These are the locations of the left side, bottom side, right
    # side, and top, as fractions of the canvas.
    legend = ROOT.TLegend(0.6, 0.72, 0.92, 0.92)
    # Make the text a nice fond, and big enough
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    # A few more formatting things .....
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(3)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    # Draw each histogram.
    # You really shouldn't put two histograms with different
    # x axes on the same plot - I'm only doing it here
    # to show you how to draw multiple plots on the same
    # canvas.
    for histo, line in zip(histos, legendLines):

        index = histos.index(histo)
        colour = goodColours[index]

        # Set up marker to look nice
        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(1)
        histo.SetMarkerStyle(20 + index)

        # Set up line to look nice
        histo.SetLineColor(colour)
        histo.SetLineWidth(3)
        histo.SetLineStyle(1)

        # Make sure we don't get a fill
        histo.SetFillColor(0)

        # Label my axes!!
        histo.GetXaxis().SetTitle("ADC")
        histo.GetYaxis().SetTitle("Fraction")
        # Move the label around if you want
        histo.GetYaxis().SetTitleOffset(1.5)

        # Set the limit for the axes
        # histo.GetXaxis().SetLimits(xRange[0],xRange[1])
        histo.GetYaxis().SetRangeUser(0.001, 0.15)

        if index == 0:
            histo.Draw("H")  # Draw data points (you'll get error bars by default)
        else:
            histo.Draw("H SAME")  # SAME means don't get rid of the previous stuff on the canvas

        # Fill entry into legend
        # "PL" means both the line and point style
        # will show up in the legend.
        legend.AddEntry(histo, line, "PL")

    # Actually draw the legend
    legend.Draw()

    # This is one way to draw text on the plot
    myLatex = ROOT.TLatex()
    myLatex.SetTextColor(ROOT.kBlack)
    myLatex.SetNDC()

    # Put an VLAST-P Internal label
    # I think it has to be Helvetica
    myLatex.SetTextSize(0.05)
    myLatex.SetTextFont(72)
    # These are the x and y coordinates of the bottom left corner of the text
    # as fractions of the canvas
    myLatex.DrawLatex(0.18, 0.88, "VLAST-P")
    # Now we switch back to normal font for the "Internal"
    myLatex.SetTextFont(42)
    myLatex.DrawLatex(0.35, 0.88, "Pedestal")

    # Update the canvas
    c.Update()

    # Save the output as a .eps, a .C, and a .root
    c.SaveAs("Plat.png")



if __name__ == "__main__":
    _do_sig_fit = SIG_FIT_ENABLED and ("--no-sig-fit" not in sys.argv)
    _do_ped_fit = "--no-ped-fit" not in sys.argv
    if "--print-mpv" in sys.argv:
        idx = sys.argv.index("--print-mpv")
        name = sys.argv[idx + 1] if idx + 1 < len(sys.argv) else "hbh"
        hist_path = "hist_calo.root"
        if "--hist" in sys.argv:
            hidx = sys.argv.index("--hist")
            if hidx + 1 < len(sys.argv):
                hist_path = sys.argv[hidx + 1]
        print("# crystal  MPV(ADC)  [%s from %s]" % (name, hist_path))
        for i, mpv in get_mpv_list(name, hist_path):
            print(i, mpv if mpv is not None else "nan")
        sys.exit(0)
    if os.path.isfile("temp_hist.root"):
        plot_ct()
    plot_all()
    plot_ped("ped_mh", do_fit=_do_ped_fit)
    plot_ped("ped_ml", do_fit=_do_ped_fit)
    plot_ped("ped_bh", do_fit=_do_ped_fit)
    plot_ped("ped_bl", do_fit=_do_ped_fit)
    plot_sig("hmh", do_fit=_do_sig_fit)
    plot_sig("hml", do_fit=_do_sig_fit)
    plot_sig("hbh", do_fit=_do_sig_fit)
    plot_sig("hbl", do_fit=_do_sig_fit)
    plot_ratio("main_hlr")
    plot_ratio("back_hlr")
    # plot_ratio("high_mbr")
    # plot_ratio("low_mbr")

    # chns = [2,12,14,24]
    # chns = [12]
    # for i in chns:
    #     plot_single("hmh_" + str(i), "Main_HG", str(i), 3e2, 5e4)
    #     plot_single("hml_" + str(i), "Main_LG", str(i), 30, 5e3)
    #     plot_single("hbh_" + str(i), "Back_HG", str(i), 3e2, 5e4)
    #     plot_single("hbl_" + str(i), "Back_LG", str(i), 30, 5e3)
    #     plot_2dsingle("main_hlr_" + str(i), "Main_HLRatio", str(i), 1, 1)
    #     plot_2dsingle("back_hlr_" + str(i), "Back_HLRatio", str(i), 1, 1)
    #     plot_2dsingle("high_mbr_" + str(i), "High_MBRatio", str(i), 1, 1)
    #     plot_2dsingle("low_mbr_" + str(i), "Low_MBRatio", str(i), 1, 1)
