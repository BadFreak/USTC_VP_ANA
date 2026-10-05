// 从 hist_calo.root 读取单通道直方图，做 Landau×Gauss 拟合并出图
// 画图风格参考 cosmic_ped/scripts/VPlot_calo.py::plot_sig
//
// 用法:
//   root -l -b -q 'draw_one_mip.cxx("hist_calo.root")'
//   root -l -b -q 'draw_one_mip.cxx("hist_calo.root", "hmh_13", "hmh_13.png")'

#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iostream>

// 与 VPlot_calo.py::_langaufun 一致
static Double_t langaufun(Double_t* x, Double_t* par)
{
	if (par[3] <= 0.)
		return 0.;
	constexpr Double_t invsq2pi = 0.3989422804014;
	constexpr Double_t mpshift = -0.22278298;
	constexpr Double_t np = 100.;
	constexpr Double_t sc = 5.;
	const Double_t mpc = par[1] - mpshift * par[0];
	const Double_t xlow = x[0] - sc * par[3];
	const Double_t xupp = x[0] + sc * par[3];
	const Double_t step = (xupp - xlow) / np;
	Double_t sum = 0.;
	for (int i = 1; i <= int(np / 2); ++i) {
		Double_t xx = xlow + (i - 0.5) * step;
		Double_t fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);
		xx = xupp - (i - 0.5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);
	}
	return par[2] * step * sum * invsq2pi / par[3];
}

struct FitRange {
	double fit_xlo;
	double fit_xhi;
	double sv_width;
	double sv_sig;
	double pllo_mp;
	double plhi_mp;
	double pllo_w;
	double plhi_w;
	double pllo_s;
	double plhi_s;
};

static FitRange fitRangeForHist(const char* hname)
{
	FitRange r{};
	const TString name(hname);
	if (name.BeginsWith("hmh") || name.BeginsWith("hbh")) {
		r.fit_xlo = 5500.;
		r.fit_xhi = 30000.;
		r.sv_width = 800.;
		r.sv_sig = 400.;
		r.pllo_mp = 5500.;
		r.plhi_mp = 20000.;
		r.pllo_w = 100.;
		r.plhi_w = 3000.;
		r.pllo_s = 50.;
		r.plhi_s = 2000.;
	} else {
		r.fit_xlo = 600.;
		r.fit_xhi = 3000.;
		r.sv_width = 150.;
		r.sv_sig = 80.;
		r.pllo_mp = 600.;
		r.plhi_mp = 2000.;
		r.pllo_w = 30.;
		r.plhi_w = 500.;
		r.pllo_s = 10.;
		r.plhi_s = 300.;
	}
	return r;
}

constexpr double kPlotXmin = 5000.;
constexpr double kPlotXmax = 30000.;
constexpr double kPlotYmin = 0.;
constexpr double kPlotYmax = 120.;

static void styleHistLikeVPlot(TH1D* h)
{
	const Color_t colour = kCyan + 2;
	h->SetMarkerColor(colour);
	h->SetMarkerSize(0);
	h->SetMarkerStyle(0);
	h->SetLineColor(colour);
	h->SetLineWidth(3);
	h->SetLineStyle(1);
	h->SetFillColor(0);
}

static void stylePlotAxes(TAxis* xax, TAxis* yax)
{
	xax->SetTitle("ADC");
	xax->CenterTitle(true);
	xax->SetTitleSize(0.045);
	xax->SetLabelSize(0.042);
	xax->SetNdivisions(505, kFALSE);

	yax->SetTitle("Entries");
	yax->SetTitleOffset(1.2);
	yax->CenterTitle(true);
	yax->SetTitleSize(0.045);
	yax->SetLabelSize(0.042);
}

static bool FitLandauGauss(TH1D* h, const FitRange& fr, TF1*& ffit, double& mpv)
{
	if (!h || h->Integral() <= 0)
		return false;

	const int nbins = h->GetNbinsX();
	int blo = std::max(1, std::min(h->GetXaxis()->FindBin(fr.fit_xlo), nbins));
	int bhi = std::max(1, std::min(h->GetXaxis()->FindBin(fr.fit_xhi), nbins));
	if (blo > bhi)
		std::swap(blo, bhi);

	const double integral_fit = h->Integral(blo, bhi);
	if (integral_fit <= 0)
		return false;

	int imax = blo;
	for (int b = blo; b <= bhi; ++b) {
		if (h->GetBinContent(b) > h->GetBinContent(imax))
			imax = b;
	}
	double mp_init = h->GetBinCenter(imax);
	mp_init = std::max(fr.pllo_mp, std::min(fr.plhi_mp, mp_init));

	ffit = new TF1(Form("langau_%s", h->GetName()), langaufun, fr.fit_xlo, fr.fit_xhi, 4);
	ffit->SetParameters(fr.sv_width, mp_init, integral_fit, fr.sv_sig);
	ffit->SetParNames("Width", "MP", "Area", "GSigma");
	ffit->SetParLimits(0, fr.pllo_w, fr.plhi_w);
	ffit->SetParLimits(1, fr.pllo_mp, fr.plhi_mp);
	ffit->SetParLimits(2, integral_fit * 0.01, 1e10);
	ffit->SetParLimits(3, fr.pllo_s, fr.plhi_s);

	const TFitResultPtr res = h->Fit(ffit, "RB0QS");
	mpv = ffit->GetParameter(1);
	bool ok = (mpv >= fr.pllo_mp && mpv <= fr.plhi_mp);
	if (res.Get() && res->Status() != 0 && res->Status() != 1 && res->Status() != 3)
		ok = false;
	return ok;
}

void draw_one_mip(const char* histfile = "hist_calo.root",
                  const char* histname = "hmh_13",
                  const char* outpng = "hmh_13.png")
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto* fin = TFile::Open(histfile, "READ");
	if (!fin || fin->IsZombie()) {
		std::cerr << "cannot open " << histfile << std::endl;
		return;
	}

	auto* hist = dynamic_cast<TH1D*>(fin->Get(histname));
	if (!hist) {
		std::cerr << "histogram not found: " << histname << " in " << histfile << std::endl;
		fin->Close();
		return;
	}

	auto* h = dynamic_cast<TH1D*>(hist->Clone(Form("%s_draw", histname)));
	h->SetDirectory(0);
	fin->Close();

	const FitRange fr = fitRangeForHist(histname);
	styleHistLikeVPlot(h);

	auto* c = new TCanvas("canvas", "", 0, 0, 1024, 768);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	gPad->SetLogx(false);
	gPad->SetLogy(false);
	gPad->SetGridx(true);
	gPad->SetGridy(true);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFrameLineWidth(2);

	// DrawFrame 可精确设定坐标范围；SetRangeUser 会贴到直方图 bin 边界（对数分 bin 时可能是 4955）
	gPad->DrawFrame(kPlotXmin, kPlotYmin, kPlotXmax, kPlotYmax);
	if (auto* frame = dynamic_cast<TH1F*>(gPad->GetPrimitive("hframe"))) {
		frame->SetLineWidth(2);
		stylePlotAxes(frame->GetXaxis(), frame->GetYaxis());
	}
	h->Draw("HIST SAME");

	TF1* ffit = nullptr;
	double mpv = 0.;
	if (FitLandauGauss(h, fr, ffit, mpv)) {
		auto* ffit_draw = dynamic_cast<TF1*>(ffit->Clone(Form("langau_norm_%s", histname)));
		ffit_draw->SetLineColor(kRed);
		ffit_draw->SetLineWidth(2);
		ffit_draw->Draw("LSAME");

		auto* pt = new TPaveText(0.58, 0.72, 0.88, 0.82, "NDC");
		pt->SetBorderSize(0);
		pt->SetFillStyle(0);
		pt->SetTextFont(42);
		pt->SetTextSize(0.045);
		pt->SetTextColor(kBlack);
		pt->AddText(Form("#bf{#it{MPV: %.0f}}", mpv));
		pt->Draw();
		std::cout << histname << " MPV = " << mpv << std::endl;
	} else {
		std::cerr << "Landau#timesGauss fit failed for " << histname << std::endl;
	}

	c->Update();
	c->SaveAs(outpng);
	std::cout << "saved " << outpng << std::endl;
}
