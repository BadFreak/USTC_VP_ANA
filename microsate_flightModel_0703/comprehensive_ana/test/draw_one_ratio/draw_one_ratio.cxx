// 从 hist_calo.root 读取单通道 HL ratio TGraph（draw_calo 写的 main_hlg_* / back_hlg_*），
// 画散点并做线性拟合
//
// 用法:
//   root -l -b -q 'draw_one_ratio.cxx("hist_calo.root")'
//   root -l -b -q 'draw_one_ratio.cxx("hist_calo.root", "main_hlg_12", "main_hlg_12.png")'
//   # 兼容旧名 main_hlr_12 -> 自动映射为 main_hlg_12

#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>

constexpr double kPlotXmin = 0.;
constexpr double kPlotXmax = 2e4;
constexpr double kPlotYmin = 0.;
constexpr double kPlotYmax = 2e3;
constexpr double kFitXmin = 200.;
constexpr double kFitXmax = 16000.;

struct LinearFitResult {
	TF1* ffit = nullptr;
	double slope = 0.;
	double slopeErr = 0.;
	bool ok = false;
};

static TString graphNameFromArg(const char* name)
{
	TString gname(name);
	if (gname.BeginsWith("main_hlr_"))
		gname.Replace(0, 9, "main_hlg_");
	else if (gname.BeginsWith("back_hlr_"))
		gname.Replace(0, 9, "back_hlg_");
	return gname;
}

static void styleGraphLikeVPlot(TGraph* g)
{
	const Color_t colour = kCyan + 2;
	g->SetMarkerColor(colour);
	g->SetMarkerSize(2);
	g->SetMarkerStyle(3);
	g->SetLineColor(colour);
	g->SetLineWidth(1);
	g->SetFillColor(0);
}

static void stylePlotAxes(TAxis* xax, TAxis* yax)
{
	xax->SetTitle("HighADC");
	xax->CenterTitle(true);
	xax->SetTitleSize(0.045);
	xax->SetLabelSize(0.042);
	xax->SetNdivisions(505, kFALSE);

	yax->SetTitle("LowADC");
	yax->SetTitleOffset(1.2);
	yax->CenterTitle(true);
	yax->SetTitleSize(0.045);
	yax->SetLabelSize(0.042);
}

static void addSlopeLegend(double slope)
{
	auto* pt = new TPaveText(0.18, 0.72, 0.48, 0.82, "NDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextFont(42);
	pt->SetTextSize(0.045);
	pt->SetTextColor(kBlack);
	pt->AddText(Form("#bf{#it{slope: %.3f}}", slope));
	pt->Draw();
}

static void setupRatioPad()
{
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	gPad->SetRightMargin(0.14);
	gPad->SetLogx(false);
	gPad->SetLogy(false);
	gPad->SetGridx(true);
	gPad->SetGridy(true);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFrameLineWidth(2);
}

static LinearFitResult FitLinear(TGraph* g, const char* tag)
{
	LinearFitResult out{};
	if (!g || g->GetN() <= 0)
		return out;

	out.ffit = new TF1(Form("pol1_%s", tag), "pol1", kFitXmin, kFitXmax);
	out.ffit->SetLineColor(kRed);
	out.ffit->SetLineWidth(3);

	const TFitResultPtr res = g->Fit(out.ffit, "RB0QS");
	out.slope = out.ffit->GetParameter(1);
	out.slopeErr = out.ffit->GetParError(1);
	out.ok = (res.Get() && res->Status() == 0);
	if (!out.ok && res.Get())
		out.ok = (res->Status() == 1 || res->Status() == 3);
	return out;
}

void draw_one_ratio(const char* histfile = "hist_calo.root",
                    const char* graphname = "main_hlg_12",
                    const char* outpng = "main_hlg_12.png")
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto* fin = TFile::Open(histfile, "READ");
	if (!fin || fin->IsZombie()) {
		std::cerr << "cannot open " << histfile << std::endl;
		return;
	}

	const TString gname = graphNameFromArg(graphname);
	auto* graph = dynamic_cast<TGraph*>(fin->Get(gname.Data()));
	if (!graph) {
		std::cerr << "TGraph not found: " << gname << " in " << histfile
		          << " (need rerun draw_calo to create main_hlg_* / back_hlg_*)" << std::endl;
		fin->Close();
		return;
	}

	auto* g = dynamic_cast<TGraph*>(graph->Clone(Form("%s_draw", gname.Data())));
	fin->Close();

	styleGraphLikeVPlot(g);

	const LinearFitResult fit = FitLinear(g, gname.Data());
	if (fit.ok)
		std::cout << gname << " n=" << g->GetN() << " slope = " << fit.slope << " +/- " << fit.slopeErr
		          << std::endl;
	else
		std::cerr << "linear fit failed for " << gname << std::endl;

	auto* c = new TCanvas("canvas_graph", "", 0, 0, 1024, 768);
	setupRatioPad();

	gPad->DrawFrame(kPlotXmin, kPlotYmin, kPlotXmax, kPlotYmax);
	if (auto* frame = dynamic_cast<TH1F*>(gPad->GetPrimitive("hframe"))) {
		frame->SetLineWidth(2);
		stylePlotAxes(frame->GetXaxis(), frame->GetYaxis());
	}

	g->Draw("P SAME");
	if (fit.ok) {
		fit.ffit->Draw("LSAME");
		addSlopeLegend(fit.slope);
	}

	c->Update();
	c->SaveAs(outpng);
	std::cout << "saved " << outpng << std::endl;
}
