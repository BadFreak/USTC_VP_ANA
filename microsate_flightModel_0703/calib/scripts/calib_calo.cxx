#include <iostream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TPad.h"
#include "TStyle.h"

// 局部极大 bin 中心作为峰位（不做高斯拟合）；标定谱随 DAC 单调升高
static const int kMergeBins = 1;

struct CaloPeakFindCfg {
	double peakSeed;
	double minAdcStep;
	double monoTol;
};

static const CaloPeakFindCfg kPeakCfgHigh = {25., 200., 50.};
static const CaloPeakFindCfg kPeakCfgLow = kPeakCfgHigh;
// Calo DAC 输入码（去掉原 0x80 前缀；首点 0x010 而非 0x8010）
static const unsigned int kCalibInputX[40] = {
	// 0x010, 0x013, 0x016, 0x019, 0x01c, 0x01f, 0x022, 0x025, 0x028, 0x02b, 
	// 0x02e, 0x031, 0x034, 0x037, 0x03a, 0x03d, 0x040, 0x05b, 0x076, 0x091,
	0x03, 0x07, 0x0B, 0x0F, 0x13, 0x17, 0x1B, 0x1F, 0x23, 0x27,
    0x2B, 0x2F, 0x33, 0x37, 0x3B, 0x3F, 0x43, 0x5B, 0x76, 0x91,
	0x0ac, 0x0c7, 0x0e2, 0x0fd, 0x118, 0x133, 0x14e, 0x169, 0x184, 0x19f, 
	0x1ba, 0x1d5, 0x1f0, 0x20b, 0x226, 0x241, 0x25c, 0x277, 0x292, 0x2ad
};

static std::vector<float> findPeak(TH1D* h, const CaloPeakFindCfg& cfg, const char* label = "") {
	std::vector<float> ret;
	if (!h)
		return ret;
	const int nx = h->GetNbinsX();

	std::vector<int> cand;
	for (int ib = 1; ib <= nx; ++ib) {
		const double c = h->GetBinContent(ib);
		if (c < cfg.peakSeed)
			continue;
		const double cl = (ib > 1) ? h->GetBinContent(ib - 1) : 0.;
		const double cr = (ib < nx) ? h->GetBinContent(ib + 1) : 0.;
		const double cm = std::max({cl, c, cr});
		if (c == cm)
			cand.push_back(ib);
	}

	std::vector<int> merged;
	for (size_t k = 0; k < cand.size(); ++k) {
		if (merged.empty()) {
			merged.push_back(cand[k]);
			continue;
		}
		int last = merged.back();
		if (cand[k] - last <= kMergeBins) {
			if (h->GetBinContent(cand[k]) > h->GetBinContent(last))
				merged.back() = cand[k];
		} else {
			merged.push_back(cand[k]);
		}
	}

	std::vector<int> monoMerged;
	for (size_t k = 0; k < merged.size(); ++k) {
		const int ib = merged[k];
		const double xc = h->GetBinCenter(ib);
		if (!monoMerged.empty()) {
			const double lastXc = h->GetBinCenter(monoMerged.back());
			if (xc < lastXc + cfg.minAdcStep)
				continue;
		}
		monoMerged.push_back(ib);
	}
	merged = monoMerged;

	for (size_t k = 0; k < merged.size(); ++k) {
		const int ib = merged[k];
		const double peakAdc = h->GetBinCenter(ib);
		const double peakCnt = h->GetBinContent(ib);
		std::cout << label << " peak#" << k << ": adc=" << peakAdc << " cnt=" << peakCnt
		          << std::endl;
		ret.emplace_back(static_cast<float>(peakAdc));
	}

	std::cout << label << " [" << (h->GetName() ? h->GetName() : "?") << "] integral="
	          << h->Integral() << " entries=" << h->GetEntries() << " n_peaks=" << ret.size()
	          << std::endl;
	return ret;
}

static std::vector<float> findPeakHigh(TH1D* h, const char* label = "") {
	return findPeak(h, kPeakCfgHigh, label);
}

static std::vector<float> findPeakLow(TH1D* h, const char* label = "") {
	return findPeak(h, kPeakCfgLow, label);
}

static std::vector<float> selectMonotonicPeaks(const std::vector<float>& peaks, int maxN,
                                               double monoTol) {
	std::vector<float> out;
	out.reserve(peaks.size());
	for (float p : peaks) {
		if (!out.empty() && p < out.back() - monoTol)
			continue;
		out.push_back(p);
		if (static_cast<int>(out.size()) >= maxN)
			break;
	}
	return out;
}

/** 5×5 每格一个晶体：标定峰位 vs DAC → Calib*.png；xRangeMin<xRangeMax 时固定横轴 */
static void saveCalibCanvas(TGraph* g[25], const char* pngPath,
                            double xRangeMin = -1., double xRangeMax = -1.) {
	TCanvas c("calib_canvas", "", 1500, 1500);
	c.Divide(5, 5);
	gStyle->SetOptStat(0);
	const bool fixX = xRangeMin < xRangeMax;
	for (int i = 0; i < 25; ++i) {
		c.cd(i + 1);
		TGraph* gr = g[i];
		if (!gr || gr->GetN() < 1)
			continue;
		gr->SetMarkerColor(kCyan + 2);
		gr->SetMarkerSize(0.5);
		gr->SetMarkerStyle(8);
		gr->SetLineColor(kCyan + 2);
		gr->SetLineWidth(1);
		gr->SetFillColor(0);
		gr->GetXaxis()->SetTitle("DAC value");
		gr->GetYaxis()->SetTitle("Output ADC");
		gr->GetYaxis()->SetTitleOffset(1.5);
		if (fixX)
			gr->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
		gr->SetMinimum(0.);
		gr->SetMaximum(66000.);
		gPad->SetLogx(false);
		gPad->SetLogy(false);
		gr->Draw("APL");
	}
	c.SaveAs(pngPath);
}

static void fillPeakGraph(TGraph* gr, const char* name, const std::vector<float>& peaks,
                          double monoTol) {
	gr->SetName(name);
	if (peaks.empty()) {
		return;
	}

	const int nTarget = 40;
	const float y_fill = *std::max_element(peaks.begin(), peaks.end());
	double prevY = 0.;
	for (int i = 0; i < nTarget; ++i) {
		const double x = static_cast<double>(kCalibInputX[i]);
		double y = static_cast<double>((i < static_cast<int>(peaks.size())) ? peaks[i] : y_fill);
		if (i > 0 && y < prevY - monoTol)
			y = prevY;
		prevY = y;
		gr->SetPoint(i, x, y);
	}
}

void calib_calo(){
	TGraph *gh[25], *gl[25], *gmh[25], *gml[25];
	auto f = TFile::Open("hist.root","READ");
	if (!f || f->IsZombie()) {
		std::cerr << "calib: cannot open hist.root" << std::endl;
		return;
	}
	for (int i = 0; i < 25; ++i) {
		gh[i] = new TGraph();
		gl[i] = new TGraph();
		gmh[i] = new TGraph();
		gml[i] = new TGraph();
		auto hbh = (TH1D*)f->Get(Form("hbh_%d", i));
		auto hbl = (TH1D*)f->Get(Form("hbl_%d", i));
		auto hmh = (TH1D*)f->Get(Form("hmh_%d", i));
		auto hml = (TH1D*)f->Get(Form("hml_%d", i));
		if (hbh) {
			TString tagh = TString::Format("cry%d hbh HG", i);
			fillPeakGraph(gh[i], Form("gh_%d", i),
			              selectMonotonicPeaks(findPeakHigh(hbh, tagh.Data()), 40,
			                                   kPeakCfgHigh.monoTol),
			              kPeakCfgHigh.monoTol);
		} else {
			gh[i]->SetName(Form("gh_%d", i));
		}
		if (hbl) {
			TString tagl = TString::Format("cry%d hbl LG", i);
			fillPeakGraph(gl[i], Form("gl_%d", i),
			              selectMonotonicPeaks(findPeakLow(hbl, tagl.Data()), 40,
			                                   kPeakCfgLow.monoTol),
			              kPeakCfgLow.monoTol);
		} else {
			gl[i]->SetName(Form("gl_%d", i));
		}
		if (hmh) {
			TString tagmh = TString::Format("cry%d hmh HG", i);
			fillPeakGraph(gmh[i], Form("gmh_%d", i),
			              selectMonotonicPeaks(findPeakHigh(hmh, tagmh.Data()), 40,
			                                   kPeakCfgHigh.monoTol),
			              kPeakCfgHigh.monoTol);
		} else {
			gmh[i]->SetName(Form("gmh_%d", i));
		}
		if (hml) {
			TString tagml = TString::Format("cry%d hml LG", i);
			fillPeakGraph(gml[i], Form("gml_%d", i),
			              selectMonotonicPeaks(findPeakLow(hml, tagml.Data()), 40,
			                                   kPeakCfgLow.monoTol),
			              kPeakCfgLow.monoTol);
		} else {
			gml[i]->SetName(Form("gml_%d", i));
		}
	}
	f->Close();
	std::cout << "Peaks found (hbh/hbl/hmh/hml)" << std::endl;
	auto fout = new TFile("peaks_calo.root", "RECREATE");
	fout->cd();
	for (int i = 0; i < 25; ++i) {
		gh[i]->Write();
		gl[i]->Write();
		gmh[i]->Write();
		gml[i]->Write();
	}
	fout->Close();

	saveCalibCanvas(gh, "Calibh.png", 0., 100.);
	saveCalibCanvas(gl, "Calibl.png");
	saveCalibCanvas(gmh, "Calibmh.png", 0., 100.);
	saveCalibCanvas(gml, "Calibml.png");
}
