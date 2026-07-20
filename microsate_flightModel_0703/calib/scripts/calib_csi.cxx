#include <iostream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

static const int kNCsiCh = 8;
static const int kNCsiCalibPoints = 20;
// CsI 标定 DAC：0x10 起公差 0x8，共 20 点（与 Calo 40 点表无关）
static const unsigned int kCalibInputXCsi[kNCsiCalibPoints] = {
	// 0x10, 0x18, 0x20, 0x28, 0x30, 0x38, 0x40, 0x48, 0x50, 0x58,
	// 0x60, 0x68, 0x70, 0x78, 0x80, 0x88, 0x90, 0x98, 0xa0, 0xa8

	0x03, 0x0C, 0x15, 0x1E, 0x27, 0x30, 0x39, 0x42, 0x4B, 0x54,
    0x5D, 0x66, 0x6F, 0x78, 0x81, 0x8A, 0x93, 0x9C, 0xA5, 0xAE
};


static const double kPeakSeed = 25.;
static const int kMergeBins = 1;
static const double kMinPeakAdcStep = 200.;
static const double kMonoPeakTol = 50.;

std::vector<float> findPeak(TH1D* h, const char* label = "") {
	std::vector<float> ret;
	if (!h)
		return ret;
	const int nx = h->GetNbinsX();

	std::vector<int> cand;
	for (int ib = 1; ib <= nx; ++ib) {
		const double c = h->GetBinContent(ib);
		if (c < kPeakSeed)
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
			if (xc < lastXc + kMinPeakAdcStep)
				continue;
		}
		monoMerged.push_back(ib);
	}
	merged = monoMerged;

	for (size_t k = 0; k < merged.size(); ++k) {
		const int ib = merged[k];
		const double mean = h->GetBinCenter(ib);
		const double halfw = 400.;
		TString fname = TString::Format("gaus_cal_%s_%zu", h->GetName(), k);
		auto* f1 = new TF1(fname.Data(), "[0]*TMath::Gaus(x,[1],[2])", mean - halfw, mean + halfw);
		f1->SetParameters(std::max(50., h->GetBinContent(ib)), mean, 80.);
		f1->SetParLimits(2, 8., 300.);
		h->Fit(fname.Data(), "RQN", "", mean - halfw, mean + halfw);
		std::cout << label << " peak#" << k << ": mean=" << f1->GetParameter(1)
		          << " rms=" << f1->GetParameter(2) << " amp=" << f1->GetParameter(0) << std::endl;
		ret.emplace_back(static_cast<float>(f1->GetParameter(1)));
		delete f1;
	}

	std::cout << label << " [" << (h->GetName() ? h->GetName() : "?") << "] integral=" << h->Integral()
	          << " entries=" << h->GetEntries() << " n_peaks=" << ret.size() << std::endl;
	return ret;
}

static std::vector<float> selectMonotonicPeaks(const std::vector<float>& peaks, int maxN) {
	std::vector<float> out;
	for (float p : peaks) {
		if (!out.empty() && p < out.back() - kMonoPeakTol)
			continue;
		out.push_back(p);
		if (static_cast<int>(out.size()) >= maxN)
			break;
	}
	return out;
}

static void fillPeakGraphCsi(TGraph* gr, const char* name, const std::vector<float>& peaks) {
	gr->SetName(name);
	if (peaks.empty())
		return;
	const float y_fill = *std::max_element(peaks.begin(), peaks.end());
	double prevY = 0.;
	for (int i = 0; i < kNCsiCalibPoints; ++i) {
		const double x = static_cast<double>(kCalibInputXCsi[i]);
		double y = static_cast<double>((i < static_cast<int>(peaks.size())) ? peaks[i] : y_fill);
		if (i > 0 && y < prevY - kMonoPeakTol)
			y = prevY;
		prevY = y;
		gr->SetPoint(i, x, y);
	}
	if (static_cast<int>(peaks.size()) != kNCsiCalibPoints)
		std::cout << name << ": expected " << kNCsiCalibPoints << " peaks, got " << peaks.size()
		          << std::endl;
}

static void saveCalibCanvas8(TGraph* g[kNCsiCh], const char* pngPath) {
	TCanvas c("calib_csi_canvas", "", 1200, 600);
	c.Divide(4, 2);
	gStyle->SetOptStat(0);
	for (int i = 0; i < kNCsiCh; ++i) {
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
		gr->SetMinimum(0.);
		gr->SetMaximum(66000.);
		gPad->SetLogx(false);
		gPad->SetLogy(false);
		gr->Draw("APL");
	}
	c.SaveAs(pngPath);
}

void calib_csi() {
	TGraph* g[kNCsiCh];
	auto f = TFile::Open("hist.root", "READ");
	if (!f || f->IsZombie()) {
		std::cerr << "calib_csi: cannot open hist.root" << std::endl;
		return;
	}
	for (int i = 0; i < kNCsiCh; ++i) {
		g[i] = new TGraph();
		auto h = (TH1D*)f->Get(Form("hch_%d", i));
		if (h) {
			TString tag = TString::Format("ch%d", i);
			fillPeakGraphCsi(g[i], Form("g_%d", i),
			                 selectMonotonicPeaks(findPeak(h, tag.Data()), kNCsiCalibPoints));
		} else {
			g[i]->SetName(Form("g_%d", i));
		}
	}
	f->Close();

	std::cout << "CsI peaks found (8 ch, same mode)" << std::endl;
	auto fout = new TFile("peaks_csi.root", "RECREATE");
	fout->cd();
	for (int i = 0; i < kNCsiCh; ++i)
		g[i]->Write();
	fout->Close();

	saveCalibCanvas8(g, "Calib.png");
}
