#include <iostream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"

// 与 CsI_inet/calibration 一致：先局部极大再拟合，避免 while(i++) 漏峰
static const double kPeakSeed = 25.;
static const int kMergeBins = 1;

std::vector<float> findPeak(TH1D* h, const char* label = ""){
	std::vector<float> ret;
	if (!h) return ret;
	const int nx = h->GetNbinsX();

	std::vector<int> cand;
	for (int ib = 1; ib <= nx; ++ib) {
		const double c = h->GetBinContent(ib);
		if (c < kPeakSeed) continue;
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

	for (size_t k = 0; k < merged.size(); ++k) {
		const int ib = merged[k];
		const double mean = h->GetBinCenter(ib);
		const double halfw = 400.;
		TString fname = TString::Format("gaus_cal_%s_%zu", h->GetName(), k);
		auto *f1 = new TF1(fname.Data(), "[0]*TMath::Gaus(x,[1],[2])", mean - halfw, mean + halfw);
		f1->SetParameters(std::max(50., h->GetBinContent(ib)), mean, 80.);
		f1->SetParLimits(2, 8., 300.);
		h->Fit(fname.Data(), "RQN", "", mean - halfw, mean + halfw);
		const double output_mean = f1->GetParameter(1);
		const double output_rms = f1->GetParameter(2);
		const double output_mag = f1->GetParameter(0);
		std::cout << label << " peak#" << k << ": mean=" << output_mean << " rms=" << output_rms
		          << " amp=" << output_mag << std::endl;
		ret.emplace_back(static_cast<float>(output_mean));
		delete f1;
	}

	std::cout << label << " [" << (h->GetName() ? h->GetName() : "?") << "] integral=" << h->Integral()
	          << " entries=" << h->GetEntries() << " n_peaks=" << ret.size() << std::endl;
	return ret;
}

void calib(){
	TGraph *gh[25],*gl[25];
	auto f = TFile::Open("hist.root","READ");
	if (!f || f->IsZombie()) {
		std::cerr << "calib: cannot open hist.root" << std::endl;
		return;
	}
	for(int i=0;i<25;i++){
		gh[i] = new TGraph();
		gh[i]->SetName(Form("gh_%d",i));
		gl[i] = new TGraph();
		gl[i]->SetName(Form("gl_%d",i));
		auto hbh = (TH1D*)f->Get(Form("hbh_%d",i));
		auto hbl = (TH1D*)f->Get(Form("hbl_%d",i));
		if (!hbh || !hbl) continue;
		TString tagh = TString::Format("cry%d hbh", i);
		TString tagl = TString::Format("cry%d hbl", i);
		auto hpeaks = findPeak(hbh, tagh.Data());
		auto lpeaks = findPeak(hbl, tagl.Data());
		int index=0;
		for(auto h:hpeaks){
			gh[i]->SetPoint(index,index+1,h);
			index++;
		}
		index=0;
		for(auto l:lpeaks){
			gl[i]->SetPoint(index,index+1,l);
			index++;
		}
	}
	std::cout<<"Peaks found"<<std::endl;
	auto fout = new TFile("peaks.root","RECREATE");
	fout->cd();
	for(int i=0;i<25;i++){
		gh[i]->Write();
		gl[i]->Write();
	}
	fout->Close();
}
