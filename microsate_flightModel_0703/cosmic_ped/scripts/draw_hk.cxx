#include <set>
#include <cmath>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TAxis.h"
#include "ROOT/RDataFrame.hxx"
void draw_temp(const std::string& fname){
	TFile *fin = TFile::Open(fname.c_str(), "READ");
	if (!fin || !fin->Get("hkTree")) {
		if (fin) fin->Close();
		std::cout << " draw_temp: no hkTree in " << fname << ", skip." << std::endl;
		return;
	}
	fin->Close();

	auto fg=[](const TString& name,const TString& title)->TGraph*{
		TGraph *g=new TGraph();
		g->SetName(name);
		g->SetTitle(title);
		return g;
	};
	// 遥测 FEE：2=CsITK，5–8=Calo 四块板（与 hkTree TelFEE / m_statusFeID 一致）
	static const int kTelFee[5] = {2, 5, 6, 7, 8};
	auto slotTel = [](int tel)->int {
		for (int s = 0; s < 5; s++) {
			if (kTelFee[s] == tel) return s;
		}
		return -1;
	};
	TGraph *gc[5][3], *gt[5][4];
	for (int s = 0; s < 5; s++) {
		int fid = kTelFee[s];
		for (int ic = 0; ic < 3; ic++) {
			gc[s][ic] = fg(TString::Format("f%dc%d", fid, ic),
				       TString::Format("TelFEE %d Current_%d", fid, ic));
		}
		for (int it = 0; it < 4; it++) {
			gt[s][it] = fg(TString::Format("f%dt%d", fid, it),
				       TString::Format("TelFEE %d Temperature_%d", fid, it));
		}
	}
	int nc[5][3] = {};
	int nt[5][4] = {};
	int seq = 0;
	ROOT::RDataFrame df("hkTree", fname);
	df.Foreach([&](const int& telFee,
			const float& c0, const float& c1, const float& c2,
			const float& t0, const float& t1, const float& t2, const float& t3){
		double x = static_cast<double>(seq++);
		int s = slotTel(telFee);
		if (s < 0) return;
		if (std::isfinite(c0)) { gc[s][0]->SetPoint(nc[s][0], x, c0); nc[s][0]++; }
		if (std::isfinite(c1)) { gc[s][1]->SetPoint(nc[s][1], x, c1); nc[s][1]++; }
		if (std::isfinite(c2)) { gc[s][2]->SetPoint(nc[s][2], x, c2); nc[s][2]++; }
		if (std::isfinite(t0)) { gt[s][0]->SetPoint(nt[s][0], x, t0); nt[s][0]++; }
		if (std::isfinite(t1)) { gt[s][1]->SetPoint(nt[s][1], x, t1); nt[s][1]++; }
		if (std::isfinite(t2)) { gt[s][2]->SetPoint(nt[s][2], x, t2); nt[s][2]++; }
		if (std::isfinite(t3)) { gt[s][3]->SetPoint(nt[s][3], x, t3); nt[s][3]++; }
	},{"FEEID","C0","C1","C2","T0","T1","T2","T3"});
	for (int s = 0; s < 5; s++) {
		int fid = kTelFee[s];
		for (int ic = 0; ic < 3; ic++) {
			gc[s][ic]->SetTitle(TString::Format("TelFEE %d current ch %d", fid, ic));
			gc[s][ic]->GetXaxis()->SetTitle("hkTree row index (decode order)");
			gc[s][ic]->GetYaxis()->SetTitle("I [mA]");
			gc[s][ic]->SetMarkerStyle(20 + ic);
			gc[s][ic]->SetMarkerSize(0.6);
			gc[s][ic]->SetLineColor(2 + ic);
			gc[s][ic]->SetMarkerColor(2 + ic);
		}
		for (int it = 0; it < 4; it++) {
			gt[s][it]->SetTitle(TString::Format("TelFEE %d temperature sensor %d", fid, it));
			gt[s][it]->GetXaxis()->SetTitle("hkTree row index (decode order)");
			gt[s][it]->GetYaxis()->SetTitle("T [^{o}C]");
			gt[s][it]->SetMarkerStyle(20 + it);
			gt[s][it]->SetMarkerSize(0.6);
			gt[s][it]->SetLineColor(1 + it);
			gt[s][it]->SetMarkerColor(1 + it);
		}
	}
	TFile *fout=new TFile("temp_hist.root","RECREATE");
	fout->cd();
	for (int s = 0; s < 5; s++) {
		gc[s][0]->Write();
		gc[s][1]->Write();
		gc[s][2]->Write();
		gt[s][0]->Write();
		gt[s][1]->Write();
		gt[s][2]->Write();
		gt[s][3]->Write();
	}

	// 直接绘图并导出到当前目录（分 pad 显示，不叠加）
	for (int s = 0; s < 5; s++) {
		int fid = kTelFee[s];
		auto cTemp = new TCanvas(TString::Format("c_temp_f%d", fid),
		                         TString::Format("Temperature TelFEE %d", fid),
		                         1200, 700);
		cTemp->Divide(2, 2);
		for (int it = 0; it < 4; it++) {
			cTemp->cd(it + 1);
			gt[s][it]->SetTitle(TString::Format("TelFEE %d T%d;hkTree row index (decode order);T [^{o}C]", fid, it));
			gt[s][it]->Draw("APL");
		}
		cTemp->Write();
		cTemp->SaveAs(TString::Format("temp_f%d.png", fid));

		auto cCur = new TCanvas(TString::Format("c_current_f%d", fid),
		                        TString::Format("Current TelFEE %d", fid),
		                        1300, 400);
		cCur->Divide(3, 1);
		for (int ic = 0; ic < 3; ic++) {
			cCur->cd(ic + 1);
			gc[s][ic]->SetTitle(TString::Format("TelFEE %d C%d;hkTree row index (decode order);I [mA]", fid, ic));
			gc[s][ic]->Draw("APL");
		}
		cCur->Write();
		cCur->SaveAs(TString::Format("current_f%d.png", fid));
	}
	fout->Close();
}
// nCrystalHit: 只保留击中该个数晶体的事例（<0 不筛选）；hitThreshold: 击中判据 ADC-pedestal > hitThreshold
// useLogBin: true 用对数分 bin，false 用线性等宽分 bin
void draw_hk(const std::string& filename){
    std::cout << " Draw " << filename << std::endl;
    TString fname = TString(filename);
    draw_temp(filename);
}
