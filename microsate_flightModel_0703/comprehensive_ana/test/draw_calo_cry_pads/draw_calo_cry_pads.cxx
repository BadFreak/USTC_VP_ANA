// 每个晶体 2×2：Main/Back × HG/LG；(CellADC-CellPLAT) vs Entry$
// 注意：caloTree 的 CellID/CellADC 是 vector，不能用 TTree::Draw 的 cut 按元素筛选
#include <iostream>
#include <vector>
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

static int channelPadIndex(int cellid) {
	const int mbid = (cellid % 10000) / 1000;
	const int hg = (cellid % 1000) / 100;
	if (mbid == 1 && hg == 1)
		return 0; // Main HG
	if (mbid == 1 && hg == 0)
		return 1; // Main LG
	if (mbid == 0 && hg == 1)
		return 2; // Back HG
	if (mbid == 0 && hg == 0)
		return 3; // Back LG
	return -1;
}

void draw_calo_cry_pads(const char* infile = "result_VLAST_data_260526.root",
                        const char* outdir = "calo_cry_pads") {
	auto* fin = TFile::Open(infile, "READ");
	if (!fin || fin->IsZombie()) {
		std::cerr << "cannot open " << infile << std::endl;
		return;
	}
	auto* tree = dynamic_cast<TTree*>(fin->Get("caloTree"));
	if (!tree) {
		std::cerr << "no caloTree in " << infile << std::endl;
		fin->Close();
		return;
	}

	std::vector<int>* cellId = nullptr;
	std::vector<int>* cellAdc = nullptr;
	std::vector<int>* cellPlat = nullptr;
	tree->SetBranchAddress("CellID", &cellId);
	tree->SetBranchAddress("CellADC", &cellAdc);
	tree->SetBranchAddress("CellPLAT", &cellPlat);

	const Long64_t nent = tree->GetEntries();
	std::vector<double> xs[25][4];
	std::vector<double> ys[25][4];
	for (int c = 0; c < 25; ++c)
		for (int p = 0; p < 4; ++p) {
			xs[c][p].reserve(static_cast<size_t>(nent / 10));
			ys[c][p].reserve(static_cast<size_t>(nent / 10));
		}

	for (Long64_t e = 0; e < nent; ++e) {
		tree->GetEntry(e);
		if (!cellId || !cellAdc || !cellPlat || cellId->size() != cellAdc->size() ||
		    cellId->size() != cellPlat->size())
			continue;
		for (size_t i = 0; i < cellId->size(); ++i) {
			const int id = cellId->at(i);
			const int cry = id / 100000;
			if (cry < 1 || cry > 25)
				continue;
			const int pad = channelPadIndex(id);
			if (pad < 0)
				continue;
			xs[cry - 1][pad].push_back(static_cast<double>(e));
			ys[cry - 1][pad].push_back(static_cast<double>(cellAdc->at(i) - cellPlat->at(i)));
		}
	}

	gSystem->mkdir(outdir, true);
	gStyle->SetOptStat(0);

	static const char* kPadTitle[] = {"Main HG", "Main LG", "Back HG", "Back LG"};
	static const Color_t kPadColor[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};

	for (int cry = 1; cry <= 25; ++cry) {
		TCanvas c(Form("Cry%d", cry), Form("Crystal %d", cry), 1400, 1200);
		c.Divide(2, 2);
		for (int p = 0; p < 4; ++p) {
			c.cd(p + 1);
			gPad->SetGrid();
			const auto& xv = xs[cry - 1][p];
			const auto& yv = ys[cry - 1][p];
			std::cout << "crystal: " << cry << " " << kPadTitle[p] << ": " << yv.size() << " points\n";
			int ycount = 0;
			std::for_each(yv.begin(), yv.end(), [&ycount](double y) {
				std::cout << y << " ";
				++ycount;
				if(ycount%800==0) std::cout << std::endl;
			});
			std::cout << std::endl;
			if (xv.empty()) {
				std::cout << "Cry" << cry << " " << kPadTitle[p] << ": no points\n";
				continue;
			}
			auto* gr = new TGraph(static_cast<int>(xv.size()), xv.data(), yv.data());
			gr->SetTitle(Form("%s;Entry;CellADC-CellPLAT", kPadTitle[p]));
			gr->SetLineColor(kPadColor[p]);
			gr->SetMarkerColor(kPadColor[p]);
			gr->SetMarkerStyle(20);
			gr->SetMarkerSize(0.35);
			gr->Draw("AP");
			gr->GetYaxis()->SetRangeUser(0., 70000.);
		}
		const TString outpng = Form("%s/Cry%02d.png", outdir, cry);
		c.SaveAs(outpng);
		std::cout << "saved " << outpng << " (pts: "
		          << xs[cry - 1][0].size() << "," << xs[cry - 1][1].size() << ","
		          << xs[cry - 1][2].size() << "," << xs[cry - 1][3].size() << ")\n";
	}
	fin->Close();
}
