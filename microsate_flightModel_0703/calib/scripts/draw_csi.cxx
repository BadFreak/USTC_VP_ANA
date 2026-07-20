// CsITK 8 路，无高低增益/正反面区分；CellID = (chn+1)*100000
void draw_ped(std::string fname, const float& dataFraction) {
	const int kNCh = 8;
	const int n_bins = 1000;
	const int bin_min = 5e3, bin_max = 7e4;

	ROOT::RDataFrame df("csiTree", fname);
	std::unordered_map<int, TH1D*> hch;
	for (int i = 0; i < kNCh; ++i) {
		hch[i] = new TH1D(Form("hch_%d", i), Form("Ch_%d_ADC", i), n_bins, bin_min, bin_max);
	}

	const float n_total = float(*df.Count());
	const int n_actual = int(n_total * dataFraction);
	auto dff = df.Range(n_actual);
	dff.Foreach(
	    [&hch](const std::vector<int>& id, const std::vector<int>& adc) {
		    for (int i = 0; i < int(id.size()); ++i) {
			    const int ch = id.at(i) / 100000 - 1;
			    if (ch < 0 || ch >= kNCh)
				    continue;
			    hch.at(ch)->Fill(adc.at(i));
		    }
	    },
	    {"CellID", "CellADC"});

	// 追加写入，勿 RECREATE：run_one 中 draw_calo 已写入 Calo 直方图
	auto fout = new TFile("hist.root", "UPDATE");
	fout->cd();
	for (int i = 0; i < kNCh; ++i)
		hch[i]->Write();
	fout->Close();
}

void draw_csi(const float& dataFraction = 1., const std::string& filename = "") {
	draw_ped(filename, dataFraction);
}
