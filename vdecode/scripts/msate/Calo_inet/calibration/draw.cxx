void draw_ped(std::string fname,const float& dataFraction){
    int pmin=1e0,pmax=1e5;
    ROOT::RDataFrame df("caloTree",fname);
    TH1D *hpmh = new TH1D("hpmh","hplat",100,pmin,pmax);
    TH1D *hpbh = new TH1D("hpbh","hplat",100,pmin,pmax);
    TH1D *hpml = new TH1D("hpml","hplat",100,pmin,pmax);
    TH1D *hpbl = new TH1D("hpbl","hplat",100,pmin,pmax);
    std::unordered_map<int,TH2D*> umain_ratio; // main,cryid, TH2D* ratio
    std::unordered_map<int,TH2D*> uback_ratio; // back,cryid, TH2D* ratio
    std::unordered_map<int,TH2D*> uhigh_ratio; //main and back
    std::unordered_map<int,TH2D*> ulow_ratio; // main and back
    for(int i=0;i<25;i++){
	    umain_ratio[i] = new TH2D(Form("hmainr_%d",i),Form("Low_High_Ratio_%d",i),100,1,1e3,100,1,1e4);
	    uback_ratio[i] = new TH2D(Form("hbackr_%d",i),Form("Low_High_Ratio_%d",i),100,1,1e3,100,1,1e4);
	    uhigh_ratio[i] = new TH2D(Form("hhighr_%d",i),Form("High_MBRatio_%d",i),100,1,1e4,100,1,1e4);
	    ulow_ratio[i] = new TH2D(Form("hlowr_%d",i),Form("Low_MBRatio_%d",i),100,1,1e3,100,1,1e3);
    }
    auto f = [](const int& feeid)->TGraph*{
	    TGraph *g=new TGraph();
	    g->SetName(TString::Format("gfee%s",std::to_string(feeid).c_str()));
	    g->SetTitle(TString::Format("FEE_%s",std::to_string(feeid).c_str()));
	    return g;
    };
    TGraph *gfee1 = f(1);
    TGraph *gfee2 = f(2);
    TGraph *gfee3 = f(3);
    TGraph *gfee4 = f(4);
    const int n_bins = 1000;
    const int bin_min = 5e0;
    const int bin_max = 7e4;
    //double bins[n_bins + 1]; // 存储 bin 边界
    //double log_min = 3; // 最小值的对数
    //double log_max = 4.7;
    //double log_bin_width = (log_max - log_min) / n_bins; // 对数 bin 宽度
    //for (int i = 0; i <= n_bins; i++) {
    //    bins[i] = TMath::Power(10, log_min + i * log_bin_width);
    //}
    TH1D *hmh = new TH1D("hmh","hE",n_bins,bin_min,bin_max);
    TH1D *hml = new TH1D("hml","hE",n_bins,bin_min,bin_max);
    TH1D *hbh = new TH1D("hbh","hE",n_bins,bin_min,bin_max);
    TH1D *hbl = new TH1D("hbl","hE",n_bins,bin_min,bin_max);
    TH1D *hh = new TH1D("hh","hE",n_bins,bin_min,bin_max);
    TH1D *hl = new TH1D("hl","hE",n_bins,bin_min,bin_max);
    std::unordered_map<int,TH1D*> umap_hpmh;
	std::unordered_map<int,TH1D*> umap_hpml;
	std::unordered_map<int,TH1D*> umap_hpbh;
	std::unordered_map<int,TH1D*> umap_hpbl;
    std::unordered_map<int,TH1D*> umap_hmh;
	std::unordered_map<int,TH1D*> umap_hml;
	std::unordered_map<int,TH1D*> umap_hbh;
	std::unordered_map<int,TH1D*> umap_hbl;
	const int N=25;
	for(int i=0;i<N;i++){
		umap_hpmh[i]=new TH1D(Form("hpmh_%d",i),Form("Pedestal_Main_High_%d",i),1000,pmin,pmax);
		umap_hpml[i]=new TH1D(Form("hpml_%d",i),Form("Pedestal_Main_Low_%d",i),1000,pmin,pmax);
		umap_hpbh[i]=new TH1D(Form("hpbh_%d",i),Form("Pedestal_Back_High_%d",i),1000,pmin,pmax);
		umap_hpbl[i]=new TH1D(Form("hpbl_%d",i),Form("Pedestal_Back_Low_%d",i),1000,pmin,pmax);
		umap_hmh[i]=new TH1D(Form("hmh_%d",i),Form("Main_High_%d",i),n_bins,bin_min,bin_max);
		umap_hml[i]=new TH1D(Form("hml_%d",i),Form("Main_Low_%d",i),n_bins,bin_min,bin_max);
		umap_hbh[i]=new TH1D(Form("hbh_%d",i),Form("Back_High_%d",i),n_bins,bin_min,bin_max);
		umap_hbl[i]=new TH1D(Form("hbl_%d",i),Form("Back_Low_%d",i),n_bins,bin_min,bin_max);
	}
	float N_total = float(*df.Count());
	int N_actual = int(N_total * dataFraction);
	auto dff = df.Range(N_actual);
    dff.Foreach([&uhigh_ratio,&ulow_ratio,&umain_ratio,&uback_ratio,hpmh,hpbh,hpml,hpbl,&umap_hpmh,&umap_hpml,&umap_hpbh,&umap_hpbl,&hmh,&hml,&hbh,&hbl,&hh,&hl,&umap_hmh,&umap_hml,&umap_hbh,&umap_hbl,gfee1,gfee2,gfee3,gfee4](const std::vector<int>& id,const std::vector<int>& adc,const std::vector<int>& plat,const int& eventid){
        float mh=0.;
        float ml=0.;
        float bh=0.;
        float bl=0.;
        float high=0.;
        float low=0.;
	int nfee1=0;
	int nfee2=0;
	int nfee3=0;
	int nfee4=0;
	std::unordered_map<int,float> umain_cryid_high;
	std::unordered_map<int,float> umain_cryid_low;
	std::unordered_map<int,float> uback_cryid_high;
	std::unordered_map<int,float> uback_cryid_low;
			for(int i=0;i<id.size();i++){
				int cellid = id.at(i);
				int celladc = adc.at(i);
				int cellplat = plat.at(i);
				// hmh/hml/… 与 umap_* 信号：用原始 CellADC；hpmh/… 台基仍为 CellPLAT
				int cryid = cellid/100000;
				int feeid = (cellid%100000)/10000;
				int mbid = (cellid%10000)/1000;
				if(cryid==0){
					std::cout<<cellid<<" "<<cryid<<std::endl;
					continue;
				}
				switch(feeid){
					case 1:gfee1->SetPoint(nfee1,eventid*0.17/60.,cellplat);nfee1++;break;
					case 2:gfee2->SetPoint(nfee2,eventid*0.17/60.,cellplat);nfee1++;break;
					case 3:gfee3->SetPoint(nfee3,eventid*0.17/60.,cellplat);nfee1++;break;
					case 4:gfee4->SetPoint(nfee4,eventid*0.17/60.,cellplat);nfee1++;break;
				}
				if((cellid%1000)/100==1){ // High Gain
					if(mbid==1){
                        mh=(celladc);
                        hpmh->Fill(cellplat);
                        umap_hpmh[cryid-1]->Fill(cellplat);
		    	umap_hmh[cryid-1]->Fill(celladc);
			umain_cryid_high[cryid-1]=mh;
					}
					else{
                        bh=(celladc);
                        hpbh->Fill(cellplat);
                        umap_hpbh[cryid-1]->Fill(cellplat);
		    	umap_hbh[cryid-1]->Fill(celladc);
			uback_cryid_high[cryid-1]=bh;
					}
                    high += (celladc);
				}
				else{ // Low Gain
					if(mbid==1){
                        ml=(celladc);
                        hpml->Fill(cellplat);
                        umap_hpml[cryid-1]->Fill(cellplat);
		    	umap_hml[cryid-1]->Fill(celladc);
			umain_cryid_low[cryid-1]=ml;
					}
					else{
                        bl=(celladc);
                        hpbl->Fill(cellplat);
                        umap_hpbl[cryid-1]->Fill(cellplat);
		    	umap_hbl[cryid-1]->Fill(celladc);
			uback_cryid_low[cryid-1]=bl;
					}
                    low += (celladc);
				}
			}
			for(int i=0;i<25;i++){
				umain_ratio.at(i)->Fill(umain_cryid_low[i],umain_cryid_high[i]);
				uback_ratio.at(i)->Fill(uback_cryid_low[i],uback_cryid_high[i]);
				uhigh_ratio.at(i)->Fill(umain_cryid_high[i],uback_cryid_high[i]);
				ulow_ratio.at(i)->Fill(umain_cryid_low[i],uback_cryid_low[i]);
			}
            hmh->Fill(mh);
            hml->Fill(ml);
            hbh->Fill(bh);
            hbl->Fill(bl);
            hh->Fill(high);
            hl->Fill(low);
			},{"CellID","CellADC","CellPLAT","EventID"});
    std::string oname = fname;
    oname = oname.substr(0,oname.find_last_of("."));
    oname = oname.substr(oname.find_first_of('_')+1);
    //auto fout = new TFile(TString::Format("%s.root",oname.c_str()),"RECREATE");
    auto fout = new TFile("hist.root","RECREATE");
    fout->cd();
    hpmh->Write();
    hpbh->Write();
    hpml->Write();
    hpbl->Write();
    hmh->Write();
    hml->Write();
    hbh->Write();
    hbl->Write();
    hh->Write();
    hl->Write();
    gfee1->Write();
    gfee2->Write();
    gfee3->Write();
    gfee4->Write();
    for(int i=0;i<N;i++){
        umap_hpmh[i]->Write();
        umap_hpml[i]->Write();
        umap_hpbh[i]->Write();
        umap_hpbl[i]->Write();
        umap_hmh[i]->Write();
        umap_hml[i]->Write();
        umap_hbh[i]->Write();
        umap_hbl[i]->Write();
	umain_ratio[i]->Write();
	uback_ratio[i]->Write();
	uhigh_ratio[i]->Write();
	ulow_ratio[i]->Write();
    }
    fout->Close();

}
void draw_temp(const std::string& fname){
	auto fg=[](const TString& name,const TString& title)->TGraph*{
		TGraph *g=new TGraph();
		g->SetName(name);
		g->SetTitle(title);
		return g;
	};
	TGraph *gc[4][3],*gt[4][4];//Graph for current and temperature
	for(int i=0;i<4;i++){
		for(int ic=0;ic<3;ic++){
			gc[i][ic]=fg(TString::Format("f%dc%d",i+1,ic),
				     TString::Format("FEE_%d Current_%d",i+1,ic));
		}
		for(int it=0;it<4;it++){
			gt[i][it]=fg(TString::Format("f%dt%d",i+1,it),
				     TString::Format("FEE_%d Temperature_%d",i+1,it));
		}
	}
	int nc[4][3]={};//FEEID point it
	int nt[4][4]={};
	ROOT::RDataFrame df("caloTree",fname);
	df.Foreach([&](  const int& id,
			const std::vector<float>& c0,
			const std::vector<float>& c1,
			const std::vector<float>& c2,
			const std::vector<float>& t0,
			const std::vector<float>& t1,
			const std::vector<float>& t2,
			const std::vector<float>& t3
			){
		for(int i=0;i<4;i++){ //FEEID
			gc[i][0]->SetPoint(nc[i][0],id,c0[i]);nc[i][0]++;
			gc[i][1]->SetPoint(nc[i][1],id,c1[i]);nc[i][1]++;
			gc[i][2]->SetPoint(nc[i][2],id,c2[i]);nc[i][2]++;
			gt[i][0]->SetPoint(nt[i][0],id,t0[i]);nt[i][0]++;
			gt[i][1]->SetPoint(nt[i][1],id,t1[i]);nt[i][1]++;
			gt[i][2]->SetPoint(nt[i][2],id,t2[i]);nt[i][2]++;
			gt[i][3]->SetPoint(nt[i][3],id,t3[i]);nt[i][3]++;
		}

	},{"TPointID","C0","C1","C2","T0","T1","T2","T3"});
	TFile *fout=new TFile("temp_hist.root","RECREATE");
	fout->cd();
	for(int i=0;i<4;i++){
		gc[i][0]->Write();
		gc[i][1]->Write();
		gc[i][2]->Write();
		gt[i][0]->Write();
		gt[i][1]->Write();
		gt[i][2]->Write();
		gt[i][3]->Write();
	}
	fout->Close();
}
void draw(const float& dataFraction=1.,const std::string& filename=""){
    draw_ped(filename,dataFraction);
    //draw_temp("temp_Calo_cosmic_20250528_1138.root");
}
