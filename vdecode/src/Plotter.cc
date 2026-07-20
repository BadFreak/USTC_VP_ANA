#include "Plotter.hh"
void Plotter::draw_ped(std::string fname,const float& dataFraction,const std::string& det){
    int pmin=0.,pmax=50000;
    ROOT::RDataFrame df(det+std::string("Tree"),fname);

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
    const int n_bins = 200;
    double bins[n_bins + 1]; // 存储 bin 边界
    double log_min = 0; // 最小值的对数
    double log_max = 7;
    double log_bin_width = (log_max - log_min) / n_bins; // 对数 bin 宽度
    for (int i = 0; i <= n_bins; i++) {
        bins[i] = TMath::Power(10, log_min + i * log_bin_width);
    }
    TH1D *hmh = new TH1D("hmh","hE",n_bins,bins);
    TH1D *hml = new TH1D("hml","hE",n_bins,bins);
    TH1D *hbh = new TH1D("hbh","hE",n_bins,bins);
    TH1D *hbl = new TH1D("hbl","hE",n_bins,bins);
    TH1D *hh = new TH1D("hh","hE",n_bins,bins);
    TH1D *hl = new TH1D("hl","hE",n_bins,bins);
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
		umap_hmh[i]=new TH1D(Form("hmh_%d",i),Form("Main_High_%d",i),n_bins,bins);
		umap_hml[i]=new TH1D(Form("hml_%d",i),Form("Main_Low_%d",i),n_bins,bins);
		umap_hbh[i]=new TH1D(Form("hbh_%d",i),Form("Back_High_%d",i),n_bins,bins);
		umap_hbl[i]=new TH1D(Form("hbl_%d",i),Form("Back_Low_%d",i),n_bins,bins);
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
                        mh=(celladc-cellplat);
                        hpmh->Fill(cellplat);
                        umap_hpmh[cryid-1]->Fill(cellplat);
		    	umap_hmh[cryid-1]->Fill(celladc-cellplat);
			umain_cryid_high[cryid-1]=mh;
					}
					else{
                        bh=(celladc-cellplat);
                        hpbh->Fill(cellplat);
                        umap_hpbh[cryid-1]->Fill(cellplat);
		    	umap_hbh[cryid-1]->Fill(celladc-cellplat);
			uback_cryid_high[cryid-1]=bh;
					}
                    high += (celladc-cellplat);
				}
				else{ // Low Gain
					if(mbid==1){
                        ml=(celladc-cellplat);
                        hpml->Fill(cellplat);
                        umap_hpml[cryid-1]->Fill(cellplat);
		    	umap_hml[cryid-1]->Fill(celladc-cellplat);
			umain_cryid_low[cryid-1]=ml;
					}
					else{
                        bl=(celladc-cellplat);
                        hpbl->Fill(cellplat);
                        umap_hpbl[cryid-1]->Fill(cellplat);
		    	umap_hbl[cryid-1]->Fill(celladc-cellplat);
			uback_cryid_low[cryid-1]=bl;
					}
                    low += (celladc-cellplat);
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
    TString oname = TString("hist_");
    oname+=TString(det.c_str());
    oname+=TString(".root");
    auto fout = new TFile(oname,"RECREATE");
    // auto fout = new TFile("hist.root","RECREATE");
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
    delete hpmh;delete hpbh;delete hpml;delete hpbl;delete hmh;delete hml;delete hbh;delete hbl;delete hh;delete hl;delete gfee1;delete gfee2;delete gfee3;delete gfee4;
    for(int i=0;i<N;i++){
        delete umap_hpmh[i];delete umap_hpml[i];delete umap_hpbh[i];delete umap_hpbl[i];delete umap_hmh[i];delete umap_hml[i];delete umap_hbh[i];delete umap_hbl[i];
        delete umain_ratio[i];delete uback_ratio[i];delete uhigh_ratio[i];delete ulow_ratio[i];
    }

}
void Plotter::draw_calib(const std::string& fname,const float& dataFraction,const std::string& det){
    int pmin=1e0,pmax=1e5;
    ROOT::RDataFrame dfo(det+std::string("Tree"),fname);
    auto df = dfo.Filter("PackageID==8");
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
                        mh=(celladc-cellplat);
                        hpmh->Fill(cellplat);
                        umap_hpmh[cryid-1]->Fill(cellplat);
		    	umap_hmh[cryid-1]->Fill(celladc-cellplat);
			umain_cryid_high[cryid-1]=mh;
					}
					else{
                        bh=(celladc-cellplat);
                        hpbh->Fill(cellplat);
                        umap_hpbh[cryid-1]->Fill(cellplat);
		    	umap_hbh[cryid-1]->Fill(celladc-cellplat);
			uback_cryid_high[cryid-1]=bh;
					}
                    high += (celladc-cellplat);
				}
				else{ // Low Gain
					if(mbid==1){
                        ml=(celladc-cellplat);
                        hpml->Fill(cellplat);
                        umap_hpml[cryid-1]->Fill(cellplat);
		    	umap_hml[cryid-1]->Fill(celladc-cellplat);
			umain_cryid_low[cryid-1]=ml;
					}
					else{
                        bl=(celladc-cellplat);
                        hpbl->Fill(cellplat);
                        umap_hpbl[cryid-1]->Fill(cellplat);
		    	umap_hbl[cryid-1]->Fill(celladc-cellplat);
			uback_cryid_low[cryid-1]=bl;
					}
                    low += (celladc-cellplat);
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
    auto fout = new TFile(TString::Format("hist_calib_%s.root",det.c_str()),"RECREATE");
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
    delete hpmh;delete hpbh;delete hpml;delete hpbl;delete hmh;delete hml;delete hbh;delete hbl;delete hh;delete hl;delete gfee1;delete gfee2;delete gfee3;delete gfee4;
    for(int i=0;i<N;i++){
        delete umap_hpmh[i];delete umap_hpml[i];delete umap_hpbh[i];delete umap_hpbl[i];delete umap_hmh[i];delete umap_hml[i];delete umap_hbh[i];delete umap_hbl[i];
        delete umain_ratio[i];delete uback_ratio[i];delete uhigh_ratio[i];delete ulow_ratio[i];
    }
}
void Plotter::draw_temp(const std::string& fname){
	auto fg=[](const TString& name,const TString& title)->TGraph*{
		TGraph *g=new TGraph();
		g->SetName(name);
		g->SetTitle(title);
		return g;
	};
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
	ROOT::RDataFrame df("hkTree",fname);
	df.Foreach([&](const int& telFee,
			const float& c0, const float& c1, const float& c2,
			const float& t0, const float& t1, const float& t2, const float& t3){
		double x = static_cast<double>(seq++);
		int s = slotTel(telFee);
		if (s < 0) return;
		gc[s][0]->SetPoint(nc[s][0], x, c0); nc[s][0]++;
		gc[s][1]->SetPoint(nc[s][1], x, c1); nc[s][1]++;
		gc[s][2]->SetPoint(nc[s][2], x, c2); nc[s][2]++;
		gt[s][0]->SetPoint(nt[s][0], x, t0); nt[s][0]++;
		gt[s][1]->SetPoint(nt[s][1], x, t1); nt[s][1]++;
		gt[s][2]->SetPoint(nt[s][2], x, t2); nt[s][2]++;
		gt[s][3]->SetPoint(nt[s][3], x, t3); nt[s][3]++;
	},{"TelFEE","C0","C1","C2","T0","T1","T2","T3"});
	TFile *fout=new TFile("hist_temp.root","RECREATE");
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
	fout->Close();
}
void Plotter::draw(const float& dataFraction=1.,const std::string& filename=""){
    draw_ped(filename,dataFraction,std::string("calo"));
    draw_ped(filename,dataFraction,std::string("csi"));
    draw_calib(filename,dataFraction,std::string("calo"));
    draw_calib(filename,dataFraction,std::string("csi"));
    draw_temp(filename);
}
