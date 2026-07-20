#include <set>
#include <cmath>
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"
#include "TH2D.h"
#include "TPad.h"
// nCrystalHit: 只保留“恰好击中该个数晶体”的事例（<0 表示不筛选）
// hitThreshold: 击中判据 ADC - pedestal > hitThreshold
// useLogBin: true 用对数分 bin，false 用线性等宽分 bin
void draw_ped(std::string fname,
              const float& dataFraction,
              const std::string tname,
              int nCrystalHit=-1,
              double hitThreshold=50.,
              bool useLogBin=false){
    TH1I *hitcell = new TH1I("hitcell","hitcell",27,-1,26);
    ROOT::RDataFrame df(tname,fname);
    // 根据 nCrystalHit 决定是否做击中晶体数筛选
    ROOT::RDF::RNode df_work = df;
    if (nCrystalHit >= 0) {
        df_work = df.Filter(
            [nCrystalHit, hitThreshold, &hitcell](const std::vector<int>& CellID,
                                                           const std::vector<int>& CellADC,
                                                           const std::vector<int>& CellPLAT){
                std::set<int> hit_crystals;

                for (size_t i = 0; i < CellID.size(); i++) {
                    if (CellID[i] <= 0) {
                        std::cout << "CellID " << CellID[i] << " <= 0" << std::endl;
                        continue;
                    }
                    int cry_id = CellID[i] / 100000;
                    if (double(CellADC[i]) - double(CellPLAT[i]) > hitThreshold) {
                        hit_crystals.insert(cry_id);
                    }
                }
                int hitcellNo = int(hit_crystals.size());
                hitcell->Fill(hitcellNo);

                return hitcellNo == nCrystalHit;
            },
            {"CellID", "CellADC", "CellPLAT"});
    }
    
    double pmin = 500., pmax = 1800.;  // Ped 直方图横轴范围（全局 + 5x5）
    TH1D *hpmh = new TH1D("hpmh","hplat",100,pmin,pmax);
    TH1D *hpbh = new TH1D("hpbh","hplat",100,pmin,pmax);
    TH1D *hpml = new TH1D("hpml","hplat",100,pmin,pmax);
    TH1D *hpbl = new TH1D("hpbl","hplat",100,pmin,pmax);
    std::unordered_map<int,TH2D*> umain_ratio; // main,cryid, TH2D* ratio
    std::unordered_map<int,TH2D*> uback_ratio; // back,cryid, TH2D* ratio
    std::unordered_map<int,TH2D*> uhigh_ratio; //main and back
    std::unordered_map<int,TH2D*> ulow_ratio; // main and back
    std::unordered_map<int,TGraph*> umain_ratio_g; // 每通道 Main HG/LG 逐事例散点
    std::unordered_map<int,TGraph*> uback_ratio_g; // 每通道 Back HG/LG 逐事例散点
    for(int i=0;i<25;i++){
	    umain_ratio[i] = new TH2D(Form("main_hlr_%d",i),Form("Main_HLRatio_%d",i),100,1,2e4,100,1,2e3);
	    uback_ratio[i] = new TH2D(Form("back_hlr_%d",i),Form("Back_HLRatio_%d",i),100,1,2e4,100,1,2e3);
	    uhigh_ratio[i] = new TH2D(Form("high_mbr_%d",i),Form("High_MBRatio_%d",i),100,1,2e4,100,1,2e4);
	    ulow_ratio[i] = new TH2D(Form("low_mbr_%d",i),Form("Low_MBRatio_%d",i),100,1,2e3,100,1,2e3);
	    umain_ratio_g[i] = new TGraph();
	    umain_ratio_g[i]->SetName(Form("main_hlg_%d",i));
	    umain_ratio_g[i]->SetTitle(Form("Main_HLRatio_%d",i));
	    uback_ratio_g[i] = new TGraph();
	    uback_ratio_g[i]->SetName(Form("back_hlg_%d",i));
	    uback_ratio_g[i]->SetTitle(Form("Back_HLRatio_%d",i));
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
    const int n_bins = 500;
    const double e_min = 400.;    // Sig_hmh / Sig_hbh 等的横轴下限
    const double e_max = 15000.;  // 横轴上限
    static double log_bins[n_bins + 1]; // 对数 bin 边界（如需）

    TH1D *hmh = nullptr;
    TH1D *hml = nullptr;
    TH1D *hbh = nullptr;
    TH1D *hbl = nullptr;
    TH1D *hh  = nullptr;
    TH1D *hl  = nullptr;

    if (useLogBin) {
        double log_min = 0; // 10^0 = 1
        double log_max = 4.4; // 10^4.2 = 15848.9
        double log_bin_width = (log_max - log_min) / n_bins;
        for (int i = 0; i <= n_bins; i++) {
            log_bins[i] = TMath::Power(10, log_min + i * log_bin_width);
        }
        hmh = new TH1D("hmh","hE",n_bins,log_bins);
        hml = new TH1D("hml","hE",n_bins,log_bins);
        hbh = new TH1D("hbh","hE",n_bins,log_bins);
        hbl = new TH1D("hbl","hE",n_bins,log_bins);
        hh  = new TH1D("hh","hE",n_bins,log_bins);
        hl  = new TH1D("hl","hE",n_bins,log_bins);
    } else {
        hmh = new TH1D("hmh","hE",n_bins,e_min,e_max);
        hml = new TH1D("hml","hE",n_bins,e_min,e_max);
        hbh = new TH1D("hbh","hE",n_bins,e_min,e_max);
        hbl = new TH1D("hbl","hE",n_bins,e_min,e_max);
        hh  = new TH1D("hh","hE",n_bins,e_min,e_max);
        hl  = new TH1D("hl","hE",n_bins,e_min,e_max);
    }
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
		umap_hpmh[i]=new TH1D(Form("ped_mh_%d",i),Form("Pedestal_Main_High_%d",i),1000,pmin,pmax);
		umap_hpml[i]=new TH1D(Form("ped_ml_%d",i),Form("Pedestal_Main_Low_%d",i),1000,pmin,pmax);
		umap_hpbh[i]=new TH1D(Form("ped_bh_%d",i),Form("Pedestal_Back_High_%d",i),1000,pmin,pmax);
		umap_hpbl[i]=new TH1D(Form("ped_bl_%d",i),Form("Pedestal_Back_Low_%d",i),1000,pmin,pmax);

        if (useLogBin) {
            umap_hmh[i]=new TH1D(Form("hmh_%d",i),Form("Main_High_%d",i),n_bins,log_bins);
            umap_hml[i]=new TH1D(Form("hml_%d",i),Form("Main_Low_%d",i),n_bins,log_bins);
            umap_hbh[i]=new TH1D(Form("hbh_%d",i),Form("Back_High_%d",i),n_bins,log_bins);
            umap_hbl[i]=new TH1D(Form("hbl_%d",i),Form("Back_Low_%d",i),n_bins,log_bins);
        } else {
            umap_hmh[i]=new TH1D(Form("hmh_%d",i),Form("Main_High_%d",i),n_bins/4,e_min,e_max);
            umap_hml[i]=new TH1D(Form("hml_%d",i),Form("Main_Low_%d",i),n_bins/10,e_min,e_max/10);
            umap_hbh[i]=new TH1D(Form("hbh_%d",i),Form("Back_High_%d",i),n_bins/4,e_min,e_max);
            umap_hbl[i]=new TH1D(Form("hbl_%d",i),Form("Back_Low_%d",i),n_bins/10,e_min,e_max/10);
        }
	}
	float N_total = float(*df_work.Count());
	std::cout << "N_total :" << N_total << std::endl; 
	int N_actual = int(N_total * dataFraction);
	auto dff = df_work.Range(N_actual);
    const bool fillSigHitOnly = (nCrystalHit >= 0);
    dff.Foreach([&uhigh_ratio,&ulow_ratio,&umain_ratio,&uback_ratio,&umain_ratio_g,&uback_ratio_g,hpmh,hpbh,hpml,hpbl,&umap_hpmh,&umap_hpml,&umap_hpbh,&umap_hpbl,&hmh,&hml,&hbh,&hbl,&hh,&hl,&umap_hmh,&umap_hml,&umap_hbh,&umap_hbl,gfee1,gfee2,gfee3,gfee4,hitThreshold,fillSigHitOnly](const std::vector<int>& id,const std::vector<int>& adc,const std::vector<int>& plat,const int& eventid){
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
	std::set<int> hit_crystals;
	std::unordered_map<int,float> umain_cryid_high;
	std::unordered_map<int,float> umain_cryid_low;
	std::unordered_map<int,float> uback_cryid_high;
	std::unordered_map<int,float> uback_cryid_low;
			for(int i=0;i<id.size();i++){
				int cellid = id.at(i);
				int celladc = adc.at(i);
				int cellplat = plat.at(i);
				// 信号 hmh/hml/… 与 umap_*：填 CellADC-CellPLAT；台基直方图仍用 CellPLAT
				int cryid = cellid/100000;
				int feeid = (cellid%100000)/10000;
				int mbid = (cellid%10000)/1000;
				float sig = float(celladc) - float(cellplat);
				if(cryid==0){
					std::cout<<cellid<<" "<<cryid<<std::endl;
					continue;
				}
				if (sig > hitThreshold)
					hit_crystals.insert(cryid - 1);
				switch(feeid){
					case 1:gfee1->SetPoint(nfee1,eventid*0.17/60.,cellplat);nfee1++;break;
					case 2:gfee2->SetPoint(nfee2,eventid*0.17/60.,cellplat);nfee1++;break;
					case 3:gfee3->SetPoint(nfee3,eventid*0.17/60.,cellplat);nfee1++;break;
					case 4:gfee4->SetPoint(nfee4,eventid*0.17/60.,cellplat);nfee1++;break;
				}
				if((cellid%1000)/100==1){ // High Gain
					if(mbid==1){
                        mh=sig;
                        hpmh->Fill(cellplat);
                        umap_hpmh[cryid-1]->Fill(cellplat);
			umain_cryid_high[cryid-1]=mh;
                        if (!fillSigHitOnly)
		    	umap_hmh[cryid-1]->Fill(sig);
					}
					else{
                        bh=sig;
                        hpbh->Fill(cellplat);
                        umap_hpbh[cryid-1]->Fill(cellplat);
			uback_cryid_high[cryid-1]=bh;
                        if (!fillSigHitOnly)
		    	umap_hbh[cryid-1]->Fill(sig);
					}
                    if (!fillSigHitOnly)
                    high += sig;
				}
				else{ // Low Gain
					if(mbid==1){
                        ml=sig;
                        hpml->Fill(cellplat);
                        umap_hpml[cryid-1]->Fill(cellplat);
			umain_cryid_low[cryid-1]=ml;
                        if (!fillSigHitOnly)
		    	umap_hml[cryid-1]->Fill(sig);
					}
					else{
                        bl=sig;
                        hpbl->Fill(cellplat);
                        umap_hpbl[cryid-1]->Fill(cellplat);
			uback_cryid_low[cryid-1]=bl;
                        if (!fillSigHitOnly)
		    	umap_hbl[cryid-1]->Fill(sig);
					}
                    if (!fillSigHitOnly)
                    low += sig;
				}
			}
			if (fillSigHitOnly) {
				mh = ml = bh = bl = 0.;
				high = low = 0.;
				for (int ci : hit_crystals) {
					if (umain_cryid_high.count(ci)) {
						umap_hmh[ci]->Fill(umain_cryid_high[ci]);
						mh = umain_cryid_high[ci];
						high += umain_cryid_high[ci];
					}
					if (umain_cryid_low.count(ci)) {
						umap_hml[ci]->Fill(umain_cryid_low[ci]);
						ml = umain_cryid_low[ci];
						low += umain_cryid_low[ci];
					}
					if (uback_cryid_high.count(ci)) {
						umap_hbh[ci]->Fill(uback_cryid_high[ci]);
						bh = uback_cryid_high[ci];
						high += uback_cryid_high[ci];
					}
					if (uback_cryid_low.count(ci)) {
						umap_hbl[ci]->Fill(uback_cryid_low[ci]);
						bl = uback_cryid_low[ci];
						low += uback_cryid_low[ci];
					}
				}
			}
			for (int i : hit_crystals) {
				if (umain_cryid_high.count(i) && umain_cryid_low.count(i)) {
					const double xh = umain_cryid_high[i];
					const double yl = umain_cryid_low[i];
					umain_ratio.at(i)->Fill(xh, yl);
					auto* gm = umain_ratio_g.at(i);
					gm->SetPoint(gm->GetN(), xh, yl);
				}
				if (uback_cryid_high.count(i) && uback_cryid_low.count(i)) {
					const double xh = uback_cryid_high[i];
					const double yl = uback_cryid_low[i];
					uback_ratio.at(i)->Fill(xh, yl);
					auto* gb = uback_ratio_g.at(i);
					gb->SetPoint(gb->GetN(), xh, yl);
				}
				if (umain_cryid_high.count(i) && uback_cryid_high.count(i))
					uhigh_ratio.at(i)->Fill(umain_cryid_high[i], uback_cryid_high[i]);
				if (umain_cryid_low.count(i) && uback_cryid_low.count(i))
					ulow_ratio.at(i)->Fill(umain_cryid_low[i], uback_cryid_low[i]);
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
    auto fout = new TFile("hist_calo.root","RECREATE");
    fout->cd();
    hitcell->Write();
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
	umain_ratio_g[i]->Write();
	uback_ratio_g[i]->Write();
    }

    // 直接输出 5x5 图：4 张 Ped + 4 张 Sig
    auto draw_grid_5x5 = [&](const char* cname,
                             const char* ctitle,
                             std::unordered_map<int,TH1D*>& hmap,
                             const char* png_name,
                             bool setLogX = false) {
        auto c = new TCanvas(cname, ctitle, 1800, 1200);
        c->Divide(5, 5);
        for (int i = 0; i < N; i++) {
            c->cd(i + 1);
            if (setLogX) gPad->SetLogx();
            auto h = hmap[i];
            if (!h) continue;
            h->SetLineWidth(2);
            h->Draw("hist");
        }
        c->Write();
        c->SaveAs(png_name);
    };

    // draw_grid_5x5("c_ped_mh_5x5", "Pedestal Main High (5x5)", umap_hpmh, "ped_mh_5x5.png", false);
    // draw_grid_5x5("c_ped_ml_5x5", "Pedestal Main Low (5x5)",  umap_hpml, "ped_ml_5x5.png", false);
    // draw_grid_5x5("c_ped_bh_5x5", "Pedestal Back High (5x5)", umap_hpbh, "ped_bh_5x5.png", false);
    // draw_grid_5x5("c_ped_bl_5x5", "Pedestal Back Low (5x5)",  umap_hpbl, "ped_bl_5x5.png", false);

    // draw_grid_5x5("c_sig_mh_5x5", "Signal Main High (5x5)", umap_hmh, "sig_mh_5x5.png", useLogBin);
    // draw_grid_5x5("c_sig_ml_5x5", "Signal Main Low (5x5)",  umap_hml, "sig_ml_5x5.png", useLogBin);
    // draw_grid_5x5("c_sig_bh_5x5", "Signal Back High (5x5)", umap_hbh, "sig_bh_5x5.png", useLogBin);
    // draw_grid_5x5("c_sig_bl_5x5", "Signal Back Low (5x5)",  umap_hbl, "sig_bl_5x5.png", useLogBin);

    auto draw_grid_5x5_h2 = [&](const char* cname, const char* ctitle,
                                std::unordered_map<int, TH2D*>& hmap, const char* png_name) {
        auto c = new TCanvas(cname, ctitle, 1800, 1200);
        c->Divide(5, 5);
        for (int i = 0; i < N; i++) {
            c->cd(i + 1);
            gPad->SetRightMargin(0.12);
            TH2D* h = hmap[i];
            if (!h) continue;
            h->Draw("COLZ");
        }
        c->Write();
        c->SaveAs(png_name);
    };
    // draw_grid_5x5_h2("c_main_hlr_5x5", "Main High/Low ratio (5x5)", umain_ratio, "main_hlr_5x5.png");
    // draw_grid_5x5_h2("c_back_hlr_5x5", "Back High/Low ratio (5x5)", uback_ratio, "back_hlr_5x5.png");

    fout->Close();

}

// nCrystalHit: 只保留击中该个数晶体的事例（<0 不筛选）；hitThreshold: 击中判据 ADC-pedestal > hitThreshold
// useLogBin: true 用对数分 bin，false 用线性等宽分 bin
void draw_calo(const float& dataFraction=1.,
              const std::string& filename="",
              int nCrystalHit=-1,
              double hitThreshold=50.,
              bool useLogBin=false) {
    std::cout << " Draw " << filename;
    if (nCrystalHit >= 0) std::cout << " [filter: nCrystalHit==" << nCrystalHit << ", threshold=" << hitThreshold << "]";
    std::cout << (useLogBin ? " [log bins]" : " [linear bins]") << std::endl;
    TString fname = TString(filename);
    
    draw_ped(filename,dataFraction,"caloTree",nCrystalHit,hitThreshold,useLogBin);

}
