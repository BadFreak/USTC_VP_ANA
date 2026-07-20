#include <set>

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
    int LoopNu = 0;
    ROOT::RDataFrame df(tname,fname);
    // 根据 nCrystalHit 决定是否做击中晶体数筛选
    ROOT::RDF::RNode df_work = df;
    if (nCrystalHit >= 0) {
        df_work = df.Filter(
            [nCrystalHit, hitThreshold, &hitcell, &LoopNu](const std::vector<int>& CellID,
                                                           const std::vector<int>& CellADC,
                                                           const std::vector<int>& CellPLAT){
                std::set<int> hit_crystals;
                for (size_t i = 0; i < CellID.size(); i++) {
                    if (CellID[i] <= 0) continue;
                    int cry_id = CellID[i] / 100000;
                    if (double(CellADC[i]) - double(CellPLAT[i]) > hitThreshold)
                        hit_crystals.insert(cry_id);
                }
                int hitcellNo = int(hit_crystals.size());
                // std::cout << "hitcellNo : " << hitcellNo << std::endl;
                hitcell->Fill(hitcellNo);
                LoopNu++;
                // std::cout << "LoopNu : " << LoopNu
                //           << " histEntry: " << hitcell->GetEntries() << std::endl;

                // 这里用“恰好 == nCrystalHit”，如需 <= 可改为 <=
                return hitcellNo <= nCrystalHit;
            },
            {"CellID", "CellADC", "CellPLAT"});
    }
    
    int pmin=200.,pmax=2000;
    TH1D *hpmh = new TH1D("hpmh","hplat",100,pmin,pmax);
    TH1D *hpbh = new TH1D("hpbh","hplat",100,pmin,pmax);
    TH1D *hpml = new TH1D("hpml","hplat",100,pmin,pmax);
    TH1D *hpbl = new TH1D("hpbl","hplat",100,pmin,pmax);
    std::unordered_map<int,TH2D*> umain_ratio; // main,cryid, TH2D* ratio
    std::unordered_map<int,TH2D*> uback_ratio; // back,cryid, TH2D* ratio
    std::unordered_map<int,TH2D*> uhigh_ratio; //main and back
    std::unordered_map<int,TH2D*> ulow_ratio; // main and back
    for(int i=0;i<25;i++){
	    umain_ratio[i] = new TH2D(Form("main_hlr_%d",i),Form("Main_HLRatio_%d",i),100,1,2e4,100,1,2e3);
	    uback_ratio[i] = new TH2D(Form("back_hlr_%d",i),Form("Back_HLRatio_%d",i),100,1,2e4,100,1,2e3);
	    uhigh_ratio[i] = new TH2D(Form("high_mbr_%d",i),Form("High_MBRatio_%d",i),100,1,2e4,100,1,2e4);
	    ulow_ratio[i] = new TH2D(Form("low_mbr_%d",i),Form("Low_MBRatio_%d",i),100,1,2e3,100,1,2e3);
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
        double log_max = 7; // 10^7
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
            umap_hmh[i]=new TH1D(Form("hmh_%d",i),Form("Main_High_%d",i),n_bins,e_min,e_max);
            umap_hml[i]=new TH1D(Form("hml_%d",i),Form("Main_Low_%d",i),n_bins,e_min,e_max);
            umap_hbh[i]=new TH1D(Form("hbh_%d",i),Form("Back_High_%d",i),n_bins,e_min,e_max);
            umap_hbl[i]=new TH1D(Form("hbl_%d",i),Form("Back_Low_%d",i),n_bins,e_min,e_max);
        }
	}
	float N_total = float(*df_work.Count());
	std::cout << "N_total :" << N_total << std::endl; 
	int N_actual = int(N_total * dataFraction);
	auto dff = df_work.Range(N_actual);
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
				umain_ratio.at(i)->Fill(umain_cryid_high[i],umain_cryid_low[i]);
				uback_ratio.at(i)->Fill(uback_cryid_high[i],uback_cryid_low[i]);
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
    }
    fout->Close();

}
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
	TGraph *gc[4][3],*gt[4][4];//Graph for current and temperature (MATLAB: FEE 1-4, Current 0-2, Temperature 0-3)
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
	int nc[4][3]={};
	int nt[4][4]={};
	ROOT::RDataFrame df("hkTree", fname);
	df.Foreach([&](  const int& id,
			const std::vector<float>& c0,
			const std::vector<float>& c1,
			const std::vector<float>& c2,
			const std::vector<float>& t0,
			const std::vector<float>& t1,
			const std::vector<float>& t2,
			const std::vector<float>& t3
			){
		for(int i=0;i<4;i++){
			gc[i][0]->SetPoint(nc[i][0],id,c0[i]);nc[i][0]++;
			gc[i][1]->SetPoint(nc[i][1],id,c1[i]);nc[i][1]++;
			gc[i][2]->SetPoint(nc[i][2],id,c2[i]);nc[i][2]++;
			gt[i][0]->SetPoint(nt[i][0],id,t0[i]);nt[i][0]++;
			gt[i][1]->SetPoint(nt[i][1],id,t1[i]);nt[i][1]++;
			gt[i][2]->SetPoint(nt[i][2],id,t2[i]);nt[i][2]++;
			gt[i][3]->SetPoint(nt[i][3],id,t3[i]);nt[i][3]++;
		}
	},{"TPoint","C0","C1","C2","T0","T1","T2","T3"});
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
// nCrystalHit: 只保留击中该个数晶体的事例（<0 不筛选）；hitThreshold: 击中判据 ADC-pedestal > hitThreshold
// useLogBin: true 用对数分 bin，false 用线性等宽分 bin
void draw_new(const float& dataFraction=1.,
              const std::string& filename="",
              int nCrystalHit=-1,
              double hitThreshold=50.,
              bool useLogBin=false) {
    std::cout << " Draw " << filename;
    if (nCrystalHit >= 0) std::cout << " [filter: nCrystalHit==" << nCrystalHit << ", threshold=" << hitThreshold << "]";
    std::cout << (useLogBin ? " [log bins]" : " [linear bins]") << std::endl;
    TString fname = TString(filename);
  
    draw_ped(filename,dataFraction,"caloTree",nCrystalHit,hitThreshold,useLogBin);
    
    draw_temp(filename);
}
