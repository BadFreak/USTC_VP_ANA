#include "Caliber.hh"

Caliber::Caliber() {
    // Constructor implementation
}

Caliber::~Caliber() {
    // Destructor implementation
}

std::vector<float> Caliber::findPeak(TH1D *h) {
    std::vector<float> ret;
    if(h->GetEntries()<1000){
        std::cout<<"Less events for calibration: "<<h->GetEntries()<<std::endl;
        return ret;
    }
	for(int i=0;i<h->GetNbinsX();i++){
		if(h->GetBinContent(i+1)>100.){
			float mean = h->GetBinCenter(i+1);
			TF1 *f1=new TF1("f1","[0]*TMath::Gaus(x,[1],[2])",mean-200,mean+200);
			f1->SetParameters(100,mean,100);
			f1->SetParLimits(2,30,100);
			h->Fit("f1","+q","",mean-200,mean+200);
			float output_mean = f1->GetParameter(1);
			float output_rms = f1->GetParameter(2);
			float output_mag = f1->GetParameter(0);
			// std::cout<<output_mean<<" "<<output_rms<<" "<<output_mag<<std::endl;
			ret.emplace_back(output_mean);
			while(h->GetBinContent(i+1)>50.)i++;
			delete f1;
		}
	}
	return ret;
}

void Caliber::calib(const std::string& det) {
    TGraph *gh[25],*gl[25];
	auto f = TFile::Open(TString::Format("hist_calib_%s.root",det.c_str()),"READ");
	int bhchns[4]={2,3,4,9};
	int blchns[4]={4,9,14,19};
	for(int i=0;i<25;i++){
		gh[i] = new TGraph();
		gh[i]->SetName(Form("gh_%d",i));
		gl[i] = new TGraph();
		gl[i]->SetName(Form("gl_%d",i));
		auto hbh = (TH1D*)f->Get(Form("hbh_%d",i));
		auto hbl = (TH1D*)f->Get(Form("hbl_%d",i));
		auto hpeaks = findPeak(hbh);
		auto lpeaks = findPeak(hbl);
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
	auto fout = new TFile(Form("peaks_%s.root",det.c_str()),"RECREATE");
	fout->cd();
	for(int i=0;i<25;i++){
		gh[i]->Write();
		gl[i]->Write();
	}
	fout->Close();
}
