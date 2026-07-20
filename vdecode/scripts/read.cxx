#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
int ihist=0;
constexpr int NChn=26;
constexpr int NSam=128;
std::vector<std::vector<TH1D*>> vvth1d;
TFile *fout=new TFile("result.root","RECREATE");
TTree *tout=new TTree("tree","tree");
float chn3,chn8,c3p,c8p;
bool find_head(const std::ifstream *file){
	char buffer[1];
	while(file->rdbuf()->sgetn(buffer,sizeof(buffer))>0){
		unsigned char byte = static_cast<unsigned char>(buffer[0]);
		if(byte == 0x55){
			char buffer2[3];
			file->rdbuf()->sgetn(buffer2,sizeof(buffer2));
			unsigned char b1 = static_cast<unsigned char>(buffer2[0]);
			unsigned char b2 = static_cast<unsigned char>(buffer2[1]);
			unsigned char b3 = static_cast<unsigned char>(buffer2[2]);
			if(b1 == 0xaa && b2 == 0xeb && b3 == 0x90){
				return true;
			}
		}
	}
	return false;
}

bool read_rest(const std::ifstream *file){
	char btime[10];
	file->rdbuf()->sgetn(btime,sizeof(btime));
	char rest[6662];
	file->rdbuf()->sgetn(rest,sizeof(rest));
	unsigned char b1 = static_cast<unsigned char>(rest[6660]);
	unsigned char b2 = static_cast<unsigned char>(rest[6661]);
	if(b1==0x5a && b2==0xa5){
		//Create hist
		std::vector<TH1D*> vth1d;
		chn3=0.;chn8=0.;
		for(size_t chn_i=1;chn_i<=26;chn_i++){
			if(chn_i!=3 && chn_i!=8)continue;
			TString name = TString::Format("chn%d",chn_i);
			//vth1d.emplace_back(new TH1D(name,name,128,0,128));
			for(size_t byte_i=0;byte_i<128;byte_i++){
				int init_pos = 256*(chn_i-1) + byte_i*2;
				auto d1 = static_cast<unsigned char>(rest[init_pos]);
				auto d2 = static_cast<unsigned char>(rest[init_pos+1]);
				unsigned short fullvalue = (static_cast<unsigned short>(d1)<<8) | d2;
				if(chn_i==3)chn3 = fullvalue > chn3 ? fullvalue:chn3;
				if(chn_i==8)chn8 = fullvalue > chn8 ? fullvalue:chn8;
				if(chn_i==3 && byte_i==0)c3p=fullvalue;
				if(chn_i==8 && byte_i==0)c8p=fullvalue;
				//if(init_pos==0)std::cout<<std::hex<<static_cast<int>(d1)<<" "<<static_cast<int>(d2)<<" "<<fullvalue<<std::dec<<" "<<fullvalue<<std::endl;
				//vth1d[vth1d.size()-1]->SetBinContent(byte_i+1,fullvalue);
			}
		}
		chn3-=c3p;
		chn8-=c8p;
		tout->Fill();
		//vvth1d.emplace_back(vth1d);
		ihist++;
		return true;
	}
	return false;
}

int main(int argc,char** argv) {
	for(int i=0;i<NChn;i++){
		TString bname =TString("chn"+std::to_string(i));
	}
tout->Branch("chn3",&chn3);
tout->Branch("chn8",&chn8);
tout->Branch("c3p",&c3p);
tout->Branch("c8p",&c8p);
    auto file = new std::ifstream(std::string(argv[1]), std::ios::binary);
    if (!file->is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }


    while(!file->eof()){
	    if(find_head(file)){
		    if(read_rest(file)){
		    }
		    else{
			    continue;
		    }
	    }
	    else{break;}
    }
    file->close();
    delete file;
    fout->cd();
    tout->Write();
    //for(int i=0;i<vvth1d.size();i++){
    //        fout->mkdir(TString::Format("event%d",i));
    //        fout->cd(TString::Format("event%d",i));
    //        for(auto h:vvth1d.at(i)){
    //    	    h->Write();
    //        }
    //}
    fout->Close();
    return 0;
}
