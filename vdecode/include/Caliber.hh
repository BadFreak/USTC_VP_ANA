#ifndef CALIBER_HH
#define CALIBER_HH
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TGraph.h>
#include <unordered_map>
#include "ROOT/RDataFrame.hxx"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

class Caliber{
    public:
        Caliber();
        ~Caliber();
        void calib(const std::string& det);
    private:
        std::vector<float> findPeak(TH1D *h);
        
};
#endif