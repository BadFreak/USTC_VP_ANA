#ifndef PLOTTER_HH
#define PLOTTER_HH

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
#include <TTree.h>
#include <TGraph.h>
#include <unordered_map>
#include "ROOT/RDataFrame.hxx"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

class Plotter {
public:
    Plotter() {}
    ~Plotter() {}
    void draw(const float& dataFraction,const std::string& filename);
private:
    void draw_ped(std::string fname,const float& dataFraction,const std::string& det);
    void draw_calib(const std::string& fname,const float& dataFraction,const std::string& det);
    void draw_temp(const std::string& fname);
};

#endif // PLOTTER_HH