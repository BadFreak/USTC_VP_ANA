#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TPad.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TFitResultPtr.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <sys/stat.h>

static void WriteCanvasToPlotAll(TCanvas* c, const TString& saveDir, const TString& keyName)
{
    if (!c)
        return;
    const TString path = saveDir + "/plot_all.root";
    struct stat st;
    const bool exists = (stat(path.Data(), &st) == 0);
    TFile* fout = TFile::Open(path, exists ? "UPDATE" : "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot open " << path << " for write" << std::endl;
        return;
    }
    fout->cd();
    c->Write(keyName, TObject::kOverwrite);
    fout->Close();
    delete fout;
    std::cout << "Wrote " << keyName << " -> " << path << std::endl;
}

// Landau×Gaussian，与 VPlot_calo.py::_langaufun / ROOT langaus.C 一致
static Double_t langaufun(Double_t* x, Double_t* par)
{
    if (par[3] <= 0.)
        return 0.;
    constexpr Double_t invsq2pi = 0.3989422804014;
    constexpr Double_t mpshift = -0.22278298;
    constexpr Double_t np = 100.;
    constexpr Double_t sc = 5.;
    const Double_t mpc = par[1] - mpshift * par[0];
    const Double_t xlow = x[0] - sc * par[3];
    const Double_t xupp = x[0] + sc * par[3];
    const Double_t step = (xupp - xlow) / np;
    Double_t sum = 0.;
    for (int i = 1; i <= int(np / 2); ++i) {
        Double_t xx = xlow + (i - 0.5) * step;
        Double_t fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);
        xx = xupp - (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);
    }
    return par[2] * step * sum * invsq2pi / par[3];
}

static bool FitLandauGauss(TH1D* h, TF1*& ffit, double& mpv)
{
    if (!h || h->Integral() <= 0)
        return false;

    constexpr double fit_xlo = 1000.;
    constexpr double fit_xhi = 14000.;
    constexpr double sv_width = 300.;
    constexpr double sv_sig = 150.;
    constexpr double pllo_mp = 1500.;
    constexpr double plhi_mp = 12000.;
    constexpr double pllo_w = 50.;
    constexpr double plhi_w = 1500.;
    constexpr double pllo_s = 20.;
    constexpr double plhi_s = 500.;

    const int nbins = h->GetNbinsX();
    int blo = std::max(1, std::min(h->GetXaxis()->FindBin(fit_xlo), nbins));
    int bhi = std::max(1, std::min(h->GetXaxis()->FindBin(fit_xhi), nbins));
    if (blo > bhi)
        std::swap(blo, bhi);

    const double integral_fit = h->Integral(blo, bhi);
    if (integral_fit <= 0)
        return false;

    int imax = blo;
    for (int b = blo; b <= bhi; ++b) {
        if (h->GetBinContent(b) > h->GetBinContent(imax))
            imax = b;
    }
    double mp_init = h->GetBinCenter(imax);
    mp_init = std::max(pllo_mp, std::min(plhi_mp, mp_init));

    ffit = new TF1(Form("langau_%s", h->GetName()), langaufun, fit_xlo, fit_xhi, 4);
    ffit->SetParameters(sv_width, mp_init, integral_fit, sv_sig);
    ffit->SetParNames("Width", "MP", "Area", "GSigma");
    ffit->SetParLimits(0, pllo_w, plhi_w);
    ffit->SetParLimits(1, pllo_mp, plhi_mp);
    ffit->SetParLimits(2, integral_fit * 0.01, 1e10);
    ffit->SetParLimits(3, pllo_s, plhi_s);

    const TFitResultPtr res = h->Fit(ffit, "RB0QS");
    mpv = ffit->GetParameter(1);
    bool ok = (mpv >= pllo_mp && mpv <= plhi_mp);
    if (res.Get() && res->Status() != 0 && res->Status() != 1 && res->Status() != 3)
        ok = false;
    return ok;
}


// 支持台基扣除：如果 pedestalVarName 非空，则绘制 varName - pedestalVarName
TH1D* MakeHistFromTree(TTree* tree,
                       const char* varName,
                       int cellID,
                       int nbins,
                       double xmin,
                       double xmax,
                       const char* pedestalVarName = nullptr)
{
    TString hname = Form("h_%s_%d", varName, cellID);
    TString cut   = Form("CellID==%d", cellID);

    TH1D* h = new TH1D(hname, "", nbins, xmin, xmax);
    h->Sumw2();

    TString drawCmd;
    if (pedestalVarName && strlen(pedestalVarName) > 0) {
        drawCmd = Form("%s-%s>>%s", varName, pedestalVarName, hname.Data());
    } else {
        drawCmd = Form("%s>>%s", varName, hname.Data());
    }
    tree->Draw(drawCmd, cut, "goff");

    return h;
}

void BeautifyHist(TH1D* h, const TString& title, const TString& xtitle)
{
    if (!h) return;

    h->SetTitle(title);
    h->SetStats(0);
    h->SetLineWidth(2);
    h->SetFillStyle(0);

    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle("Counts");

    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();

    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);

    h->GetXaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetTitleOffset(1.20);
}

// 新增 pedestalVarName 参数，支持台基扣除
void DrawOneCanvas(TTree* tree,
                   const char* varName,
                   const std::vector<int>& cellIDs,
                   const TString& saveDir,
                   const TString& outName,
                   int nbins,
                   double xmin,
                   double xmax,
                   const char* pedestalVarName = nullptr)
{
    TCanvas* c = new TCanvas(outName, outName, 1800, 900);
    c->Divide(4, 2, 0.001, 0.001);

    for (int i = 0; i < 8; ++i) {
        c->cd(i + 1);

        gPad->SetGrid();
        // gPad->SetLogy();
        // gPad->SetLogx();
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.10);

        int cellID = cellIDs[i];
        TH1D* h = MakeHistFromTree(tree, varName, cellID, nbins, xmin, xmax, pedestalVarName);

        TString xtitle = varName;
        if (pedestalVarName && strlen(pedestalVarName) > 0) {
            xtitle = Form("%s - %s", varName, pedestalVarName);
        }
        TString title = Form("%s, CellID = %d", xtitle.Data(), cellID);
        BeautifyHist(h, title, xtitle);

        h->Draw("hist");

        if (pedestalVarName && strlen(pedestalVarName) > 0 && TString(varName) == "CellADC") {
            TF1* ffit = nullptr;
            double mpv = 0.;
            if (FitLandauGauss(h, ffit, mpv)) {
                ffit->SetLineColor(kRed + 1);
                ffit->SetLineWidth(2);
                ffit->Draw("same");

                auto* pt = new TPaveText(0.68, 0.82, 0.97, 0.96, "NDC");
                pt->SetBorderSize(0);
                pt->SetFillStyle(0);
                pt->SetTextFont(42);
                pt->SetTextSize(0.04);
                pt->AddText(Form("MPV: %d", int(std::lround(mpv))));
                pt->Draw();
            }
        } else if (TString(varName) == "CellPLAT") {
            // 仅对 CellPLAT 做高斯拟合，并在图例显示 mean/sigma
            TF1* fgaus = new TF1(Form("fgaus_%s_%d", varName, cellID), "gaus", 700, 1300);
            h->Fit(fgaus, "RQ0");
            fgaus->SetLineColor(kRed + 1);
            fgaus->SetLineWidth(2);
            fgaus->Draw("same");

            TLegend* leg = new TLegend(0.58, 0.72, 0.94, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.035);
            leg->AddEntry(fgaus, "gaus fit [700,1300]", "l");
            leg->AddEntry((TObject*)nullptr, Form("mean = %.2f", fgaus->GetParameter(1)), "");
            leg->AddEntry((TObject*)nullptr, Form("sigma = %.2f", fgaus->GetParameter(2)), "");
            leg->Draw();
        }
    }

    c->Update();

    TString pngName = saveDir + "/" + outName + ".png";
    c->SaveAs(pngName);
    WriteCanvasToPlotAll(c, saveDir, outName);

    std::cout << "Saved: " << pngName << std::endl;
}

void draw_csi(const char* filename = "result_Calo_cosmic_20260403_1422.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFontSize(0.045);

    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree* csiTree = dynamic_cast<TTree*>(f->Get("csiTree"));
    if (!csiTree) {
        std::cerr << "Error: cannot find tree csiTree in " << filename << std::endl;
        f->Close();
        return;
    }

    std::vector<int> cellIDs = {
        100000, 200000, 300000, 400000,
        500000, 600000, 700000, 800000
    };

    TString fullpath(filename);
    Ssiz_t pos = fullpath.Last('/');
    TString saveDir = ".";
    if (pos != kNPOS) saveDir = fullpath(0, pos);


    // CellADC 扣除台基（CellPLAT）：固定范围 1000 ~ 14000
    DrawOneCanvas(csiTree, "CellADC", cellIDs, saveDir,
                  "CellADC_4x2", 200, 1000, 14000, "CellPLAT");

    // CellPLAT: 固定范围 700 ~ 1300
    DrawOneCanvas(csiTree, "CellPLAT", cellIDs, saveDir,
                  "CellPLAT_4x2", 100, 700, 1300);

    f->Close();
}
