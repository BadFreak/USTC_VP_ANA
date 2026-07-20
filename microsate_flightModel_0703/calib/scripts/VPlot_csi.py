import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")

def plot_calib():
    colour = ROOT.kCyan + 2
    infile = ROOT.TFile.Open("peaks_csi.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError("cannot open peaks_csi.root")

    graphs = []
    for c in range(8):
        g = infile.Get("g_" + str(c))
        graphs.append(g.Clone() if g else None)
    infile.Close()

    c = ROOT.TCanvas("calib_csi", "", 1200, 600)
    c.Divide(4, 2)
    for index, histo in enumerate(graphs):
        if not histo or histo.GetN() < 1:
            continue
        histo.SetMarkerColor(colour)
        histo.SetMarkerSize(0.5)
        histo.SetMarkerStyle(8)
        histo.SetLineColor(colour)
        histo.SetLineWidth(1)
        histo.SetFillColor(0)
        histo.GetXaxis().SetTitle("DAC value")
        histo.GetYaxis().SetTitle("Output ADC")
        histo.GetYaxis().SetTitleOffset(1.5)
        histo.SetMinimum(0.0)
        histo.SetMaximum(66000.0)
        c.cd(index + 1)
        histo.Draw("APL")
    c.SaveAs("Calib.png")


if __name__ == "__main__":
    plot_calib()
