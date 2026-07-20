import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")

def _load_graphs(infile, prefix, n=25):
    graphs = []
    for c in range(n):
        g = infile.Get(prefix + str(c))
        graphs.append(g.Clone() if g else None)
    return graphs

def _save_calib_canvas(graphs, png_path, x_range=None):
    colour = ROOT.kCyan + 2
    c = ROOT.TCanvas("calib_" + png_path.replace(".", "_"), "", 0, 0, 1500, 1500)
    c.Divide(5, 5)
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
        if x_range is not None:
            histo.GetXaxis().SetRangeUser(x_range[0], x_range[1])
        histo.SetMinimum(0.0)
        histo.SetMaximum(66000.0)
        c.cd(index + 1)
        histo.Draw("APL")
    c.SaveAs(png_path)

def plot_calib():
    infile = ROOT.TFile.Open("peaks_calo.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError("cannot open peaks_calo.root")
    x_h_mh = (0, 100)
    _save_calib_canvas(_load_graphs(infile, "gh_"), "Calibh.png", x_h_mh)
    _save_calib_canvas(_load_graphs(infile, "gl_"), "Calibl.png")
    _save_calib_canvas(_load_graphs(infile, "gmh_"), "Calibmh.png", x_h_mh)
    _save_calib_canvas(_load_graphs(infile, "gml_"), "Calibml.png")
    infile.Close()

if __name__ == "__main__":
    plot_calib()
