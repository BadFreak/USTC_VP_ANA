#!/usr/bin/env python3
"""vdecode 结束后：caloTree / csiTree 各画 PackageCount、TriggerID vs Entry(AP)、TimeCode(hist)。"""
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# 全量 Draw("...:Entry$", "", "AP") 在百万级条目下 PyROOT 会 segfault，需子采样
MAX_GRAPH_PTS = 50000


def _draw_vs_entry(tree, ybranch, title, outpath):
    nent = tree.GetEntries()
    step = 1 if nent <= MAX_GRAPH_PTS else (nent + MAX_GRAPH_PTS - 1) // MAX_GRAPH_PTS
    cut = "" if step == 1 else f"Entry$%{step}==0"
    c = ROOT.TCanvas(os.path.basename(outpath).replace(".", "_"), "", 900, 600)
    drawn = tree.Draw(f"{ybranch}:Entry$", cut, "AP")
    if drawn <= 0:
        c.Close()
        raise RuntimeError(f"Draw failed: {ybranch}:Entry$ cut={cut!r}")
    prim = c.GetListOfPrimitives()
    if prim and prim.GetSize() > 0 and prim.At(0).InheritsFrom("TGraph"):
        prim.At(0).SetTitle(title)
        prim.At(0).SetMarkerStyle(20)
        prim.At(0).SetMarkerSize(0.5)
    c.SaveAs(outpath)
    c.Close()
    note = f"(entries={nent}, plotted={drawn}"
    if step > 1:
        note += f", step={step}"
    note += ")"
    print("saved", outpath, note)


def _draw_timecode_hist(tree, title, outpath):
    nent = tree.GetEntries()
    c = ROOT.TCanvas(os.path.basename(outpath).replace(".", "_"), "", 900, 600)
    tree.Draw("TimeCode", "", "hist")
    if tree.GetHistogram():
        tree.GetHistogram().SetTitle(title)
    c.SaveAs(outpath)
    c.Close()
    print("saved", outpath, f"(entries={nent})")


def plot_tree_after_decode(tree, prefix, outdir):
    """prefix 如 calo / csi，输出 caloPackageCount_vs_Entry.png 等。"""
    nent = tree.GetEntries()
    if nent < 1:
        print(f"{tree.GetName()} empty, skip")
        return

    _draw_vs_entry(
        tree,
        "PackageCount",
        f"{prefix} PackageCount vs Entry;Entry;PackageCount",
        os.path.join(outdir, f"{prefix}PackageCount_vs_Entry.png"),
    )
    _draw_vs_entry(
        tree,
        "TriggerID",
        f"{prefix} TriggerID vs Entry;Entry;TriggerID",
        os.path.join(outdir, f"{prefix}TriggerID_vs_Entry.png"),
    )
    _draw_timecode_hist(
        tree,
        f"{prefix} TimeCode",
        os.path.join(outdir, f"{prefix}TimeCode.png"),
    )


def plot_after_decode(rootfile, outdir):
    fin = ROOT.TFile.Open(rootfile, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError(f"cannot open {rootfile}")

    os.makedirs(outdir, exist_ok=True)

    for tname, prefix in [("caloTree", "calo"), ("csiTree", "csi")]:
        tree = fin.Get(tname)
        if not tree:
            print(f"no {tname} in {rootfile}, skip")
            continue
        print(f"plot {tname} -> {prefix}+*.png")
        plot_tree_after_decode(tree, prefix, outdir)

    fin.Close()


if __name__ == "__main__":
    plot_after_decode(
        sys.argv[1] if len(sys.argv) > 1 else "result.root",
        sys.argv[2] if len(sys.argv) > 2 else ".",
    )
