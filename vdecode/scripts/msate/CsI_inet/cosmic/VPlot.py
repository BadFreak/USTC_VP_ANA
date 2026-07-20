# An example python script for making a tidy plot with several lines on it
# I'm using TProfiles because that's what is in the TIDA output, but 
# the standard histogram is a "TH1F" and the exact same commands work for it.
# If you want to use graphs instead of histograms you can check out
# the other example script in here.

import ROOT
import math

# Imports ATLAS style for plotting
# You have to have set it up first (see README for instructions)
# You can run it without this but it will have an ugly stats box and so on
# that you'd have to turn off manually.
#import VPStyle
#ROOT.SetVPStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
def plot_ct(cort):
  # Load some histos from the example file
  # (Gaussian limits from 2016 TLA conf)
  infile = ROOT.TFile.Open("temp_hist.root","READ")
  histos = []
  legendLines = []

  mini=9999.
  maxi=0.
  cort_index=3 if cort=="c" else 4
  for i in range(0,4):
    for c in range(0,cort_index): #Current plots
      hname = "f"+str(i+1)+str(cort)+str(c)
      histo = infile.Get(hname)
      histos.append(histo)
  # Close the input file
  infile.Close()

  # Make a canvas to put the plot on.
  # We don't want log axes for this plot, 
  # but if you do you can control them here.
  c = ROOT.TCanvas("canvas",'',0,0,1500,1500)
  c.Divide(4,cort_index)
  c.SetLogx(False)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what x and y range to use in the display.
  xRange = [1,1e5]
  yRange = [0.9,1.1]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kOrange,ROOT.kBlack]

  # Make a legend.
  # These are the locations of the left side, bottom side, right
  # side, and top, as fractions of the canvas.
  

  # Draw each histogram.
  # You really shouldn't put two histograms with different
  # x axes on the same plot - I'm only doing it here
  # to show you how to draw multiple plots on the same
  # canvas.
  for histo in histos :

    index = histos.index(histo)
    colour = goodColours[0]

    # Set up marker to look nice
    histo.SetMarkerColor(colour)
    histo.SetMarkerSize(0)  
    histo.SetMarkerStyle(0)

    # Set up line to look nice
    histo.SetLineColor(colour)
    histo.SetLineWidth(3)
    histo.SetLineStyle(1)

    # Make sure we don't get a fill
    histo.SetFillColor(0)

    # Label my axes!!
    histo.GetXaxis().SetTitle("Time[s]")
    histo.GetYaxis().SetTitle("Temperature[Celsius]")
    # Move the label around if you want
    histo.GetYaxis().SetTitleOffset(1.5)
    maxi=-9999.
    mini=9999.
    for i in range(histo.GetN()):
        maxi = histo.GetPointY(i) if histo.GetPointY(i) > maxi else maxi
        mini = histo.GetPointY(i) if histo.GetPointY(i) < mini else mini
    hrange = maxi-mini
    mini = (mini - hrange) if cort == "t" else (mini - 2.*hrange)
    maxi = (maxi + hrange) if cort == "t" else (maxi + 2.*hrange)


    # Set the limit for the axes
    #histo.GetXaxis().SetRangeUser(xRange[0],xRange[1])
    histo.GetYaxis().SetRangeUser(mini,maxi)
    c.cd(index+1)
    histo.Draw("") # Draw data points (you'll get error bars by default)

    # Fill entry into legend
    # "PL" means both the line and point style
    # will show up in the legend.
  oname = "Current.png" if cort=="c" else "Temperature.png"
  c.SaveAs(oname)
def plot_sig(name):
  # Load some histos from the example file
  # (Gaussian limits from 2016 TLA conf)
  infile = ROOT.TFile.Open("hist.root","READ")
  histos = []
  legendLines = []
  hname = []

  mini=9999.
  maxi=0.
  for i in range(0,25):
    histo = infile.Get(name+"_"+str(i))
    histo.SetDirectory(0)
    if(histo.Integral()>0.):
        histo.Scale(1./histo.Integral())
    histos.append(histo)
  # Close the input file
  infile.Close()

  # Make a canvas to put the plot on.
  # We don't want log axes for this plot, 
  # but if you do you can control them here.
  c = ROOT.TCanvas("canvas",'',0,0,1500,1500)
  c.Divide(5,5)
  c.SetLogx(True)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what x and y range to use in the display.
  xRange = [1,1e5]
  yRange = [0.9,1.1]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kOrange,ROOT.kBlack]

  # Make a legend.
  # These are the locations of the left side, bottom side, right
  # side, and top, as fractions of the canvas.
  

  # Draw each histogram.
  # You really shouldn't put two histograms with different
  # x axes on the same plot - I'm only doing it here
  # to show you how to draw multiple plots on the same
  # canvas.
  for histo in histos :

    index = histos.index(histo)
    colour = goodColours[0]

    # Set up marker to look nice
    histo.SetMarkerColor(colour)
    histo.SetMarkerSize(0)  
    histo.SetMarkerStyle(0)

    # Set up line to look nice
    histo.SetLineColor(colour)
    histo.SetLineWidth(3)
    histo.SetLineStyle(1)

    # Make sure we don't get a fill
    histo.SetFillColor(0)

    # Label my axes!!
    histo.GetXaxis().SetTitle("ADC")
    histo.GetYaxis().SetTitle("Fraction")
    # Move the label around if you want
    histo.GetYaxis().SetTitleOffset(1.5)

    # Set the limit for the axes
    histo.GetXaxis().SetRangeUser(xRange[0],xRange[1])
    histo.GetYaxis().SetRangeUser(0.,histo.GetMaximum())
    c.cd(index+1)
    ROOT.gPad.SetLogx(True)
    ROOT.gPad.SetLogy(False)

    if index==0 :
      histo.Draw("H") # Draw data points (you'll get error bars by default)
    else :
      histo.Draw("H") # SAME means don't get rid of the previous stuff on the canvas
    legend = ROOT.TLegend(0.6,0.72,0.92,0.92)
    # Make the text a nice fond, and big enough
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    # A few more formatting things .....
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(3)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.AddEntry(histo,str(index))
    legend.Draw()
    # Fill entry into legend
    # "PL" means both the line and point style
    # will show up in the legend.
  c.SaveAs("Sig_"+name+".png")

def plot_sub(name):
  # Load some histos from the example file
  # (Gaussian limits from 2016 TLA conf)
  infile = ROOT.TFile.Open("hist.root","READ")
  histos = []
  legendLines = []
  hname = []

  mini=9999.
  maxi=0.
  for i in range(0,25):
    histo = infile.Get(name+"_"+str(i))
    histo.SetDirectory(0)
    if(histo.Integral()>0.):
        histo.Scale(1./histo.Integral())
    histos.append(histo)
  # Close the input file
  infile.Close()

  # Make a canvas to put the plot on.
  # We don't want log axes for this plot, 
  # but if you do you can control them here.
  c = ROOT.TCanvas("canvas",'',0,0,1500,1500)
  c.Divide(5,5)
  c.SetLogx(False)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what x and y range to use in the display.
  xRange = [-3.2,3.2]
  yRange = [0.9,1.1]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kOrange,ROOT.kBlack]

  # Make a legend.
  # These are the locations of the left side, bottom side, right
  # side, and top, as fractions of the canvas.
  

  # Draw each histogram.
  # You really shouldn't put two histograms with different
  # x axes on the same plot - I'm only doing it here
  # to show you how to draw multiple plots on the same
  # canvas.
  for histo in histos :

    index = histos.index(histo)
    colour = goodColours[0]

    # Set up marker to look nice
    histo.SetMarkerColor(colour)
    histo.SetMarkerSize(0)  
    histo.SetMarkerStyle(0)

    # Set up line to look nice
    histo.SetLineColor(colour)
    histo.SetLineWidth(3)
    histo.SetLineStyle(1)

    # Make sure we don't get a fill
    histo.SetFillColor(0)

    # Label my axes!!
    histo.GetXaxis().SetTitle("ADC")
    histo.GetYaxis().SetTitle("Fraction")
    # Move the label around if you want
    histo.GetYaxis().SetTitleOffset(1.5)

    # Set the limit for the axes
    #histo.GetXaxis().SetLimits(xRange[0],xRange[1])
    histo.GetYaxis().SetRangeUser(0.,0.5)
    c.cd(index+1)

    if index==0 :
      histo.Draw("H") # Draw data points (you'll get error bars by default)
    else :
      histo.Draw("H") # SAME means don't get rid of the previous stuff on the canvas
    legend = ROOT.TLegend(0.6,0.72,0.92,0.92)
    # Make the text a nice fond, and big enough
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    # A few more formatting things .....
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(3)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.AddEntry(histo,str(index))
    legend.Draw()
    # Fill entry into legend
    # "PL" means both the line and point style
    # will show up in the legend.
  c.SaveAs("Plat_"+name+".png")
def plot_2dsingle(name,leg,txt,r1,r2):
  # Load some histos from the example file
  # (Gaussian limits from 2016 TLA conf)
  infile = ROOT.TFile.Open("hist.root","READ")
  histos = []
  legendLines = []
  if "14" in name:
      r1=3e2
      r2=1e5

  mini=9999.
  maxi=0.

  histo = infile.Get(name)
  histo.SetDirectory(0)
  histos.append(histo)
  legendLines.append(leg)
  # Close the input file
  infile.Close()

  # Make a canvas to put the plot on.
  # We don't want log axes for this plot, 
  # but if you do you can control them here.
  c = ROOT.TCanvas("canvas",'',0,0,1024,768)
  c.SetLogx(False)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what x and y range to use in the display.
  xRange = [-3.2,3.2]
  yRange = [0.9,1.1]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kOrange,ROOT.kBlack]

  # Make a legend.
  # These are the locations of the left side, bottom side, right
  # side, and top, as fractions of the canvas.
  legend = ROOT.TLegend(0.6,0.72,0.92,0.92)
  # Make the text a nice fond, and big enough
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  # A few more formatting things .....
  legend.SetBorderSize(0)
  legend.SetLineColor(0)
  legend.SetLineStyle(1)
  legend.SetLineWidth(3)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

  # Draw each histogram.
  # You really shouldn't put two histograms with different
  # x axes on the same plot - I'm only doing it here
  # to show you how to draw multiple plots on the same
  # canvas.
  for histo, line in zip(histos,legendLines) :

    index = histos.index(histo)
    colour = goodColours[index]

    # Set up marker to look nice
    histo.SetMarkerColor(colour)
    histo.SetMarkerSize(0)  
    histo.SetMarkerStyle(20+index)

    # Set up line to look nice
    histo.SetLineColor(colour)
    histo.SetLineWidth(3)
    histo.SetLineStyle(1)

    # Make sure we don't get a fill
    histo.SetFillColor(0)

    # Label my axes!!
    histo.GetXaxis().SetTitle("ADC")
    histo.GetYaxis().SetTitle("Fraction")
    # Move the label around if you want
    histo.GetYaxis().SetTitleOffset(1.5)

    # Set the limit for the axes
    #r1=histo.GetBinCenter(histo.GetMinimumBin())
    #r2=histo.GetBinCenter(histo.GetMaximumBin())*50.
    #histo.GetXaxis().SetRangeUser(r1,r2)
    #histo.GetYaxis().SetRangeUser(r1,r2)

    if index==0 :
      histo.Draw("COLZ") # Draw data points (you'll get error bars by default)
    else :
      histo.Draw("H SAME") # SAME means don't get rid of the previous stuff on the canvas

    # Fill entry into legend
    # "PL" means both the line and point style
    # will show up in the legend.
    legend.AddEntry(histo,line,"PL")

  # Actually draw the legend
  legend.Draw()

  # This is one way to draw text on the plot
  myLatex = ROOT.TLatex()
  myLatex.SetTextColor(ROOT.kBlack)
  myLatex.SetNDC()

  # Put an VLAST-P Internal label
  # I think it has to be Helvetica
  myLatex.SetTextSize(0.05)
  myLatex.SetTextFont(72)
  # These are the x and y coordinates of the bottom left corner of the text
  # as fractions of the canvas
  myLatex.DrawLatex(0.58,0.68,"VLAST-P")
  # Now we switch back to normal font for the "Internal"
  myLatex.SetTextFont(42)
  myLatex.DrawLatex(0.75,0.68,txt)

  # Update the canvas
  c.Update()

  # Save the output as a .eps, a .C, and a .root
  c.SaveAs(name+leg+txt+".png")

def plot_single(name,leg,txt,r1,r2):
  # Load some histos from the example file
  # (Gaussian limits from 2016 TLA conf)
  infile = ROOT.TFile.Open("hist.root","READ")
  histos = []
  legendLines = []
  if "14" in name:
      r1=3e2
      r2=1e5

  mini=9999.
  maxi=0.

  histo = infile.Get(name)
  histo.SetDirectory(0)
  if(histo.Integral()>0.):
      histo.Scale(1./histo.Integral())
  histos.append(histo)
  legendLines.append(leg)
  # Close the input file
  infile.Close()

  # Make a canvas to put the plot on.
  # We don't want log axes for this plot, 
  # but if you do you can control them here.
  c = ROOT.TCanvas("canvas",'',0,0,1024,768)
  c.SetLogx(True)
  c.SetLogy(False)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what x and y range to use in the display.
  xRange = [-3.2,3.2]
  yRange = [0.9,1.1]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kOrange,ROOT.kBlack]

  # Make a legend.
  # These are the locations of the left side, bottom side, right
  # side, and top, as fractions of the canvas.
  legend = ROOT.TLegend(0.6,0.72,0.92,0.92)
  # Make the text a nice fond, and big enough
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  # A few more formatting things .....
  legend.SetBorderSize(0)
  legend.SetLineColor(0)
  legend.SetLineStyle(1)
  legend.SetLineWidth(3)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

  # Draw each histogram.
  # You really shouldn't put two histograms with different
  # x axes on the same plot - I'm only doing it here
  # to show you how to draw multiple plots on the same
  # canvas.
  for histo, line in zip(histos,legendLines) :

    index = histos.index(histo)
    colour = goodColours[index]

    # Set up marker to look nice
    histo.SetMarkerColor(colour)
    histo.SetMarkerSize(0)  
    histo.SetMarkerStyle(20+index)

    # Set up line to look nice
    histo.SetLineColor(colour)
    histo.SetLineWidth(3)
    histo.SetLineStyle(1)

    # Make sure we don't get a fill
    histo.SetFillColor(0)

    # Label my axes!!
    histo.GetXaxis().SetTitle("ADC")
    histo.GetYaxis().SetTitle("Fraction")
    # Move the label around if you want
    histo.GetYaxis().SetTitleOffset(1.5)

    # Set the limit for the axes
    #r1=histo.GetBinCenter(histo.GetMinimumBin())
    #r2=histo.GetBinCenter(histo.GetMaximumBin())*50.
    histo.GetXaxis().SetRangeUser(r1,r2)
    histo.GetYaxis().SetRangeUser(0.0,histo.GetMaximum()*1.4)

    if index==0 :
      histo.Draw("H") # Draw data points (you'll get error bars by default)
    else :
      histo.Draw("H SAME") # SAME means don't get rid of the previous stuff on the canvas

    # Fill entry into legend
    # "PL" means both the line and point style
    # will show up in the legend.
    legend.AddEntry(histo,line,"PL")

  # Actually draw the legend
  legend.Draw()

  # This is one way to draw text on the plot
  myLatex = ROOT.TLatex()
  myLatex.SetTextColor(ROOT.kBlack)
  myLatex.SetNDC()

  # Put an VLAST-P Internal label
  # I think it has to be Helvetica
  myLatex.SetTextSize(0.05)
  myLatex.SetTextFont(72)
  # These are the x and y coordinates of the bottom left corner of the text
  # as fractions of the canvas
  myLatex.DrawLatex(0.58,0.68,"VLAST-P")
  # Now we switch back to normal font for the "Internal"
  myLatex.SetTextFont(42)
  myLatex.DrawLatex(0.75,0.68,txt)

  # Update the canvas
  c.Update()

  # Save the output as a .eps, a .C, and a .root
  c.SaveAs(name+leg+txt+".png")
def plot_all():
  # Load some histos from the example file
  # (Gaussian limits from 2016 TLA conf)
  infile = ROOT.TFile.Open("hist.root","READ")
  histos = []
  legendLines = ["Main HG","Main LG","Back HG","Back LG"]
  hname = ["hpmh","hpml","hpbh","hpbl"]

  mini=9999.
  maxi=0.

  for i in hname:
    histo = infile.Get(i)
    histo.SetDirectory(0)
    if(histo.Integral()>0.):
        histo.Scale(1./histo.Integral())
    histos.append(histo)
  # Close the input file
  infile.Close()

  # Make a canvas to put the plot on.
  # We don't want log axes for this plot, 
  # but if you do you can control them here.
  c = ROOT.TCanvas("canvas",'',0,0,1024,768)
  c.SetLogx(True)
  c.SetLogy(True)
  c.SetGridx(0)
  c.SetGridy(0)

  # Decide what x and y range to use in the display.
  xRange = [-3.2,3.2]
  yRange = [0.9,1.1]

  # Decide what colours to use.
  # These ones look decent, but obviously use
  # whatever you like best.
  goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kOrange,ROOT.kBlack]

  # Make a legend.
  # These are the locations of the left side, bottom side, right
  # side, and top, as fractions of the canvas.
  legend = ROOT.TLegend(0.6,0.72,0.92,0.92)
  # Make the text a nice fond, and big enough
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  # A few more formatting things .....
  legend.SetBorderSize(0)
  legend.SetLineColor(0)
  legend.SetLineStyle(1)
  legend.SetLineWidth(3)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

  # Draw each histogram.
  # You really shouldn't put two histograms with different
  # x axes on the same plot - I'm only doing it here
  # to show you how to draw multiple plots on the same
  # canvas.
  for histo, line in zip(histos,legendLines) :

    index = histos.index(histo)
    colour = goodColours[index]

    # Set up marker to look nice
    histo.SetMarkerColor(colour)
    histo.SetMarkerSize(1)  
    histo.SetMarkerStyle(20+index)

    # Set up line to look nice
    histo.SetLineColor(colour)
    histo.SetLineWidth(3)
    histo.SetLineStyle(1)

    # Make sure we don't get a fill
    histo.SetFillColor(0)

    # Label my axes!!
    histo.GetXaxis().SetTitle("ADC")
    histo.GetYaxis().SetTitle("Fraction")
    # Move the label around if you want
    histo.GetYaxis().SetTitleOffset(1.5)

    # Set the limit for the axes
    #histo.GetXaxis().SetLimits(xRange[0],xRange[1])
    histo.GetYaxis().SetRangeUser(0.001,0.5)

    if index==0 :
      histo.Draw("H") # Draw data points (you'll get error bars by default)
    else :
      histo.Draw("H SAME") # SAME means don't get rid of the previous stuff on the canvas

    # Fill entry into legend
    # "PL" means both the line and point style
    # will show up in the legend.
    legend.AddEntry(histo,line,"PL")

  # Actually draw the legend
  legend.Draw()

  # This is one way to draw text on the plot
  myLatex = ROOT.TLatex()
  myLatex.SetTextColor(ROOT.kBlack)
  myLatex.SetNDC()

  # Put an VLAST-P Internal label
  # I think it has to be Helvetica
  myLatex.SetTextSize(0.05)
  myLatex.SetTextFont(72)
  # These are the x and y coordinates of the bottom left corner of the text
  # as fractions of the canvas
  myLatex.DrawLatex(0.18,0.88,"VLAST-P")
  # Now we switch back to normal font for the "Internal"
  myLatex.SetTextFont(42)
  myLatex.DrawLatex(0.35,0.88,"Pedestal")

  # Update the canvas
  c.Update()

  # Save the output as a .eps, a .C, and a .root
  c.SaveAs("Plat.png")


if __name__ == "__main__":
  #plot_ct("c")
  #plot_ct("t")
  plot_all()
  plot_sub("hpmh")
  plot_sub("hpml")
  plot_sub("hpbh")
  plot_sub("hpbl")
  plot_sig("hmh")
  plot_sig("hml")
  plot_sig("hbh")
  plot_sig("hbl")
  chns = [2,12,14,24]
  for i in chns:
    plot_single("hmh_"+str(i),"Main_HG",str(i),3e2,5e4)
    plot_single("hml_"+str(i),"Main_LG",str(i),30,5e3)
    plot_single("hbh_"+str(i),"Back_HG",str(i),3e2,5e4)
    plot_single("hbl_"+str(i),"Back_LG",str(i),30,5e3)
    plot_2dsingle("hmainr_"+str(i),"Main_HLRatio",str(i),1,1)
    plot_2dsingle("hbackr_"+str(i),"Back_HLRatio",str(i),1,1)
    plot_2dsingle("hhighr_"+str(i),"High_MBRatio",str(i),1,1)
    plot_2dsingle("hlowr_"+str(i),"Low_MBRatio",str(i),1,1)
