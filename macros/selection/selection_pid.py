import os
import sys
import pandas as pd
from tabulate import tabulate
from ROOT import gStyle, TFile, TCanvas, TPad, TGaxis
from ROOT import TCanvas, TH2F, TLegend

from pathlib import Path
macroDir = Path(__file__).parent.parent.absolute()
libDir   = os.path.join(macroDir, 'lib')
sys.path.append(os.path.join(macroDir, 'lib'))
from processDict import production

def normHist(h):
  h.Scale(1.0 / h.Integral())

def drawHist(c, h, logy=False, option=""):
  c.cd()
  if logy:
    c.SetLogy()
  h.Draw(option)
  c.Update()

def addHist(file, hpath, htype, hsum):
  h = file.Get(f'{hpath}_{htype}')
  hsum.Add(h)

def to_latex(s):
    # Convert the TString to a Python string
    s = str(s)
    label = f"{s[0]}#bar{{{s[0]}}}"

    # Return the LaTeX representation
    return label

def styleHist(h, isty):
  # normHist(h)
  h.SetFillColor(0)
  h.SetLineWidth(3)
  h.SetLineStyle(isty)  # Different line style for each histogram


def styleLegend(leg):
  leg.SetLineColor(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.0001)
  leg.SetMargin(0.8)

def stylePad(pad, top=0, buttom=0, left=0, right=0):
  # pad.SetGrid(1,1)
  pad.SetTopMargin(top)
  pad.SetBottomMargin(buttom)
  pad.SetLeftMargin(left)
  pad.SetRightMargin(right)
  pad.Draw()
  pad.cd()

def main():

  inDir = os.path.join(macroDir, '..', 'rootfiles', 'merged')
  # inDir = os.path.join(macroDir, '..', 'rootfiles', 'offline')

  chirals    = ["eLpR"]
  processes  = ["P2f_z_h"]
  qqbars     = ["dd", "uu", "ss", "cc", "bb"]
  categories = ["PID"]

  files = {}
  for process in processes:
    for chiral in chirals:
      chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
      prochiTuple = (process, chiral)
      processID  = production[process, chiral][0]
      filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
      files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))


  legend = TLegend(0.60,0.75,0.88,0.88)

  for process in processes:

    if process == "P2f_z_h": # signal

      for category in categories:

        hprefix = f'h_{process}_{category}'
        h_sum = TH2F(f"{hprefix}_sum", f"{hprefix}_sum", 4, 0, 4, 4, 0, 4)

        for chiral in chirals:
          file = files[process, chiral]
        for qqbar in qqbars:
          hpath = f'{qqbar}/efficiency'
          hname_prefix = f'h_{qqbar}_reco_Pi'
          h = file.Get(f'{hpath}/{hname_prefix}_{category}_offset')
          h_sum.Add(h)

  TGaxis.SetMaxDigits(3)
  gStyle.SetOptStat(0)

  c = TCanvas("c", "c", 900, 900)
  c.cd()
  pad = TPad("pad", "pad", 0, 0, 1, 1)
  stylePad(pad,0.1,0.1,0.12,0.1)
  h_sum_plot = h_sum.Clone()
  # h_sum_plot.Scale(1.0 / h_sum_plot.GetEntries())
  h_sum_plot.SetTitle(";Reconstructed;Truth")
  h_sum_plot.GetXaxis().SetRangeUser(0, 3)
  h_sum_plot.GetYaxis().SetRangeUser(1, 4)

  h_sum_plot.GetXaxis().SetBinLabel(1, "#pi^{#pm}")
  h_sum_plot.GetXaxis().SetBinLabel(2, "K^{#pm}")
  h_sum_plot.GetXaxis().SetBinLabel(3, "p (#bar{p})")

  h_sum_plot.GetYaxis().SetBinLabel(4, "#pi^{#pm}")
  h_sum_plot.GetYaxis().SetBinLabel(3, "K^{#pm}")
  h_sum_plot.GetYaxis().SetBinLabel(2, "p (#bar{p})")

  h_sum_plot.GetXaxis().SetLabelSize(0.05)
  h_sum_plot.GetYaxis().SetLabelSize(0.05)

  h_sum_plot.GetXaxis().SetTitleOffset(1.2)
  h_sum_plot.GetYaxis().SetTitleOffset(1.2)

  h_sum_plot.Draw("col text")
  c.SaveAs(f"plots/c_PID.pdf")


  # purity
  total_gen = [0,0,0] # K, pi, p
  for xbin in range(3):
    for ybin in reversed(range(4)):
      total_gen[xbin] += h_sum.GetBinContent(xbin+1, ybin+1)


  # print(total_gen)
  # print(h_sum.GetBinContent(2,1))

  purity = [ [h_sum.GetBinContent(x+1,y+1) / total_gen[x] for y in reversed(range(4))] for x in range(3) ]

  print(tabulate(purity, tablefmt='latex',floatfmt=".3f"))
  # print(*purity,sep='\n')

  # efficiency
  total_reco = [0,0,0] # K, pi, p
  for ybin in reversed(range(1,4)):
    for xbin in range(3):
      total_reco[xbin] += h_sum.GetBinContent(xbin+1, ybin+1)

  print(total_reco)

  efficiency = [ h_sum.GetBinContent(i+1,4-i) / total_reco[i] for i in range(3) ]

  print(f"efficiency of pi: {efficiency[0]}")
  print(f"efficiency of K: {efficiency[1]}")
  print(f"efficiency of p: {efficiency[2]}")

if __name__ == "__main__":
  main()