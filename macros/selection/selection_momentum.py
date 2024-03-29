import os
import sys
import pandas as pd
from ROOT import gStyle, TFile, TCanvas, TPad
from ROOT import TCanvas, TH1F, THStack, TLine, TLegend

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
  pad.SetGrid(1,1)
  pad.SetTopMargin(top)
  pad.SetBottomMargin(buttom)
  pad.SetLeftMargin(left)
  pad.SetRightMargin(right)
  pad.Draw()
  pad.cd()

def main():

  inDir = os.path.join(macroDir, '..', 'rootfiles', 'merged')

  chirals    = ["eLpR"]
  processes  = ["P2f_z_h"]
  qqbars     = ["dd", "uu", "ss"]
  genreco    = ["reco"]
  categories = ["momentum_LPFO", "momentum_SPFO"]

  xnbins = 120
  xmin   = 0
  xmax   = 120

  files = {}
  for process in processes:
    for chiral in chirals:
      chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
      prochiTuple = (process, chiral)
      processID  = production[process, chiral][0]
      filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
      files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))


  legend = TLegend(0.60,0.75,0.88,0.88)

  ths_LPFO_SPFO = []

  for process in processes:

    if process == "P2f_z_h": # signal

      for igenreco in genreco:

        tmp_ths_LPFO_SPFO = THStack(f"ths_{process}_{igenreco}", f"ths_{process}_{igenreco}")

        for category in categories:

          hprefix = f'h_{process}_{igenreco}_{category}'
          h_sum = TH1F(f"{hprefix}_sum", f"{hprefix}_sum", xnbins, xmin, xmax)

          for chiral in chirals:
            file = files[process, chiral]
          for qqbar in qqbars:
            hpath = f'{qqbar}/efficiency'
            hname_prefix = f'h_{qqbar}_{igenreco}_Pi'
            h = file.Get(f'{hpath}/{hname_prefix}_{category}_nocut')
            h_sum.Add(h)

          styleHist(h_sum,1)
          # legend.AddEntry(h_sum, category.split("_")[1], "l")
          if category == "momentum_LPFO":
            legend.AddEntry(h_sum, "LPFO", "l")
          elif category == "momentum_SPFO":
            legend.AddEntry(h_sum, "non-LPFO", "l")

          h_sum.SetTitle(f"{process} {category} {igenreco}")
          tmp_ths_LPFO_SPFO.Add(h_sum)
        
        ths_LPFO_SPFO.append(tmp_ths_LPFO_SPFO)


  h_leg_gen  = TH1F("h_leg_gen", "h_leg_gen", 120, 0, 120)
  styleHist(h_leg_gen,1)
  # legend.AddEntry(h_leg_gen, "Parton level", "l")

  h_leg_reco = TH1F("h_leg_reco", "h_leg_reco", 120, 0, 120)
  styleHist(h_leg_reco,1)
  h_leg_reco.SetLineStyle(2)
  # legend.AddEntry(h_leg_reco, "Reconstructed", "l")

  c = TCanvas("c", "c", 900, 900)
  c.cd()
  pad = TPad("pad", "pad", 0, 0, 1, 1)
  stylePad(pad,0.1,0.1,0.15,0.1)
  gStyle.SetPalette(51)
  ths_LPFO_SPFO[0].SetTitle(f";Momentum [GeV];Entries")
  ths_LPFO_SPFO[0].Draw("h plc nostack")
  # ths_LPFO_SPFO[1].Draw("h plc nostack same")
  ths_LPFO_SPFO[0].GetXaxis().SetLimits(0, 120)
  legend.Draw()
  pad.SetLogy()
  c.SaveAs(f"plots/c_LPFO_SPFO.pdf")



if __name__ == "__main__":
  main()