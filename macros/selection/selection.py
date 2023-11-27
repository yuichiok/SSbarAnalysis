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
  normHist(h)
  h.SetFillColor(0)
  h.SetLineWidth(3)
  h.SetLineStyle(isty)  # Different line style for each histogram

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
  qqbars     = ["dd", "uu", "ss", "cc", "bb"]
  categories ={
    "ctag" : [100,0,1],
    "btag" : [100,0,1],
    "nvtx" : [3,-0.5,2.5],
    "momentum_LPFO" : [120,0,120],
    "momentum_SPFO" : [120,0,120],
    "LPFOacol" : [100,0,1],
    "offset_non-hyperon" : [100,0,100],
  }
  cuts       = ["nocut", "SPFO"]

  files = {}
  for process in processes:
    for chiral in chirals:
      chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
      prochiTuple = (process, chiral)
      processID  = production[process, chiral][0]
      filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
      files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))

  class PlotManager:
    def __init__(self, category, title, xrange, ylog):
      self.category = category
      self.title = title
      self.bin = xrange[0]
      self.xmin = xrange[1]
      self.xmax = xrange[2]
      self.ylog = ylog
      self.canvas = TCanvas(f"c_{category}", f"c_{category}", 900, 900)
      self.stack  = THStack(f"ths_{category}", title)
      self.legend = TLegend(0.60,0.7,0.75,0.85)
      self.styleLegend()
      self.line = None

    def styleLegend(self):
      self.legend.SetLineColor(0)
      self.legend.SetFillColor(0)
      self.legend.SetFillStyle(0)
      self.legend.SetTextSize(0.027)
      self.legend.SetMargin(0.8)

  PMs = { f"{category}_{cut}": PlotManager(f"{category}_{cut}", f";{category};Norm.", xrange, 1)
          for category, xrange in categories.items()
          for cut in cuts
        }

  for process in processes:

    if process == "P2f_z_h": # signal
      for qqbar in qqbars:
        for category, PM in PMs.items():

          hpath = f'{qqbar}/efficiency'
          hname_prefix = f'h_{qqbar}_reco_Pi'

          for chiral in chirals:
            file = files[process, chiral]

          h = file.Get(f'{hpath}/{hname_prefix}_{category}')

          styleHist(h,1)
          h.SetTitle(f"{process} {qqbar}")
          PM.stack.Add(h)
          PM.legend.AddEntry(h, to_latex(qqbar), "l")

  # Draw the histograms
  for category, PM in PMs.items():
    c = PM.canvas
    c.cd()
    pad = TPad("pad", "pad", 0, 0, 1, 1)
    stylePad(pad,0.1,0.1,0.15,0.1)
    gStyle.SetPalette(55)
    PM.stack.Draw("h plc nostack")
    PM.stack.GetXaxis().SetLimits(PM.xmin, PM.xmax)
    PM.legend.Draw()
    
    xtitle = PM.stack.GetXaxis().GetTitle()
    if PM.ylog:
      pad.SetLogy()

    c.SaveAs(f"plots/c_{category}.pdf")

if __name__ == "__main__":
  main()