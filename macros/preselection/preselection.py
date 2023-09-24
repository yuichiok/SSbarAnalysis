import os
import sys
from ROOT import TFile, TCanvas, TLegend
from ROOT import TCanvas, TH1F, THStack
import subprocess
import argparse

from pathlib import Path
macroDir = Path(__file__).parent.parent.absolute()
libDir   = os.path.join(macroDir, 'lib')
sys.path.append(os.path.join(macroDir, 'lib'))
from processDict import production

inDir = os.path.join(macroDir, '..', 'rootfiles', 'merged')

chirals   = ["eLpR", "eRpL"]
processes = ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"]
qqbars    = ["dd", "uu", "ss", "cc", "bb"]

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
  normHist(h)
  hsum.Add(h)  

files = dict()
for process in processes:
  for chiral in chirals:
    chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
    prochiTuple = (process, chiral)
    processID  = production[process, chiral][0]
    filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
    files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))

class PlotManager:
  def __init__(self, category, bin, xmin, xmax):
    self.category = category
    self.bin = bin
    self.xmin = xmin
    self.xmax = xmax
    self.canvas = TCanvas(f"c_{category}", f"c_{category}", 900, 900)
    self.stack  = THStack(f"ths_{category}", f"ths_{category}")


PMs = {
  "sinacol": PlotManager("sinacol", 100, 0, 1),
  "invM": PlotManager("invM", 100, 0, 500),
  "y23": PlotManager("y23", 50, 0, 0.25)
}
for process in processes:

  if process == "P2f_z_h": # signal
    for qqbar in qqbars:

      for category, PM in PMs.items():

        hprefix = f'h_{process}_{qqbar}_{category}'
        h_sum = TH1F(f"{hprefix}_sum", f"{hprefix}_sum", PM.bin, PM.xmin, PM.xmax)

        for chiral in chirals:
          chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"

          file = files[process, chiral]
          
          tmp_hname = f'{qqbar}/preselection/h_{qqbar}'
          addHist(file, tmp_hname, category, h_sum)

        PM.stack.Add(h_sum)

  else: # background
    for category, PM in PMs.items():

      hprefix = f'h_{process}_{category}'
      h_sum = TH1F(f"{hprefix}_sum", f"{hprefix}_sum", PM.bin, PM.xmin, PM.xmax)

      for chiral in chirals:
        chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"

        file = files[process, chiral]
        
        tmp_hname = f'bg/preselection/h_bg'
        addHist(file, tmp_hname, category, h_sum)

      PM.stack.Add(h_sum)

for category, PM in PMs.items():
  PM.canvas.cd()
  PM.stack.Draw("h nostack")
  PM.canvas.SaveAs(f"c_{category}.pdf")
