import os
import sys
from ROOT import TFile, TCanvas, TLegend
from ROOT import TCanvas, TH1F
import subprocess
import argparse

from pathlib import Path
macroDir = Path(__file__).parent.parent.absolute()
libDir   = os.path.join(macroDir, 'lib')
sys.path.append(os.path.join(macroDir, 'lib'))
print(macroDir)
from processDict import production

inDir = os.path.join(macroDir, '..', 'rootfiles', 'merged')

chirals   = ["eLpR", "eRpL"]
processes = ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"]
qqbars    = ["dd", "uu", "ss", "cc", "bb"]

c_sinacol = TCanvas("c_sinacol", "c_sinacol", 900, 900)
c_invM    = TCanvas("c_invM", "c_invM", 900, 900)
c_y23     = TCanvas("c_y23", "c_y23", 900, 900)

def drawHist(c, h, logy=False, option=""):
  c.cd()
  if logy:
    c.SetLogy()
  h.Draw(option)
  c.Update()

def addHist(file, hpath, htype, hsum):
  h = file.Get(f'{hpath}_{htype}')
  hsum.Add(h)  

  
files = dict()
for process in processes:
  for chiral in chirals:
    chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
    prochiTuple = (process, chiral)
    processID  = production[process, chiral][0]
    filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
    files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))

hists = []
for process in processes:

  if process == "P2f_z_h": # signal
    for qqbar in qqbars:
      hprefix = f'h_{process}_{qqbar}'
      h_sinacol_sum = TH1F(f"{hprefix}_sinacol_sum", f"{hprefix}_sinacol_sum", 100, 0, 1)
      h_invM_sum    = TH1F(f"{hprefix}_invM_sum", f"{hprefix}_invM_sum", 100, 0, 500)
      h_y23_sum     = TH1F(f"{hprefix}_y23_sum", f"{hprefix}_y23_sum", 50, 0, 0.25)

      for chiral in chirals:
        chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"

        file = files[process, chiral]

        tmp_hname = f'{qqbar}/preselection/h_{qqbar}'
        addHist(file, tmp_hname, "sinacol", h_sinacol_sum)
        addHist(file, tmp_hname, "invM"   , h_invM_sum)
        addHist(file, tmp_hname, "y23"    , h_y23_sum)

      hists.append((h_sinacol_sum, h_invM_sum, h_y23_sum))

drawCounter = 0
for iqq, hlist in zip(qqbars, hists):
  h_sinacol_sum, h_invM_sum, h_y23_sum = hlist
  drawOption = "h" if drawCounter == 0 else "hsame"
  
  c_sinacol.cd()
  h_sinacol_sum.Draw(drawOption)
  
  c_invM.cd()
  h_invM_sum.Draw(drawOption)

  c_y23.cd()
  h_y23_sum.Draw(drawOption)

  drawCounter += 1

c_invM.SaveAs("c_invM.pdf")



  


# parser = argparse.ArgumentParser(description='Preselection Analysis')
# parser.add_argument('--process', type=str, required=True,
#                     help='Production process (P2f_z_h, P4f_ww_h, P4_zz_h, Pe1e1h)')
# parser.add_argument('--chiral',  type=str, required=True,
#                     help='Polarization of beam (eLpR or eRpL or eLpL or eRpR)')

# args  = parser.parse_args()
# if args.chiral not in ["eLpR", "eRpL", "eLpL", "eRpR"]:
#   sys.exit("Error: chiral must be eLpR or eRpL or eLpL or eRpR")
# if args.process not in ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"]:
#   sys.exit("Error: process must be P2f_z_h, P4f_ww_h, P4f_zz_h, or Pe1e1h")

# chiralDot = "eL.pR" if args.chiral == "eLpR" else "eR.pL"

# processID  = production[args.process, args.chiral][0]
# prodIDList = production[args.process, args.chiral][1]

# # rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I402001.Pe1e1h.eL.pR.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root

# for prodID in prodIDList:
#   fileName = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{args.process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
