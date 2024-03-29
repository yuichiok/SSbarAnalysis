import os
import sys
import pandas as pd
from ROOT import gStyle, TFile, TCanvas, TPad, kBlack, kBlue
from ROOT import TCanvas, TH1F, THStack, TLine, TLegend, TGaxis

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

def styleLegend(leg):
  leg.SetLineColor(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.02)
  leg.SetMargin(0.4)

def to_latex(s):
    # Convert the TString to a Python string
    s = str(s)
    label = f"{s[0]}#bar{{{s[0]}}}"

    # Return the LaTeX representation
    return label

def styleHist(h):
  h.SetLineWidth(3)

def stylePad(pad, top=0, buttom=0, left=0, right=0):
  pad.SetGrid(1,1)
  pad.SetTopMargin(top)
  pad.SetBottomMargin(buttom)
  pad.SetLeftMargin(left)
  pad.SetRightMargin(right)
  pad.Draw()
  pad.cd()

def main():

  TGaxis.SetMaxDigits(3)
  gStyle.SetOptStat(0)

  inDir = os.path.join(macroDir, '..', 'rootfiles', 'merged')

  chirals    = ["eLpR","eRpL"]
  processes  = ["P2f_z_h"]
  LPFO_modes = ["K","Pi"]
  qqbars     = ["dd", "uu", "ss", "cc", "bb"]

  files = {}
  for process in processes:
    for chiral in chirals:
      chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
      prochiTuple = (process, chiral)
      processID  = production[process, chiral][0]
      filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
      files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))


  for LPFO_mode in LPFO_modes:

    c_gen_reco_N_cos = TCanvas(f"c_{LPFO_mode}_gen_reco_N_cos",f"{LPFO_mode}_gen_reco_N_cos",800,800)
    c_gen_reco_N_cos.cd()
    h_gen_N_cos_sum  = TH1F(f"h_{LPFO_mode}_gen_N_cos_sum",";cos#theta;Entries",100,-1,1)
    h_reco_N_cos_sum = TH1F(f"h_{LPFO_mode}_reco_N_cos_sum",";cos#theta;Entries",100,-1,1)

    for qqbar in qqbars:

      hpath = f'{qqbar}/resolution'
      hname_gen_prefix = f'h_{qqbar}_{LPFO_mode}_gen_N_cos'
      hname_reco_prefix = f'h_{qqbar}_{LPFO_mode}_reco_N_cos'

      for chiral in chirals:
        file  = files["P2f_z_h", chiral]
        hgen  = file.Get(f'{hpath}/{hname_gen_prefix}')
        hreco = file.Get(f'{hpath}/{hname_reco_prefix}')
        h_gen_N_cos_sum.Add(hgen)
        h_reco_N_cos_sum.Add(hreco)

    styleHist(h_gen_N_cos_sum)
    styleHist(h_reco_N_cos_sum)
    h_gen_N_cos_sum.SetLineColor(kBlack)
    h_reco_N_cos_sum.SetLineColor(kBlue)

    pad = TPad(f"{LPFO_mode}_pad", f"{LPFO_mode}_pad", 0, 0, 1, 1)
    stylePad(pad,0.1,0.1,0.15,0.1)
    h_gen_N_cos_sum.GetYaxis().SetRangeUser(0,1500E3)
    h_gen_N_cos_sum.Draw("h")
    h_reco_N_cos_sum.Draw("hsame")

    legend = TLegend(0.2,0.15,0.5,0.3)
    styleLegend(legend)
    legend.AddEntry(h_gen_N_cos_sum,"Number of generated PFO","l")
    legend.AddEntry(h_reco_N_cos_sum,"Number of reconstructed PFO","l")
    legend.Draw("same")
    pad.Draw()
    c_gen_reco_N_cos.Draw()
    c_gen_reco_N_cos.SaveAs(f"plots/c_{LPFO_mode}.pdf")

    c_gen_reco_cos_ratio = TCanvas(f"c_{LPFO_mode}_gen_reco_cos_ratio",f"{LPFO_mode}_gen_reco_cos_ratio",800,800)
    p_gen_reco_cos_ratio = TPad(f"p_{LPFO_mode}_gen_reco_cos_ratio", f"p_{LPFO_mode}_gen_reco_cos_ratio", 0, 0, 1, 1)
    stylePad(p_gen_reco_cos_ratio,0.1,0.1,0.15,0.1)
    heff = h_reco_N_cos_sum.Clone()
    heff.Divide(h_gen_N_cos_sum)
    heff.Draw("")
    p_gen_reco_cos_ratio.Draw()
    c_gen_reco_cos_ratio.Draw()
    c_gen_reco_cos_ratio.SaveAs(f"plots/c_{LPFO_mode}_gen_reco_cos_ratio.pdf")






if __name__ == "__main__":
  main()