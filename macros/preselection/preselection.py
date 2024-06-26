import os
import sys
import pandas as pd
from tabulate import tabulate
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

def addData(file, hpath, df, process, qqbar, chiral):
  hcategory = ["cosBF","sinacol", "invM", "y23","cosAF"]
  # Initialize dictionaries to accumulate values
  values_dict = {
      "process": process,
      "qqbar": qqbar,
      "chiral": chiral,
  }

  for category in hcategory:
    h = file.Get(f'{hpath}_{category}')
    entry = (h.Integral() / 2.0) if (category == "cosBF") or (category == "cosAF") else h.Integral()
    values_dict[category] = entry

  # Append the values to a temporary DataFrame
  temp_df = pd.DataFrame([values_dict])
  
  # Concatenate the temporary DataFrame with the main DataFrame
  df = pd.concat([df, temp_df], ignore_index=True)
  return df

def to_latex(s):
    bg_label = {
      "P4f_ww_h": "WW",
      "P4f_zz_h": "ZZ",
      "Pqqh": "qqH",
    }
    # Convert the TString to a Python string
    s = str(s)
    label = f"{s[0]}#bar{{{s[0]}}}" if len(s) == 2 else bg_label[s]
    
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

  chirals   = ["eLpR", "eRpL"]
  processes = ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"]
  qqbars    = ["dd", "uu", "ss", "cc", "bb", "rr"]

  files = {}
  for process in processes:
    for chiral in chirals:
      chiralDot = "eL.pR" if chiral == "eLpR" else "eR.pL"
      prochiTuple = (process, chiral)
      processID  = production[process, chiral][0]
      filename = f"rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root"
      files[prochiTuple] = TFile.Open(os.path.join(inDir, filename))

  class PlotManager:
    cuts = {
      "sinacol" : 0.3,
      "invM"    : 140,
      "y23"     : 0.02,
    }
    def __init__(self, category, title, bin, xmin, xmax):
      self.category = category
      self.title = title
      self.bin = bin
      self.xmin = xmin
      self.xmax = xmax
      self.canvas = TCanvas(f"c_{category}", f"c_{category}", 900, 900)
      self.stack  = THStack(f"ths_{category}", title)
      self.legend = TLegend(0.60,0.5,0.75,0.85) if category != "y23" else TLegend(0.75,0.55,0.90,0.9)
      self.styleLegend()
      self.line = None

    def styleLegend(self):
      self.legend.SetLineColor(0)
      self.legend.SetFillColor(0)
      self.legend.SetFillStyle(0)
      self.legend.SetTextSize(0.03)
      self.legend.SetMargin(0.8)

    def styleLine(self):
      self.canvas.Range(self.xmin,0,self.xmax,1)
      self.line = TLine(self.cuts[category], self.stack.GetMinimum("nostack"), self.cuts[category], self.stack.GetMaximum("nostack")*1.05)
      self.line.SetLineColor(2)
      self.line.SetLineWidth(3)
      self.line.SetLineStyle(1)

  columnDf = ["process", "qqbar", "chiral", "cosBF", "sinacol", "invM", "y23", "cosAF"]
  totDf = pd.DataFrame(columns=columnDf)
  PMs = {
    "sinacol": PlotManager("sinacol", ";sin#Psi_{acol};Norm.", 100, 0, 1),
    "invM": PlotManager("invM", ";m_{j_{1},j_{2}};Norm.", 100, 0, 500),
    "y23": PlotManager("y23", ";y_{23};Norm.", 50, 0, 0.25),
  }

  for process in processes:

    if process == "P2f_z_h": # signal
      for qqbar in qqbars:

        for chiral in chirals:
          hpath = f'{qqbar}/preselection/h_{qqbar}'
          totDf = addData(files[process, chiral], hpath, totDf, process, qqbar, chiral)

        for category, PM in PMs.items():

          hprefix = f'h_{process}_{qqbar}_{category}'
          h_sum = TH1F(f"{hprefix}_sum", f"{hprefix}_sum", PM.bin, PM.xmin, PM.xmax)

          for chiral in chirals:
            file = files[process, chiral]
            tmp_hname = f'{qqbar}/preselection/h_{qqbar}'
            addHist(file, tmp_hname, category, h_sum)

          styleHist(h_sum,1)
          h_sum.SetTitle(f"{process} {qqbar}")
          PM.stack.Add(h_sum)
          PM.legend.AddEntry(h_sum, to_latex(qqbar), "l")

    else: # background
      for chiral in chirals:
        hpath = f'bg/preselection/h_bg'
        totDf = addData(files[process, chiral], hpath, totDf, process, 'bg', chiral)

      for category, PM in PMs.items():

        hprefix = f'h_{process}_{category}'
        h_sum = TH1F(f"{hprefix}_sum", f"{hprefix}_sum", PM.bin, PM.xmin, PM.xmax)

        for chiral in chirals:
          file = files[process, chiral]
          tmp_hname = f'bg/preselection/h_bg'
          addHist(file, tmp_hname, category, h_sum)

        styleHist(h_sum,7)
        h_sum.SetTitle(f"{process}")
        PM.stack.Add(h_sum)
        PM.legend.AddEntry(h_sum, to_latex(process), "l")

  # Calculate the efficiencies
  effDf = pd.DataFrame(columns=columnDf[:2])
  cutno = 1
  for column in columnDf:
    if column == "process" or column == "qqbar" or column == "chiral":
      effDf[column] = totDf[column]
    elif column == "cosBF":
      denom = totDf[column]
      effDf['None'] = denom
      # print(totDf["process"], totDf["qqbar"], totDf["chiral"], denom)
    else:
      effDf[f'cut{cutno}'] = totDf[column] / denom * 100.
      cutno += 1
  with pd.option_context('display.float_format', '{:0.1f}'.format):
    effDf = effDf.sort_values(by=['chiral', 'process'])
    eLpRdf = effDf.head(9).T
    eRpLdf = effDf.loc[1:].T
    # print(eLpRdf)
    # print(eRpLdf)

    tmp = eRpLdf.tail(5)

    print(tabulate(tmp, tablefmt='latex',floatfmt=".1f",showindex=False))



    # file_eLpR = Path('eLpR.csv')  
    file_eLpR = os.path.join(macroDir, 'preselection', 'eLpR.csv')
    file_eRpL = os.path.join(macroDir, 'preselection', 'eRpL.csv')
    eLpRdf.to_csv(file_eLpR)
    eRpLdf.to_csv(file_eRpL)

  # Draw the histograms
  for category, PM in PMs.items():
    c = PM.canvas
    c.cd()
    pad = TPad("pad", "pad", 0, 0, 1, 1)
    stylePad(pad,0.1,0.1,0.15,0.1)
    gStyle.SetPalette(55)
    PM.stack.Draw("h plc nostack")
    PM.legend.Draw()
    PM.styleLine()
    if category == "y23":
      pad.SetLogy()
    PM.line.Draw()
    c.SaveAs(f"c_{category}.pdf")

if __name__ == "__main__":
  main()