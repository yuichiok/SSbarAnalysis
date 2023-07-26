#include <iostream>
#include "include/Styles.hh"

TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.spfox04.eff.hists.all.root","READ");

void Normalize(TH1F *h)
{
  h->Scale( 1.0 / h->GetEntries() );
}

void Mom_LPFO_SPFO()
{
  TH2F *h_reco_Pi_SLPFOPiK_mom  = (TH2F*) file->Get("h_reco_Pi_SLPFOPiK_mom");
  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  StylePad(gPad,0,0.15,0.14,0.14);
  h_reco_Pi_SLPFOPiK_mom->Draw("colz");
}

void Mom_LSPFO_diff()
{
  TH1F *h_reco_Pi_SLPFOPiK_mom_diff  = (TH1F*) file->Get("h_reco_Pi_SLPFOPiK_mom_diff");

  StyleHist(h_reco_Pi_SLPFOPiK_mom_diff,kBlue);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  StylePad(gPad,0,0.15,0,0.17);

  h_reco_Pi_SLPFOPiK_mom_diff->SetTitle("LPFO Momentum Separation;#Delta_{p} |LPFO - SPFO|; Entries");
  h_reco_Pi_SLPFOPiK_mom_diff->Draw("h");

  c0->Draw();
}

void Mom_LSPFO_diff_sigma()
{
  TH1F *h_reco_Pi_SLPFOPiK_mom_diff_sigma  = (TH1F*) file->Get("h_reco_Pi_SLPFOPiK_mom_diff_sigma");

  StyleHist(h_reco_Pi_SLPFOPiK_mom_diff_sigma,kBlue);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  StylePad(gPad,0,0.15,0,0.17);

  h_reco_Pi_SLPFOPiK_mom_diff_sigma->SetTitle("LPFO Momentum Separation;#Delta_{p} / p_{LPFO}; Entries");
  h_reco_Pi_SLPFOPiK_mom_diff_sigma->Draw("h");

  c0->Draw();
}

void Mom_diff_sigma()
{
  TH1F *h_reco_K_SLPFO_mom_diff_sigma  = (TH1F*) file->Get("h_reco_K_SLPFO_mom_diff_sigma");
  TH1F *h_reco_Pi_SLPFO_mom_diff_sigma = (TH1F*) file->Get("h_reco_Pi_SLPFO_mom_diff_sigma");
  h_reco_K_SLPFO_mom_diff_sigma->Sumw2();
  h_reco_Pi_SLPFO_mom_diff_sigma->Sumw2();

  StyleHist(h_reco_K_SLPFO_mom_diff_sigma,kBlack);
  StyleHist(h_reco_Pi_SLPFO_mom_diff_sigma,kBlue);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  StylePad(gPad,0,0.15,0,0.17);

  h_reco_Pi_SLPFO_mom_diff_sigma->SetTitle("LPFO Momentum Separation;#Delta_{p} / p_{LPFO}; Entries");
  h_reco_Pi_SLPFO_mom_diff_sigma->Draw("h");
  h_reco_K_SLPFO_mom_diff_sigma->Draw("hsame");

  TLegend *leg_trp = new TLegend(0.3,0.8,0.7,0.85);
  leg_trp->SetMargin(0.4);
  leg_trp->SetLineColor(0);
  leg_trp->AddEntry(h_reco_K_SLPFO_mom_diff_sigma,"K^{#pm}","l");
  leg_trp->AddEntry(h_reco_Pi_SLPFO_mom_diff_sigma,"#pi^{#pm}","l");
  leg_trp->Draw();

  c0->Draw();
}

void Mom_diff()
{
  TH1F *h_reco_K_SLPFO_mom_diff  = (TH1F*) file->Get("h_reco_K_SLPFO_mom_diff");
  TH1F *h_reco_Pi_SLPFO_mom_diff = (TH1F*) file->Get("h_reco_Pi_SLPFO_mom_diff");
  h_reco_K_SLPFO_mom_diff->Sumw2();
  h_reco_Pi_SLPFO_mom_diff->Sumw2();

  StyleHist(h_reco_K_SLPFO_mom_diff,kBlack);
  StyleHist(h_reco_Pi_SLPFO_mom_diff,kBlue);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  StylePad(gPad,0,0.15,0,0.17);

  h_reco_Pi_SLPFO_mom_diff->SetTitle("LPFO Momentum Separation;#Delta_{p} |LPFO - SPFO|; Entries");
  h_reco_Pi_SLPFO_mom_diff->Draw("h");
  h_reco_K_SLPFO_mom_diff->Draw("hsame");

  TLegend *leg_trp = new TLegend(0.3,0.8,0.7,0.85);
  leg_trp->SetMargin(0.4);
  leg_trp->SetLineColor(0);
  leg_trp->AddEntry(h_reco_K_SLPFO_mom_diff,"K^{#pm}","l");
  leg_trp->AddEntry(h_reco_Pi_SLPFO_mom_diff,"#pi^{#pm}","l");
  leg_trp->Draw();

  c0->Draw();
}

void SLPFO_mom_diff()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  // Mom_LPFO_SPFO();
  Mom_LSPFO_diff();
  // Mom_LSPFO_diff_sigma();
  // Mom_diff();
  // Mom_diff_sigma();


}