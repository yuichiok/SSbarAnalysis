#include <iostream>
#include "include/Styles.hh"

void Normalize(TH1F *h)
{
  h->Scale( 1.0 / h->GetEntries() );
}

void LPFO_sep_cos()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.eff.hists.all.root","READ");

  TH1F *h_reco_K_SLPFO_mom_diff  = (TH1F*) file->Get("h_lpfo_reco_K_sep_cos");
  TH1F *h_reco_Pi_SLPFO_mom_diff = (TH1F*) file->Get("h_lpfo_reco_Pi_sep_cos");
  h_reco_K_SLPFO_mom_diff->Sumw2();
  h_reco_Pi_SLPFO_mom_diff->Sumw2();

  Normalize(h_reco_K_SLPFO_mom_diff);
  Normalize(h_reco_Pi_SLPFO_mom_diff);

  StyleHist(h_reco_K_SLPFO_mom_diff,kBlack);
  StyleHist(h_reco_Pi_SLPFO_mom_diff,kBlue);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  StylePad(gPad,0,0.15,0,0.17);

  h_reco_K_SLPFO_mom_diff->GetYaxis()->SetRangeUser(0,1);
  h_reco_K_SLPFO_mom_diff->GetXaxis()->SetRangeUser(-1,-0.6);

  h_reco_K_SLPFO_mom_diff->SetTitle("LPFO cos#theta Separation;cos#theta_{sep}; Entries(Norm.)");
  h_reco_K_SLPFO_mom_diff->Draw("h");
  h_reco_Pi_SLPFO_mom_diff->Draw("hsame");

  TLegend *leg_trp = new TLegend(0.3,0.8,0.7,0.85);
  leg_trp->SetMargin(0.4);
  leg_trp->SetLineColor(0);
  leg_trp->AddEntry(h_reco_K_SLPFO_mom_diff,"K^{#pm}","l");
  leg_trp->AddEntry(h_reco_Pi_SLPFO_mom_diff,"#pi^{#pm}","l");
  leg_trp->Draw();

  c0->Draw();

}