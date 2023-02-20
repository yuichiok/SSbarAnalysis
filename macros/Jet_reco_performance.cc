#include <iostream>
#include "Styles.cc"

template <class H2>
void Normalize(H2 *h)
{
  h->Scale( 1.0 / h->GetEntries() );
}

void StyleHist(TH1F *h, Color_t col)
{
  h->SetLineWidth(3);
  h->SetLineColor(col);
  h->SetFillStyle(3002);
  h->SetFillColor(col);
}

void Jet_reco_performance()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  TH1F *h2_jet_mult_cos       = (TH1F*) file->Get("jet/h2_jet_mult_cos");
  TH1F *h2_jet_mult_cos_noISR = (TH1F*) file->Get("jet/h2_jet_mult_cos_noISR");

  // Normalize(h2_jet_mult_cos);
  // Normalize(h2_jet_mult_cos_noISR);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TPad *pad0 = new TPad("pad0", "pad0",0,0,1,1);
  StylePad(pad0,0,0.12,0.15,0.15);

  h2_jet_mult_cos->Draw("colz");
  c0->Draw();

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TPad *pad1 = new TPad("pad1", "pad1",0,0,1,1);
  StylePad(pad1,0,0.12,0.15,0.15);

  h2_jet_mult_cos_noISR->Draw("colz");
  c1->Draw();

}