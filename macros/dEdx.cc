#include <iostream>

void StylePad(TPad *pad, Float_t t, Float_t b, Float_t r, Float_t l)
{
  pad->SetGrid(1,1);
  if(t) pad->SetTopMargin(t);
  if(b) pad->SetBottomMargin(b);
  if(r) pad->SetRightMargin(r);
  if(l) pad->SetLeftMargin(l);
  pad->Draw();
  pad->cd();

}

void dEdx()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp10_pNaN.tpc0.hists.all.root","READ");

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TPad *pad0 = new TPad("pad0", "pad0",0,0,1,1);
  StylePad(pad0,0,0.15,0,0.17);

  TH1F *h2_gen_K_dEdx_p  = (TH1F*) file->Get("dEdx/gen_K_dEdx_p");
  TH1F *h2_gen_pi_dEdx_p = (TH1F*) file->Get("dEdx/gen_pi_dEdx_p");
  TH1F *h2_gen_p_dEdx_p  = (TH1F*) file->Get("dEdx/gen_p_dEdx_p");

  h2_gen_K_dEdx_p->SetMarkerColor(kRed);
  h2_gen_pi_dEdx_p->SetMarkerColor(kBlue);
  h2_gen_p_dEdx_p->SetMarkerColor(kGreen);

  h2_gen_K_dEdx_p->SetFillColor(kRed);
  h2_gen_pi_dEdx_p->SetFillColor(kBlue);
  h2_gen_p_dEdx_p->SetFillColor(kGreen);

  pad0->SetLogx();
  h2_gen_K_dEdx_p->SetTitle(";Track momentum [GeV];#frac{dE}{dx}[MeV]");
  h2_gen_K_dEdx_p->GetXaxis()->SetTitleOffset(1.5);
  h2_gen_K_dEdx_p->Draw("box");
  h2_gen_pi_dEdx_p->Draw("box same");
  h2_gen_p_dEdx_p->Draw("box same");

}