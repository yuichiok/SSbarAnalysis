#include <iostream>

void dEdx()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp10_pNaN.tpc0.hists.all.root","READ");

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  gPad->SetGrid(1,1);

  TH1F *h2_gen_K_dEdx_p  = (TH1F*) file->Get("dEdx/gen_K_dEdx_p");
  TH1F *h2_gen_pi_dEdx_p = (TH1F*) file->Get("dEdx/gen_pi_dEdx_p");
  TH1F *h2_gen_p_dEdx_p  = (TH1F*) file->Get("dEdx/gen_p_dEdx_p");


}