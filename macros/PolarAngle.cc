#include <iostream>

void Normalize(TH1F *h)
{
  // h->Scale( 1.0 / h->GetEntries() );
  // h->Scale( 1.0 / h->Integral(30,70) );
  h->Scale( 1.0 / h->Integral(12,88) );
}

void StyleHist(TH1F *h, Color_t col)
{
  h->SetLineWidth(3);
  h->SetLineColor(col);
  h->SetFillStyle(3002);
  h->SetFillColor(col);
}

void PolarAngle()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.root","READ");

  TH1F *h_gen_q_qcos  = (TH1F*) file->Get("h_gen_q_qcos");
  TH1F *h_reco_K_qcos = (TH1F*) file->Get("h_reco_K_qcos");
  h_gen_q_qcos->Sumw2();
  h_reco_K_qcos->Sumw2();

  Normalize(h_gen_q_qcos);
  Normalize(h_reco_K_qcos);

  StyleHist(h_gen_q_qcos,kBlack);
  StyleHist(h_reco_K_qcos,kBlue);


  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  gPad->SetGrid(1,1);

  h_gen_q_qcos->Draw("h");
  h_reco_K_qcos->Draw("hsame");
  h_reco_K_qcos->SetTitle(";Kaon cos#theta (no filter); Entries / 0.02");
  c0->Draw();

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gPad->SetGrid(1,1);

  TH1F *eff = (TH1F*) h_reco_K_qcos->Clone();
  eff->Sumw2();
  StyleHist(eff,kGreen+1);
  eff->Divide(h_gen_q_qcos);
  eff->Draw("h");
  eff->GetYaxis()->SetRangeUser(0,1.6);
  eff->SetTitle(";Kaon cos#theta (no filter); Efficiency / 0.02");

  c1->Draw();

  c0->Update();
  c0->Modified();
  c1->Update();


  // file->Close();

}