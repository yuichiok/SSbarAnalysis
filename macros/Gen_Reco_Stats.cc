#include <iostream>

void Normalize(TH1F *h)
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

void Gen_Reco_Stats()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.root","READ");

  TTree *t_data = (TTree*) file->Get("data");

  Int_t bin  = 21;
  Float_t xmax = 1.05;

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TH1F *h_stable = new TH1F("h_stable",";stability;a.u.",bin,0,xmax);
  TH1F *h_purity = new TH1F("h_purity",";purity;a.u.",bin,0,xmax);

  TH1F *h_stable_KK = new TH1F("h_stable_KK","h_stable",bin,0,xmax);
  TH1F *h_purity_KK = new TH1F("h_purity_KK","h_purity",bin,0,xmax);

  t_data->Draw("stability >> h_stable", "stability >= 0");
  t_data->Draw("purity >> h_purity", "purity >= 0");

  t_data->Draw("stability >> h_stable_KK", "stability >= 0 && dEdx_pdg_match == 1");
  t_data->Draw("purity >> h_purity_KK", "purity >= 0 && dEdx_pdg_match == 1");


  Normalize(h_stable);
  Normalize(h_purity);
  Normalize(h_stable_KK);
  Normalize(h_purity_KK);
  
  StyleHist(h_stable,kBlack);
  StyleHist(h_purity,kGreen+2);
  h_stable->SetLineStyle(2);
  h_purity->SetLineStyle(2);

  StyleHist(h_stable_KK,kBlack);
  StyleHist(h_purity_KK,kGreen+2);
  h_stable_KK->SetLineStyle(3004);
  h_purity_KK->SetLineStyle(3004);

  h_stable->SetTitle(";Ratio;a.u.");

  h_stable->Draw("h");
  h_purity->Draw("hsame");
  h_stable_KK->Draw("hsame");
  h_purity_KK->Draw("hsame");

  TLegend *leg = new TLegend(0.15,0.75,0.45,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h_stable,"Stability","l");
  leg->AddEntry(h_purity,"Purity","l");
  leg->AddEntry(h_stable_KK,"Stability (KK events)","l");
  leg->AddEntry(h_purity_KK,"Purity (KK events)","l");
  leg->Draw();

  gPad->SetGrid(1,1);

  c0->Draw();

}