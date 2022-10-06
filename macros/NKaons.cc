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

void NKaons()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.root","READ");

  TTree *t_data = (TTree*) file->Get("data");

  Int_t bin  = 30;
  Float_t xmax = 30.0;

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TH1F *h_N_K_Gen = new TH1F("h_N_K_Gen",";N_K_Gen;a.u.",bin,0,xmax);
  TH1F *h_N_K_PFO = new TH1F("h_N_K_PFO",";N_K_PFO;a.u.",bin,0,xmax);
  TH1F *h_N_K_PFO_KK = new TH1F("h_N_K_PFO_KK",";N_K_PFO_KK;a.u.",bin,0,xmax);
  h_N_K_Gen->Sumw2();
  h_N_K_PFO->Sumw2();
  h_N_K_PFO_KK->Sumw2();

  t_data->Draw("N_K_Gen >> h_N_K_Gen");
  t_data->Draw("N_K_PFO >> h_N_K_PFO");
  t_data->Draw("N_K_PFO >> h_N_K_PFO_KK","dEdx_pdg_match == 1");

  Normalize(h_N_K_Gen);
  Normalize(h_N_K_PFO);
  Normalize(h_N_K_PFO_KK);
  
  StyleHist(h_N_K_Gen,kBlack);
  StyleHist(h_N_K_PFO,kGreen+2);
  StyleHist(h_N_K_PFO_KK,kGreen+2);
  h_N_K_PFO_KK->SetLineStyle(2);

  h_N_K_Gen->SetTitle(";N Kaons;a.u.");
  h_N_K_Gen->GetYaxis()->SetRangeUser(0,0.26);

  h_N_K_Gen->Draw("h");
  h_N_K_PFO->Draw("hsame");
  h_N_K_PFO_KK->Draw("hsame");

  TLegend *leg = new TLegend(0.15,0.75,0.45,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h_N_K_Gen,"Generated","l");
  leg->AddEntry(h_N_K_PFO,"Reconstructed","l");
  leg->AddEntry(h_N_K_PFO_KK,"KK Reconstructed Events","l");
  leg->Draw();

  gPad->SetGrid(1,1);

  c0->Draw();

}