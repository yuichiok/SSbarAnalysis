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

void gen_K_p_cos()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.root","READ");

  TH2F *h_gen_K_p_cos  = (TH2F*) file->Get("h2_gen_K_p_cos");
  TH2F *h_reco_K_p_cos = (TH2F*) file->Get("h2_reco_K_p_cos");
  h_gen_K_p_cos->Sumw2();
  h_reco_K_p_cos->Sumw2();

  TCanvas *c0 = new TCanvas("c0","c0",1600,800);
  c0->Divide(2,1);
  gPad->SetGrid(1,1);
  c0->cd(1);
  h_gen_K_p_cos->SetTitle("Gen stable K^{#pm} p vs cos#theta");
  h_gen_K_p_cos->Draw("colz");
  c0->cd(2);
  h_reco_K_p_cos->SetTitle("Reco pdg cheat K^{#pm} p vs cos#theta");
  h_reco_K_p_cos->Draw("colz");
  c0->Draw();

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TH2F *h_K_p_cos_div = (TH2F*) h_reco_K_p_cos->Clone();
  h_K_p_cos_div->Divide(h_gen_K_p_cos);
  h_K_p_cos_div->SetTitle("Reco pdg cheat / Gen stable K^{#pm}");
  h_K_p_cos_div->SetMaximum(2.0);
  h_K_p_cos_div->Draw("colz");
  c1->Draw();

  Int_t N_gen_K  = h_gen_K_p_cos->GetEntries();
  Int_t N_reco_K = h_reco_K_p_cos->GetEntries();

  cout << "*************************\n";
  cout << " N_gen_K:  " << N_gen_K  << "\n";
  cout << " N_reco_K: " << N_reco_K << "\n";
  cout << "*************************\n";

}