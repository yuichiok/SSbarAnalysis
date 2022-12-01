#include <iostream>

void Normalize(TH1F *h)
{
  h->Scale( 1.0 / h->GetEntries() );
}

void Normalize_Integral(TH1F *h)
{
  h->Scale( 1.0 / h->Integral(1,100) );
}

void BinNormalize(TH1F *h)
{
  Int_t nbins = h->GetNbinsX();

  for ( int ibin=1; ibin <= nbins; ibin++ ){
    Float_t binc = h->GetBinContent(ibin);
    Float_t binw = h->GetBinWidth(ibin);
    Float_t bin_div = binc / binw;
    h->SetBinContent(ibin,bin_div);
  }

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
  gStyle->SetPadBorderSize(1);

  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.p5.root","READ");
  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp15_pNaN.tpc0.hists.all.root","READ");

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

  BinNormalize(h_N_K_Gen);
  BinNormalize(h_N_K_PFO);
  BinNormalize(h_N_K_PFO_KK);

  Normalize(h_N_K_Gen);
  Normalize(h_N_K_PFO);
  Normalize(h_N_K_PFO_KK);
  
  StyleHist(h_N_K_Gen,kBlack);
  StyleHist(h_N_K_PFO,kGreen+2);
  StyleHist(h_N_K_PFO_KK,kGreen+2);
  h_N_K_PFO_KK->SetLineStyle(2);

  h_N_K_Gen->SetTitle(";N Kaons;a.u.");
  h_N_K_Gen->GetXaxis()->SetRangeUser(0,25);
  h_N_K_Gen->GetYaxis()->SetRangeUser(0,1);

  h_N_K_Gen->Draw("h");
  h_N_K_PFO->Draw("hsame");
  // h_N_K_PFO_KK->Draw("hsame");

  TLegend *leg0 = new TLegend(0.6,0.75,0.85,0.85);
  leg0->SetLineColor(0);
  leg0->AddEntry(h_N_K_Gen,"Generated","l");
  leg0->AddEntry(h_N_K_PFO,"Reconstructed","l");
  // leg0->AddEntry(h_N_K_PFO_KK,"KK Reconstructed Events","l");
  leg0->Draw();

  gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.15);
  c0->Draw();


  // For differential cosθ
  TH1F *h_gen_N_K_cos  = (TH1F*) file->Get("h_gen_N_K_cos");
  TH1F *h_reco_N_K_cos = (TH1F*) file->Get("h_reco_N_K_cos");
  TH1F *h_N_K_corr_cos = (TH1F*) file->Get("h_N_K_corr_cos");
  StyleHist(h_gen_N_K_cos,kBlack);
  StyleHist(h_reco_N_K_cos,kRed+2);
  StyleHist(h_N_K_corr_cos,kBlue+2);
  
  h_gen_N_K_cos->Sumw2();
  h_reco_N_K_cos->Sumw2();
  h_N_K_corr_cos->Sumw2();

  TH1F *h_stable_cos = (TH1F*) h_N_K_corr_cos->Clone();
  h_stable_cos->Divide(h_gen_N_K_cos);
  StyleHist(h_stable_cos,kBlack);

  TH1F *h_purity_cos = (TH1F*) h_N_K_corr_cos->Clone();
  h_purity_cos->Divide(h_reco_N_K_cos);
  StyleHist(h_purity_cos,kGreen+2);


  TCanvas *c1 = new TCanvas("c1","c1",800,800);

  // BinNormalize(h_stable_cos);
  // BinNormalize(h_purity_cos);

  h_stable_cos->SetTitle(";cos#theta;Ratio");
  h_stable_cos->GetYaxis()->SetRangeUser(0,1);
  h_stable_cos->Draw("h");
  h_purity_cos->Draw("hsame");

  TLegend *leg1 = new TLegend(0.3,0.78,0.45,0.85);
  leg1->SetLineColor(0);
  leg1->AddEntry(h_stable_cos,"Stability","l");
  leg1->AddEntry(h_purity_cos,"Purity","l");
  leg1->Draw();

  gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.15);
  c1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",800,800);

  BinNormalize(h_gen_N_K_cos);
  BinNormalize(h_reco_N_K_cos);
  BinNormalize(h_N_K_corr_cos);

  h_reco_N_K_cos->SetTitle(";cos#theta;N Kaons (a.u.)");
  // h_gen_N_K_cos->GetYaxis()->SetRangeUser(1E3,2E5);
  h_reco_N_K_cos->Draw("h");
  h_gen_N_K_cos->Draw("hsame");
  h_N_K_corr_cos->Draw("hsame");

  TLegend *leg2 = new TLegend(0.3,0.75,0.6,0.85);
  leg2->SetLineColor(0);
  leg2->AddEntry(h_gen_N_K_cos,"N Generated Kaons","l");
  leg2->AddEntry(h_reco_N_K_cos,"N Reconstructed Kaons","l");
  leg2->AddEntry(h_N_K_corr_cos,"N Correct Reco Kaons","l");
  leg2->Draw();

  gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.15);
  c2->SetLogy();
  c2->Draw();

}
