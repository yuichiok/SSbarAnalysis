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

void Stability_Purity()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderSize(1);

  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.dd.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  TTree *t_data = (TTree*) file->Get("data");

  Int_t bin  = 30;
  Float_t xmax = 30.0;

  // For differential cosθ
  TH1F *h_gen_N_K_cos  = (TH1F*) file->Get("h_gen_N_K_cos2");
  TH1F *h_reco_N_K_cos = (TH1F*) file->Get("h_reco_N_K_cos2");
  TH1F *h_N_K_corr_cos = (TH1F*) file->Get("h_N_K_corr_cos2");
  StyleHist(h_gen_N_K_cos,kBlack);
  StyleHist(h_reco_N_K_cos,kRed+2);
  StyleHist(h_N_K_corr_cos,kBlue+2);
  
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

  // BinNormalize(h_gen_N_K_cos);
  // BinNormalize(h_reco_N_K_cos);
  // BinNormalize(h_N_K_corr_cos);

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

  TCanvas *c3 = new TCanvas("c3","c3",800,800);

  TH1F * st_pr = (TH1F*) h_reco_N_K_cos->Clone();
  st_pr->Divide(h_gen_N_K_cos);

  st_pr->GetYaxis()->SetRangeUser(0,3);
  st_pr->Draw("h");

  gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.15);
  c3->Draw();



}
