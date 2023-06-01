#include <iostream>
#include "include/Styles.hh"

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

void plot_stability_purity(TFile *file)
{
  Int_t bin  = 30;
  Float_t xmax = 30.0;

  // For differential cosθ
  TH1F *h_gen_N_K_cos  = (TH1F*) file->Get("h_gen_N_K_cos");
  TH1F *h_reco_N_K_cos = (TH1F*) file->Get("h_reco_N_K_cos");
  TH1F *h_N_K_corr_cos = (TH1F*) file->Get("h_N_K_corr_cos");
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

TH1F * extract_hists(TFile *file)
{
  TH1F *h_gen_N_K_cos  = (TH1F*) file->Get("h_gen_N_K_cos");
  TH1F *h_reco_N_K_cos = (TH1F*) file->Get("h_reco_N_K_cos");
  TH1F *h_N_K_corr_cos = (TH1F*) file->Get("h_N_K_corr_cos");

  TH1F *h_stable_cos = (TH1F*) h_N_K_corr_cos->Clone();
  TH1F *h_purity_cos = (TH1F*) h_N_K_corr_cos->Clone();
  h_stable_cos->Divide(h_gen_N_K_cos);
  h_purity_cos->Divide(h_reco_N_K_cos);

  TH1F *h_weight = (TH1F*) h_purity_cos->Clone();
  h_weight->Divide(h_stable_cos);

  h_weight->SetTitle(";cos#theta;P/S ratio");
  h_weight->GetYaxis()->SetRangeUser(0,1.5);

  return h_weight;

}

void plot_sp_ratio(TFile *uu_file, TFile *ss_file)
{
  TCanvas *c_plot_sp_ratio = new TCanvas("c_plot_sp_ratio","c_plot_sp_ratio",800,800);
  TPad *pad_sp_ratio = new TPad("pad_sp_ratio", "pad_sp_ratio",0,0,1,1);
  StylePad(pad_sp_ratio,0,0.12,0,0.15);

  TH1F *h_uu_weight = extract_hists(uu_file);
  TH1F *h_ss_weight = extract_hists(ss_file);

  StyleHist(h_uu_weight,kBlue+2);
  StyleHist(h_ss_weight,kRed+2);

  h_ss_weight->Draw("");
  h_uu_weight->Draw("same");

  TLegend *leg2 = new TLegend(0.2,0.15,0.35,0.23);
  leg2->SetLineColor(0);
  leg2->AddEntry(h_ss_weight,"u#bar{u}","l");
  leg2->AddEntry(h_uu_weight,"s#bar{s}","l");
  leg2->Draw();

}

void Stability_Purity()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderSize(1);

  // TFile *uu_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  TFile *uu_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *ss_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  TFile *ss_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *us_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.dd.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  TFile *ud_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.testfinal.root","READ");

  plot_stability_purity(ud_file);
  // plot_stability_purity(ss_file);
  // plot_sp_ratio(ss_file,uu_file);


}
