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

void gen_reco_K_sep_mincos()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.root","READ");
  TFile *file2 = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.p2.root","READ");
  TFile *file5 = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.p5.root","READ");
  TFile *file20 = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.hists.p20.root","READ");

  TH1F *h_gen_reco_K_sep_mincos     = (TH1F*) file->Get("h_gen_reco_K_sep_mincos");
  TH1F *h_gen_reco_K_sep_mincos_p2  = (TH1F*) file2->Get("h_gen_reco_K_sep_mincos");
  TH1F *h_gen_reco_K_sep_mincos_p5  = (TH1F*) file5->Get("h_gen_reco_K_sep_mincos");
  TH1F *h_gen_reco_K_sep_mincos_p20 = (TH1F*) file20->Get("h_gen_reco_K_sep_mincos");

  // Normalize(h_gen_reco_K_sep_mincos);
  // Normalize(h_gen_reco_K_sep_mincos_p2);
  // Normalize(h_gen_reco_K_sep_mincos_p5);
  // Normalize(h_gen_reco_K_sep_mincos_p20);

  StyleHist(h_gen_reco_K_sep_mincos,kBlack);
  StyleHist(h_gen_reco_K_sep_mincos_p2,kRed+2);
  StyleHist(h_gen_reco_K_sep_mincos_p5,kBlue+2);
  StyleHist(h_gen_reco_K_sep_mincos_p20,kGreen+2);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  gPad->SetGrid(1,1);
  gPad->SetLogy();

  h_gen_reco_K_sep_mincos->SetTitle(";Kaon 1 - cos#theta_{gen-reco}; a.u.");
  h_gen_reco_K_sep_mincos->GetXaxis()->SetRangeUser(0,1);
  h_gen_reco_K_sep_mincos->Draw("h");
  h_gen_reco_K_sep_mincos_p2->Draw("hsame");
  h_gen_reco_K_sep_mincos_p5->Draw("hsame");
  h_gen_reco_K_sep_mincos_p20->Draw("hsame");

  c0->Draw();

  TLegend *leg = new TLegend(0.15,0.75,0.45,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h_gen_reco_K_sep_mincos,"No cut","l");
  leg->AddEntry(h_gen_reco_K_sep_mincos_p2,"p > 2.5 GeV","l");
  leg->AddEntry(h_gen_reco_K_sep_mincos_p5,"p > 5.0 GeV","l");
  leg->AddEntry(h_gen_reco_K_sep_mincos_p20,"p > 20.0 GeV","l");
  leg->Draw();

}