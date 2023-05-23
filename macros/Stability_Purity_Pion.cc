#include "Styles.cc"

vector<TH1F*> extract_ps_hists(TFile *file)
{
  TH1F *h_gen_N_K_cos  = (TH1F*) file->Get("h_gen_N_K_cos");
  TH1F *h_reco_N_K_cos = (TH1F*) file->Get("h_reco_N_K_cos");
  TH1F *h_N_K_corr_cos = (TH1F*) file->Get("h_N_K_corr_cos");

  TH1F *h_stable_cos = (TH1F*) h_N_K_corr_cos->Clone();
  TH1F *h_purity_cos = (TH1F*) h_N_K_corr_cos->Clone();
  h_stable_cos->Divide(h_gen_N_K_cos);
  h_purity_cos->Divide(h_reco_N_K_cos);

  vector<TH1F*> hists;
  hists.push_back(h_stable_cos);
  hists.push_back(h_purity_cos);

  return hists;

}

void Stability_Purity(TFile *file)
{
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderSize(1);

  TCanvas *c_stability_purity = new TCanvas("c_stability_purity","c_stability_purity",800,800);
  TPad *pad_stability_purity  = new TPad("pad_stability_purity", "pad_stability_purity",0,0,1,1);
  StylePad(pad_stability_purity,0,0.12,0,0.15);

  enum {kStability,kPurity};
  vector<TH1F*> h_ps = extract_ps_hists(file);

  StyleHist(h_ps.at(kStability),kBlue+2);
  StyleHist(h_ps.at(kPurity),kRed+2);

  h_ps.at(kStability)->Draw("");
  h_ps.at(kPurity)->Draw("same");

  // TLegend *leg2 = new TLegend(0.2,0.15,0.35,0.23);
  // leg2->SetLineColor(0);
  // leg2->AddEntry(h_ss_weight,"u#bar{u}","l");
  // leg2->AddEntry(h_uu_weight,"s#bar{s}","l");
  // leg2->Draw();

}



void Stability_Purity_Pion()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderSize(1);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  Stability_Purity(file);

}