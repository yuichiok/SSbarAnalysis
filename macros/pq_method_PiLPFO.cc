#include <iostream>
#include <vector>

#include "include/Styles.hh"
#include "include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector;

void BinNormal(TH1F *h)
{
  const Int_t nbins = h->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++){
    Float_t binc = h->GetBinContent(ibin);
    Float_t binw = h->GetBinWidth(ibin);
    Float_t bin_normal = binc / binw;
    h->SetBinContent(ibin,bin_normal);
  }

}

void Normalize(TH1F *h)
{
  // h->Scale( 1.0 / h->GetEntries() );
  const Int_t nbins = h->GetNbinsX();
  Int_t nbins4 = nbins / 4;
  Int_t int_high = (nbins / 2) + nbins4;
  Int_t int_low  = (nbins / 2 + 1) - nbins4;
  h->Scale( 1.0 / h->Integral(int_low,int_high) );
}

void Normalize2Gen(TH1F *h, TH1F *h_gen)
{
	double intCosReco = h->Integral(20,80);
	double intCosGen  = h_gen->Integral(20,80);
  h->Scale( intCosGen / intCosReco );
}


void main_pq()
{
  gStyle->SetOptStat(0);

  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.dd.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.dd.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.test.root","READ");

  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  if (!file->IsOpen()) return;

  // reco and gen polar
  TH1F *h_gen_q_qcos  = (TH1F*) file->Get("h_gen_q_qcos");
  TH1F *h_reco_K_scos = (TH1F*) file->Get("h_reco_Pi_scos");
  TH1F *h_reco_K_qcos = (TH1F*) file->Get("h_reco_Pi_qcos");
  TH1F *h_cheat_Pi_qcos = (TH1F*) file->Get("h_cheat_Pi_qcos");

  cout << "========================" << endl;
  cout << "NReco Pi 30th bin = " << h_reco_K_qcos->GetBinContent(30) << endl;

  // efficiency correction
  Bool_t isEffCorr = false;
  TH1F *h_reco_K_scos_eff_corr;
  TH1F *h_reco_K_qcos_eff_corr;
  if (isEffCorr)
  {
    h_reco_K_scos_eff_corr = Efficiency_Correction(h_reco_K_scos,"scos_corr",file);
    h_reco_K_qcos_eff_corr = Efficiency_Correction(h_reco_K_qcos,"qcos_corr",file);
  }else{
    h_reco_K_scos_eff_corr = (TH1F*) h_reco_K_scos->Clone();
    h_reco_K_qcos_eff_corr = (TH1F*) h_reco_K_qcos->Clone();
  }

  // used for pq correction
  TH1F *h_acc_KK_cos  = (TH1F*) file->Get("pq/h_acc_PiPi_cos");
  TH1F *h_rej_KK_cos  = (TH1F*) file->Get("pq/h_rej_PiPi_cos");

  TH1F *h_acc_KK_cos_eff_corr;
  TH1F *h_rej_KK_cos_eff_corr;
  if (isEffCorr)
  {
    h_acc_KK_cos_eff_corr = Efficiency_Correction(h_acc_KK_cos,"acc_corr",file);
    h_rej_KK_cos_eff_corr = Efficiency_Correction(h_rej_KK_cos,"rej_corr",file);
  }else{
    h_acc_KK_cos_eff_corr = (TH1F*) h_acc_KK_cos->Clone();
    h_rej_KK_cos_eff_corr = (TH1F*) h_rej_KK_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);
  h_gen_q_qcos->SetFillStyle(0);
  h_gen_q_qcos->SetLineStyle(2);

  StyleHist(h_cheat_Pi_qcos,kOrange+1);
  h_cheat_Pi_qcos->SetFillStyle(0);
  h_cheat_Pi_qcos->SetLineStyle(2);

  StyleHist(h_reco_K_scos_eff_corr,kBlack);
  h_reco_K_scos_eff_corr->SetFillStyle(0);
  StyleHist(h_reco_K_qcos_eff_corr,kRed+2);
  StyleHist(h_acc_KK_cos_eff_corr,kRed+2);
  StyleHist(h_rej_KK_cos_eff_corr,kBlue+2);

  const Int_t nbins = h_reco_K_scos_eff_corr->GetNbinsX();

  TH1F *p_KK = new TH1F("p_KK", "p_KK", 50,0,1);
  p_KK->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_KK_cos_eff_corr, h_rej_KK_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_KK->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_KK->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_K_pq_cos = CorrectHist(h_reco_K_qcos_eff_corr, p_vec);
  StyleHist(h_reco_K_pq_cos,kBlue);

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TPad *pad0 = new TPad("pad0", "pad0",0,0,1,1);
  StylePad(pad0,0,0.12,0,0.15);

  Normalize2Gen(h_gen_q_qcos,h_reco_K_scos_eff_corr);
  Normalize2Gen(h_cheat_Pi_qcos,h_gen_q_qcos);
  // Normalize(h_reco_K_scos_eff_corr);
  // Normalize(h_reco_K_pq_cos);
  // Normalize(h_reco_K_qcos_eff_corr);

  // Fitting
  TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_reco->SetParNames("S","A");
  h_reco_K_pq_cos->Fit("f_reco","MNRS");

  // h_reco_K_pq_cos->GetYaxis()->SetRangeUser(0,50E3);
  h_reco_K_pq_cos->SetTitle(";cos#theta_{#pi^{-}};a.u.");
  h_reco_K_pq_cos->Draw("h");
  h_reco_K_qcos_eff_corr->Draw("hsame");
  h_reco_K_scos_eff_corr->Draw("hsame");
  h_gen_q_qcos->Draw("hsame");
  h_cheat_Pi_qcos->Draw("hsame");

  f_reco->Draw("same");

  TLegend *leg = new TLegend(0.2,0.76,0.7,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h_gen_q_qcos,"Generated quark angle","l");
  leg->AddEntry(h_cheat_Pi_qcos,"Cheated #pi^{-} PFO","l");
  leg->AddEntry(h_reco_K_scos_eff_corr,"Reconstructed #pi^{-} matched with quark angle","l");
  leg->AddEntry(h_reco_K_qcos_eff_corr,"Reconstructed #pi^{-}","l");
  leg->AddEntry(h_reco_K_pq_cos,"Reconstructed #pi^{-} (corrected)","l");
  leg->Draw();

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TPad *pad1 = new TPad("pad1", "pad1",0,0,1,1);
  StylePad(pad1,0,0.12,0,0.15);
  
  StyleHist(p_KK,kGreen+2);
  p_KK->SetTitle(";cos#theta_{#pi^{-}};p value");
  p_KK->GetYaxis()->SetRangeUser(0,1);
  p_KK->Draw("h");

  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  TGaxis::SetMaxDigits(3);
  gPad->SetGrid(1,1);
  h_acc_KK_cos_eff_corr->SetTitle(";cos#theta_{#pi^{-}};Entries");

  // h_acc_KK_cos_eff_corr->GetYaxis()->SetRangeUser(0,50E3);
  h_acc_KK_cos_eff_corr->Draw("h");
  h_rej_KK_cos_eff_corr->Draw("hsame");

  TH1F * acc_full = (TH1F*) h_acc_KK_cos_eff_corr->Clone();
  TH1F * acc_add  = new TH1F("acc_add","acc_add",nbins,-1,1);

  for (int i = 1; i < nbins / 2 + 1; i++)
  {
    float accepted = acc_full->GetBinContent(nbins + 1 - i);
    accepted += acc_full->GetBinContent(i);
    acc_add->SetBinContent(i,accepted);
    acc_add->SetBinContent(nbins+1-i,accepted);
  }
  StyleHist(acc_add,kRed);
  acc_add->Draw("hsame");

  TLegend *leg2 = new TLegend(0.15,0.75,0.45,0.85);
  leg2->SetLineColor(0);
  leg2->AddEntry(h_acc_KK_cos_eff_corr,"N Accepted","l");
  leg2->AddEntry(h_rej_KK_cos_eff_corr,"N Rejected","l");
  leg2->AddEntry(acc_add,"N Accepted + opp. bin","l");
  leg2->Draw();

}

void pq_method_PiLPFO()
{
  try
  {
    main_pq();
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}