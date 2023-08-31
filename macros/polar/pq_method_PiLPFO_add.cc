#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString prod_mode = "ud";
TString LPFO_mode = "Pi";

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
	double intCosReco = h->Integral(20,30);
	double intCosGen  = h_gen->Integral(20,30);
  h->Scale( intCosGen / intCosReco );
}

unordered_map<TString, TH1F*> main_pq(TFile* file, Float_t &ratio)
{
  gStyle->SetOptStat(0);

  TString prodMode = getProductionMode(file);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get("gen/h_qcos");
  TH1F *h_reco_Pi_scos  = (TH1F*) file->Get("cos/h_Pi_scos");
  TH1F *h_reco_Pi_qcos  = (TH1F*) file->Get("cos/h_Pi_qcos");

  cout << prodMode << " reco entry = " << h_reco_Pi_qcos->GetEntries() << endl;
  cout << prodMode << " gen entry  = " << h_gen_q_qcos->GetEntries() << endl;
  cout << prodMode << " reco eff   = " << (float)h_reco_Pi_qcos->GetEntries() / (float)h_gen_q_qcos->GetEntries() << endl;
  
  ratio = (Float_t) h_reco_Pi_qcos->GetEntries() / (Float_t) h_gen_q_qcos->GetEntries();

  // used for pq correction
  TH1F *h_acc_PiPi_cos  = (TH1F*) file->Get("cos/h_Pi_acc_cos");
  TH1F *h_rej_PiPi_cos  = (TH1F*) file->Get("cos/h_Pi_rej_cos");

  // efficiency correction
  Bool_t isEffCorr = true;
  TH1F *h_reco_Pi_scos_eff_corr;
  TH1F *h_reco_Pi_qcos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_reco_Pi_scos_eff = efficiencyCorrection(h_reco_Pi_scos,"Pi",file);
    TH1F* h_reco_Pi_qcos_eff = efficiencyCorrection(h_reco_Pi_qcos,"Pi",file);
    h_reco_Pi_scos_eff_corr = resolutionCorrection(h_reco_Pi_scos_eff,"Pi",file);
    h_reco_Pi_qcos_eff_corr = resolutionCorrection(h_reco_Pi_qcos_eff,"Pi",file);

  }else{
    h_reco_Pi_scos_eff_corr = (TH1F*) h_reco_Pi_scos->Clone();
    h_reco_Pi_qcos_eff_corr = (TH1F*) h_reco_Pi_qcos->Clone();
  }

  TH1F *h_acc_PiPi_cos_eff_corr;
  TH1F *h_rej_PiPi_cos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_acc_PiPi_cos_eff = efficiencyCorrection(h_acc_PiPi_cos,"Pi",file);
    TH1F* h_rej_PiPi_cos_eff = efficiencyCorrection(h_rej_PiPi_cos,"Pi",file);
    h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos_eff,"Pi",file);
    h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos_eff,"Pi",file);
  }else{
    h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
    h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);

  const Int_t nbins = h_reco_Pi_scos_eff_corr->GetNbinsX();

  // pq correction
  TString pValName = "p_PiPi_" + prodMode;
  TH1F *p_PiPi = new TH1F(pValName,pValName, 50,0,1);
  p_PiPi->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_PiPi_cos_eff_corr, h_rej_PiPi_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_PiPi->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_PiPi->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_Pi_pq_cos = CorrectHist(prodMode, h_reco_Pi_qcos_eff_corr, p_vec);
  StyleHist(h_reco_Pi_pq_cos,kBlue);

  // Fitting
  TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_gen->SetParNames("S","A");
  h_gen_q_qcos->Fit("f_gen","MNRS");
  cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

  TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_reco->SetParNames("S","A");
  h_reco_Pi_pq_cos->Fit("f_reco","MNRS");
  cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;

  // output
  unordered_map<TString, TH1F*> hmap;
  hmap["gen"] = h_gen_q_qcos;
  hmap["reco"] = h_reco_Pi_pq_cos;

  return hmap;


}

void pq_method_PiLPFO_add()
{
  try
  {
    TFile *uuFile = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root","READ");
    TFile *ddFile = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.dd.KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root","READ");
    if (!uuFile->IsOpen() || !ddFile->IsOpen()) return;

    Float_t dd_ratio;
    Float_t uu_ratio;

    unordered_map<TString, TH1F*> hmap_uu = main_pq(uuFile, uu_ratio);
    unordered_map<TString, TH1F*> hmap_dd = main_pq(ddFile, dd_ratio);

    cout << "uu ratio = " << uu_ratio << endl;
    cout << "dd ratio = " << dd_ratio << endl;
    cout << "weight   = " << dd_ratio / uu_ratio << endl;

    hmap_dd.at("gen")->Scale( dd_ratio / uu_ratio );

    TCanvas *cts_reco = new TCanvas("cts_reco","cts_reco",800,800);
    TCanvas *cts_gen  = new TCanvas("cts_gen ","cts_gen ",800,800);

    THStack *ts_reco = new THStack("ts_reco","ts_reco");
    ts_reco->Add(hmap_uu.at("reco"));
    ts_reco->Add(hmap_dd.at("reco"));

    THStack *ts_gen = new THStack("ts_gen","ts_gen");
    ts_gen->Add(hmap_uu.at("gen"));
    ts_gen->Add(hmap_dd.at("gen"));

    cts_reco->cd();
    ts_reco->Draw("h");

    cts_gen->cd();
    ts_gen->Draw("h");

    TH1F *h_reco = (TH1F*) hmap_uu.at("reco")->Clone();
    h_reco->Add(hmap_dd.at("reco"));
    StyleHist(h_reco,kBlue);
    
    TH1F *h_gen = (TH1F*) hmap_uu.at("gen")->Clone();
    h_gen->Add(hmap_dd.at("gen"));
    StyleHist(h_gen,kGreen+1);

    TCanvas *cTotal = new TCanvas("cTotal","cTotal",800,800);
    Normalize2Gen(h_gen,h_reco);
    h_gen->Draw("h");
    h_reco->Draw("h same");




  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}