#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

// TString LPFO_mode = "Pi";
TString LPFO_mode = "K";
TString ichiral = "eL.pR";
// TString ichiral = "eR.pL";
TString qq = "ss";

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
array<TString,1> processes = {"P2f_z_h"};

unordered_map<pair<TString,TString>,pair<Int_t,Int_t>, hash_pair> production = {
    {{"P2f_z_h", "eL.pR"}, {500010,4994}},
    {{"P2f_z_h", "eR.pL"}, {500012,4994}},
    {{"P4f_ww_h", "eL.pR"}, {500066,4996}},
    {{"P4f_ww_h", "eR.pL"}, {500068,5116}},
    {{"P4f_zz_h", "eL.pR"}, {500062,5052}},
    {{"P4f_zz_h", "eR.pL"}, {500064,5109}},
    {{"Pqqh", "eL.pR"}, {402011,1457}},
    {{"Pqqh", "eR.pL"}, {402012,2278}},
};

Float_t fitRange = 0.6;

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
  if(h->Integral(int_low,int_high)){
    h->Scale( 1.0 / h->Integral(int_low,int_high) );
  }
}

void Normalize2Gen(TH1F *h, TH1F *h_gen)
{
	double intCosReco = h->Integral(20,80);
	double intCosGen  = h_gen->Integral(20,80);
  h->Scale( intCosGen / intCosReco );
}

unordered_map<TString, TH1F*> main_pq(TFile* file, TString prodMode)
{
  gStyle->SetOptStat(0);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get(prodMode + "/gen/h_" + prodMode + "_qcos");
  TH1F *h_reco_LPFO_scos  = (TH1F*) file->Get(prodMode + "/cos/h_" + prodMode + "_" + LPFO_mode + "_scos");
  TH1F *h_reco_LPFO_qcos  = (TH1F*) file->Get(prodMode + "/cos/h_" + prodMode + "_" + LPFO_mode + "_qcos");
  
  // used for pq correction
  TH1F *h_acc_PiPi_cos  = (TH1F*) file->Get(prodMode + "/cos/h_" + prodMode + "_" + LPFO_mode + "_acc_cos");
  TH1F *h_rej_PiPi_cos  = (TH1F*) file->Get(prodMode + "/cos/h_" + prodMode + "_" + LPFO_mode + "_rej_cos");

  // efficiency correction
  Bool_t isEffCorr = true;
  TH1F *h_reco_LPFO_scos_eff_corr;
  TH1F *h_reco_LPFO_qcos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_reco_LPFO_scos_eff = efficiencyCorrection(h_reco_LPFO_scos,LPFO_mode,file,prodMode);
    TH1F* h_reco_LPFO_qcos_eff = efficiencyCorrection(h_reco_LPFO_qcos,LPFO_mode,file,prodMode);
    h_reco_LPFO_scos_eff_corr = resolutionCorrection(h_reco_LPFO_scos_eff,LPFO_mode,file,prodMode);
    h_reco_LPFO_qcos_eff_corr = resolutionCorrection(h_reco_LPFO_qcos_eff,LPFO_mode,file,prodMode);
    // h_reco_LPFO_scos_eff_corr = (TH1F*) h_reco_LPFO_scos_eff->Clone();
    // h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos_eff->Clone();

  }else{
    h_reco_LPFO_scos_eff_corr = (TH1F*) h_reco_LPFO_scos->Clone();
    h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos->Clone();
  }

  TH1F *h_acc_PiPi_cos_eff_corr;
  TH1F *h_rej_PiPi_cos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_acc_PiPi_cos_eff = efficiencyCorrection(h_acc_PiPi_cos,LPFO_mode,file,prodMode);
    TH1F* h_rej_PiPi_cos_eff = efficiencyCorrection(h_rej_PiPi_cos,LPFO_mode,file,prodMode);
    h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos_eff,LPFO_mode,file,prodMode);
    h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos_eff,LPFO_mode,file,prodMode);
    // h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos_eff->Clone();
    // h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos_eff->Clone();

  }else{
    h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
    h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);

  const Int_t nbins = h_reco_LPFO_scos_eff_corr->GetNbinsX();

  // pq correction
  TString pValName = "p_" + LPFO_mode + LPFO_mode + "_" + prodMode;
  TH1F *p_PiPi = new TH1F(pValName,pValName, 50,0,1);
  p_PiPi->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_PiPi_cos_eff_corr, h_rej_PiPi_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_PiPi->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_PiPi->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_LPFO_pq_cos = CorrectHist(prodMode, h_reco_LPFO_qcos_eff_corr, p_vec);
  // TH1F *h_reco_LPFO_pq_cos = (TH1F*) h_reco_LPFO_qcos_eff_corr->Clone();
  StyleHist(h_reco_LPFO_pq_cos,kBlue);

  // output
  unordered_map<TString, TH1F*> hmap;
  hmap["gen"] = h_gen_q_qcos;
  hmap["reco"] = h_reco_LPFO_pq_cos;

  return hmap;


}

void pq_method_PiLPFO_simple()
{
  gStyle->SetOptStat(0);
  try
  {
    unordered_map<TString, unordered_map<TString, TFile*> > file_map;
    for( auto process : processes ){
      for( auto chiral : chirals ){
        Int_t processID = production.at({process,chiral}).first;
        cout << process << " " << chiral << " " << processID << endl;
        TString filename = inputDir + "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I" + processID + "." + process + "." + chiral + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root";
        TFile *file = new TFile(filename,"READ");
        if( !file->IsOpen() ) throw std::runtime_error("File not found");
        file_map[process][chiral] = file;
      }
    }

    unordered_map<TString, TH1F*> h_map = main_pq(file_map["P2f_z_h"][ichiral], qq);
    TH1F* h_gen = (TH1F*) h_map.at("gen")->Clone();
    TH1F* h_reco = (TH1F*) h_map.at("reco")->Clone();

    StyleHist(h_gen,kGreen+2);
    StyleHist(h_reco,kBlue+2);
    Normalize(h_gen);
    Normalize(h_reco);
    h_reco->SetFillStyle(0);

    // Fitting
    TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
    f_gen->SetParNames("S","A");
    h_gen->Fit("f_gen","MNRS");
    f_gen->SetLineColor(kGreen+2);
    cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

    TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
    f_reco->SetParNames("S","A");
    h_reco->Fit("f_reco","MNRS");
    f_reco->SetLineColor(kRed);
    cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;




    // TLegend *legend = new TLegend(0.60,0.75,0.88,0.88);
    TLegend *legend = new TLegend(0.2,0.74,0.5,0.84);
    legend->SetMargin(0.4);
    legend->SetLineColor(0);
    legend->AddEntry(h_gen,"Parton level","f");
    legend->AddEntry(h_reco,"Reconstructed","le");


    TCanvas *c_signal = new TCanvas("c_signal","c_signal",800,800);
    TPad    *p_signal = new TPad("p_signal","p_signal",0,0,1,1);
    StylePad(p_signal,0,0.15,0,0.17);
    h_gen->SetTitle(";cos#theta_{K^{-}};Entries (norm.)");
    h_gen->GetYaxis()->SetRangeUser(0,0.07);
    h_gen->Draw("h");
    h_reco->Draw("e same");
    // f_reco->Draw("same");
    // f_gen->Draw("same");

    legend->Draw("same");


    // Draw polar angle reco/gen ratio
      TCanvas  *c_ratio_genreco = new TCanvas("c_ratio_genreco","c_ratio_genreco",700,900);
      TH1F *rGen  = (TH1F*) h_gen->Clone();
      TH1F *rReco = (TH1F*) h_reco->Clone();
      // rReco->GetYaxis()->SetRangeUser(0,800E3);
      rReco->SetTitle(";cos#theta;Entries");

      auto trp_genreco = new TRatioPlot(rReco,rGen);
      trp_genreco->SetGraphDrawOpt("P");
      trp_genreco->SetSeparationMargin(0.0);
      trp_genreco->Draw();
      trp_genreco->GetLowerRefYaxis()->SetTitle("Data / MC");
      trp_genreco->GetLowerRefGraph()->SetMinimum(0.5);
      trp_genreco->GetLowerRefGraph()->SetMaximum(1.5);
      trp_genreco->GetLowerRefYaxis()->SetLabelSize(0.02);

      trp_genreco->GetUpperPad()->cd();
      TLegend *leg_trp_genreco = new TLegend(0.51,0.74,0.89,0.89);
      leg_trp_genreco->SetMargin(0.4);
      leg_trp_genreco->SetLineColor(0);
      leg_trp_genreco->AddEntry(rReco,"Reconstructed #pi^{-}","l");
      leg_trp_genreco->AddEntry(rGen,"Generated #pi^{-}","l");
      leg_trp_genreco->Draw();


    // stat + sys
    // cout << "Stat + sys error =====\n";
    // for (int ibin = 1; ibin <= h_reco->GetNbinsX(); ibin++)
    // {
    //   cout << h_reco->GetBinError(ibin) / h_reco->GetBinContent(ibin) * 100  << "\n";
    // }



  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}