#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString LPFO_mode = "K";
// TString chiral = "eL.pR";
TString chiral = "eR.pL";

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
// array<TString,4> processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"};
array<TString,1> processes = {"P2f_z_h"};
// array<TString,6> qqbars    = {"dd", "uu", "ss", "cc", "bb", "rr"};
// array<TString,3> qqbars    = {"dd", "uu", "ss"};
array<TString,1> qqbars    = {"ss"};

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

Float_t fitRange = 0.8;

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

    // h_reco_LPFO_scos_eff_corr = efficiencyCorrection(h_reco_LPFO_scos,LPFO_mode,file,prodMode);
    // h_reco_LPFO_qcos_eff_corr = efficiencyCorrection(h_reco_LPFO_qcos,LPFO_mode,file,prodMode);

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

    // h_acc_PiPi_cos_eff_corr = efficiencyCorrection(h_acc_PiPi_cos,LPFO_mode,file,prodMode);
    // h_rej_PiPi_cos_eff_corr = efficiencyCorrection(h_rej_PiPi_cos,LPFO_mode,file,prodMode);

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

void pq_method_KLPFO_add()
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

    // TH1F *h_gen  = (TH1F*) file_map["P2f_z_h"][chiral]->Get("ss/gen/h_ss_qcos");
    // TH1F *h_reco = (TH1F*) file_map["P2f_z_h"][chiral]->Get("ss/cos/h_ss_" + LPFO_mode + "_qcos");
    // StyleHist(h_gen,kRed+2);
    // StyleHist(h_reco,kRed+2);
    // Normalize(h_gen);
    // Normalize(h_reco);
    // h_reco->SetFillStyle(0);

    // TLegend *legend = new TLegend(0.60,0.75,0.88,0.88);
    // legend->AddEntry(h_gen,"Parton level","f");
    // legend->AddEntry(h_reco,"Reconstructed","le");


    // TCanvas *c_signal = new TCanvas("c_signal","c_signal",800,800);
    // TPad    *p_signal = new TPad("p_signal","p_signal",0,0,1,1);
    // StylePad(p_signal,0.1,0.1,0.08,0.155);
    // // p_signal->Draw();
    // h_gen->SetTitle(";cos#theta;Entries");
    // h_gen->Draw("h");
    // h_reco->Draw("e same");
    // legend->Draw("same");




    unordered_map<TString, THStack*> hs_map;
    hs_map["gen"]  = new THStack("hs_gen",";cos#theta;Entries");
    hs_map["reco"] = new THStack("hs_reco",";cos#theta;Entries");
    
    TLegend *legend = new TLegend(0.60,0.75,0.88,0.88);

    // for ( auto process : processes ){
    //   if(process == "P2f_z_h"){
    //     for ( auto qq : qqbars ){
    //       unordered_map<TString, TH1F*> reco_gen_map = main_pq(file_map[process][chiral], qq);
    //       // TH1F *h_gen  = (TH1F*) file_map[process][chiral]->Get(qq + "/gen/h_" + qq + "_qcos");
    //       // TH1F *h_reco = (TH1F*) file_map[process][chiral]->Get(qq + "/cos/h_" + qq + "_" + LPFO_mode + "_qcos");
    //       TH1F * h_gen = reco_gen_map["gen"];
    //       TH1F * h_reco = reco_gen_map["reco"];
    //       Normalize(h_gen);
    //       Normalize(h_reco);
    //       h_gen->SetLineWidth(3);
    //       h_gen->SetFillStyle(3002);
    //       h_reco->SetLineWidth(3);
          
    //       legend->AddEntry(h_reco,qq,"l");

    //       hs_map.at("gen")->Add(h_gen);
    //       hs_map.at("reco")->Add(h_reco);
    //     }

    //   }
    // }

    // gStyle->SetPalette(kRainbow);
    // TCanvas * c_polar_nostack = new TCanvas("c_polar_nostack","c_polar_nostack",800,800);
    // hs_map.at("reco")->Draw("he plc nostack");
    // hs_map.at("gen")->Draw("he plc nostack same");
    // legend->Draw("same");
    // TCanvas * c_polar_stack = new TCanvas("c_polar_stack","c_polar_stack",800,800);
    // hs_map.at("gen")->Draw("he plc");
    // hs_map.at("reco")->Draw("he plc same");


    // Fitting
    unordered_map<TString, TH1F*> reco_gen_map_fit = main_pq(file_map["P2f_z_h"]["eR.pL"], "ss");
    // unordered_map<TString, TH1F*> reco_gen_map_fit = main_pq(file_map["P2f_z_h"]["eR.pL"], "ss");
    TH1F *h_gen_q_qcos     = (TH1F*) reco_gen_map_fit.at("gen")->Clone();
    TH1F *h_reco_Pi_pq_cos = (TH1F*) reco_gen_map_fit.at("reco")->Clone();
    Normalize(h_gen_q_qcos);
    Normalize(h_reco_Pi_pq_cos);

    TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-0.8,0.8);
    f_gen->SetParNames("S","A");
    h_gen_q_qcos->Fit("f_gen","MNRS");
    cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

    TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-0.8,0.8);
    f_reco->SetParNames("S","A");
    h_reco_Pi_pq_cos->Fit("f_reco","MNRS");
    cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;
    
    // output AFB
    Float_t AFB_gen  = AFB_calculation(f_gen);
    Float_t AFB_reco = AFB_calculation(f_reco);
    cout << "Gen  AFB = " << AFB_gen << endl;
    cout << "Reco AFB = " << AFB_reco << endl;


    // Draw polar angle
    TCanvas *c_polar = new TCanvas("c_polar","c_polar",800,800);
    TPad  *pad_polar = new TPad("pad_polar", "pad_polar",0,0,1,1);
    StylePad(pad_polar,0,0.12,0,0.15);

    h_reco_Pi_pq_cos->SetTitle(";cos#theta;Entries (norm.)");
    h_reco_Pi_pq_cos->GetYaxis()->SetRangeUser(0,0.07);
    h_reco_Pi_pq_cos->Draw("");
    h_gen_q_qcos->Draw("hsame");

    f_gen->SetLineWidth(4);
    f_gen->SetLineColor(kRed);
    f_reco->SetLineWidth(4);
    f_reco->SetLineColor(kBlack);
    f_gen->Draw("same");
    f_reco->Draw("same");

    TLegend *leg = new TLegend(0.59,0.65,0.89,0.85);
    leg->SetMargin(0.4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_gen_q_qcos,"Parton level","f");
    leg->AddEntry(h_reco_Pi_pq_cos,"Data","le");
    leg->AddEntry(f_gen,"Fit parton level","l");
    leg->AddEntry(f_reco,"Fit data","l");
    leg->Draw("same");


  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}