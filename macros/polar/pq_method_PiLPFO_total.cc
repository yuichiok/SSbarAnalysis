#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::array; using std::unordered_map;

// TString prod_mode = "uu";
// TString chiral    = "eR.pL";
TString LPFO_mode = "K";
Float_t TopRange = 300E3;

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
array<TString,4> processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"};
array<TString,6> qqbars    = {"dd", "uu", "ss", "cc", "bb", "rr"};

unordered_map<pair<TString,TString>,Int_t> production = {
    {{"P2f_z_h", "eL.pR"}, 500010},
    {{"P2f_z_h", "eR.pL"}, 500012},
    {{"P4f_ww_h", "eL.pR"}, 500066},
    {{"P4f_ww_h", "eR.pL"}, 500068},
    {{"P4f_zz_h", "eL.pR"}, 500062},
    {{"P4f_zz_h", "eR.pL"}, 500064},
    {{"Pe1e1h", "eL.pL"}, 402013},
    {{"Pe1e1h", "eL.pR"}, 402001},
    {{"Pe1e1h", "eR.pL"}, 402002},
    {{"Pe1e1h", "eR.pR"}, 402014}
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
  h->Scale( 1.0 / h->Integral(int_low,int_high) );
}

void Normalize2Gen(TH1F *h, TH1F *h_gen)
{
	double intCosReco = h->Integral(20,80);
	double intCosGen  = h_gen->Integral(20,80);
  h->Scale( intCosGen / intCosReco );
}

TString histLabel(TString process, TString category){

  if (category.Length() != 2) {
    std::cerr << "Error: category must be two letters" << std::endl;
    return "";
  }

  if(process=="P2f_z_h" && category!="rr"){
    std::string category_str(category.Data(), category.Length());
    return TString(category_str.substr(0, 1) + "#bar{" + category_str.substr(1, 1) + "}").Data();
  }else if (category == "rr"){
    return "Rad. Ret.";
  }else if (process=="P4f_ww_h"){
    return "WW";
  }else if (process=="P4f_zz_h"){
    return "ZZ";
  }else if (process=="Pe1e1h"){
    return "q#bar{q}H";
  }
  return "";

}

unordered_map<TString, TH1F*> main_pq(TFile* file, TString category)
{
  gStyle->SetOptStat(0);

  TString prodMode = getProductionMode(file);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get(category + "/gen/h_" + category + "_qcos");
  TH1F *h_reco_LPFO_scos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_scos");
  TH1F *h_reco_LPFO_qcos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_qcos");

  cout << prodMode << " reco entry = " << h_reco_LPFO_qcos->GetEntries() << endl;
  cout << prodMode << " gen entry  = " << h_gen_q_qcos->GetEntries() << endl;
  cout << prodMode << " reco eff   = " << (float)h_reco_LPFO_qcos->GetEntries() / (float)h_gen_q_qcos->GetEntries() << endl;
  
  // used for pq correction
  TH1F *h_acc_PiPi_cos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_acc_cos");
  TH1F *h_rej_PiPi_cos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_rej_cos");

  // efficiency correction
  Bool_t isEffCorr = true;
  TH1F *h_reco_LPFO_scos_eff_corr;
  TH1F *h_reco_LPFO_qcos_eff_corr;
  // if (isEffCorr && category!="bg")
  if (isEffCorr)
  {
    TH1F* h_reco_LPFO_scos_eff = efficiencyCorrection(h_reco_LPFO_scos,LPFO_mode,file,category);
    TH1F* h_reco_LPFO_qcos_eff = efficiencyCorrection(h_reco_LPFO_qcos,LPFO_mode,file,category);
    h_reco_LPFO_scos_eff_corr = resolutionCorrection(h_reco_LPFO_scos_eff,LPFO_mode,file,category);
    h_reco_LPFO_qcos_eff_corr = resolutionCorrection(h_reco_LPFO_qcos_eff,LPFO_mode,file,category);

  }else{
    h_reco_LPFO_scos_eff_corr = (TH1F*) h_reco_LPFO_scos->Clone();
    h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos->Clone();
  }

  TH1F *h_acc_PiPi_cos_eff_corr;
  TH1F *h_rej_PiPi_cos_eff_corr;
  // if (isEffCorr && category!="bg")
  if (isEffCorr)
  {
    TH1F* h_acc_PiPi_cos_eff = efficiencyCorrection(h_acc_PiPi_cos,LPFO_mode,file,category);
    TH1F* h_rej_PiPi_cos_eff = efficiencyCorrection(h_rej_PiPi_cos,LPFO_mode,file,category);
    h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos_eff,LPFO_mode,file,category);
    h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos_eff,LPFO_mode,file,category);
  }else{
    h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
    h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);

  const Int_t nbins = h_reco_LPFO_scos_eff_corr->GetNbinsX();

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

  TH1F *h_reco_LPFO_pq_cos;
  h_reco_LPFO_pq_cos = CorrectHist(prodMode, h_reco_LPFO_qcos_eff_corr, p_vec);
  // if( category!="bg" ){
  //   h_reco_LPFO_pq_cos = CorrectHist(prodMode, h_reco_LPFO_qcos_eff_corr, p_vec);
  // }else{
  //   h_reco_LPFO_pq_cos = (TH1F*) h_reco_LPFO_qcos_eff_corr->Clone();
  // }
  StyleHist(h_reco_LPFO_pq_cos,kBlue);

  // Fitting
  TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
  f_gen->SetParNames("S","A");
  h_gen_q_qcos->Fit("f_gen","MNRS");
  cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

  TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
  f_reco->SetParNames("S","A");
  h_reco_LPFO_pq_cos->Fit("f_reco","MNRS");
  cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;

  // output
  unordered_map<TString, TH1F*> hmap;
  hmap["gen"] = h_gen_q_qcos;
  hmap["reco"] = h_reco_LPFO_pq_cos;

  return hmap;


}

void pq_method_PiLPFO_total()
{
  try
  {
    unordered_map<TString, unordered_map<TString, TFile*> > file_map;
    unordered_map<TString, unordered_map<TString, unordered_map< TString, unordered_map<TString, TH1F*> > > > hmap;
    for( auto process : processes ){
      for( auto chiral : chirals ){
        Int_t processID = production.at({process,chiral});
        cout << process << " " << chiral << " " << processID << endl;
        TString filename = inputDir + "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I" + processID + "." + process + "." + chiral + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root";
        TFile *file = new TFile(filename,"READ");
        if( !file->IsOpen() ) throw std::runtime_error("File not found");

        if( process=="P2f_z_h" ){
          for( auto category : qqbars ){
            hmap[process][chiral][category] = main_pq(file,category);
          }
        }else{
          hmap[process][chiral]["bg"] = main_pq(file,"bg");
        }
        // file_map[process][chiral] = file;
      }
    }

    unordered_map<TString, THStack*> hs_reco;
    for( auto chiral : chirals ){
      hs_reco[chiral] = new THStack("hs_reco_" + chiral,";cos#theta;Entries");
    }

    for( auto process : processes ){
      for( auto chiral : chirals ){
        if( process=="P2f_z_h" ){
          for( auto category : qqbars ){
            // if(category=="bb" || category=="cc" || category=="ss") continue;
            if(category=="bb" || category=="cc" ) continue;
            TH1F *h = hmap[process][chiral][category]["reco"];
            h->GetYaxis()->SetRangeUser(0,TopRange);
            h->SetFillStyle(0);
            h->SetTitle(histLabel(process,category));
            hs_reco.at(chiral)->Add(h);
          }
        }else{
          TH1F *h = hmap[process][chiral]["bg"]["reco"];
          h->GetYaxis()->SetRangeUser(0,TopRange);
          h->SetFillStyle(0);
          h->SetTitle(histLabel(process,"bg"));
          h->SetLineStyle(7);
          hs_reco.at(chiral)->Add(h);
        }
      }
    }

    for( auto chiral : chirals ){
      TCanvas *c_hs_reco = new TCanvas("c_hs_reco_" + chiral,"c_hs_reco_" + chiral,900,900);
      TPad *pad_hs_reco = new TPad("pad_hs_reco_" + chiral, "pad_hs_reco_" + chiral,0,0,1,1);
      StylePad(pad_hs_reco,0,0.12,0,0.15);
      gStyle->SetHistTopMargin(0);
      gStyle->SetPalette(55);
      hs_reco.at(chiral)->Draw("h plc nostack");
      pad_hs_reco->BuildLegend();
    }

    // TCanvas *c_hs_reco = new TCanvas("c_hs_reco","c_hs_reco",900,900);
    // TPad *pad_hs_reco = new TPad("pad_hs_reco", "pad_hs_reco",0,0,1,1);
    // StylePad(pad_hs_reco,0,0.12,0,0.15);
    // gStyle->SetPalette(55);
    // hs_reco->Draw("h plc nostack");
    // pad_hs_reco->BuildLegend();


  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}