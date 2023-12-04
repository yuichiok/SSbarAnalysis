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
// Float_t TopRange = 700;
Float_t TopRange = 550;
// Float_t TopRange = 0.7;

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
// array<TString,4> processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"};
array<TString,4> processes = {"Pqqh", "P4f_zz_h", "P4f_ww_h", "P2f_z_h"};

// array<TString,6> qqbars    = {"rr", "bb", "cc", "ss", "dd", "uu"};
array<TString,6> qqbars    = {"rr", "bb", "cc", "dd", "uu", "ss"};

array<TString,4> leg_processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"};
array<TString,6> leg_qqbars    = {"dd", "uu", "ss", "cc", "bb", "rr"};



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
  }else if (process=="Pqqh"){
    return "q#bar{q}H";
  }
  return "";

}

unordered_map<TString, TH1F*> main_pq(TFile* file, TString process, TString chiral, TString category)
{
  gStyle->SetOptStat(0);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get(category + "/gen/h_" + category + "_qcos");
  TH1F *h_reco_LPFO_scos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_scos");
  TH1F *h_reco_LPFO_qcos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_qcos");

  cout << category << " reco entry = " << h_reco_LPFO_qcos->GetEntries() << endl;
  cout << category << " gen entry  = " << h_gen_q_qcos->GetEntries() << endl;
  cout << category << " reco eff   = " << (float)h_reco_LPFO_qcos->GetEntries() / (float)h_gen_q_qcos->GetEntries() << endl;
  
  // used for pq correction
  TH1F *h_acc_LL_cos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_acc_cos");
  TH1F *h_rej_LL_cos  = (TH1F*) file->Get(category + "/cos/h_" + category + "_" + LPFO_mode + "_rej_cos");

  // efficiency correction
  Bool_t isEffCorr = true;
  TH1F *h_reco_LPFO_scos_eff_corr;
  TH1F *h_reco_LPFO_qcos_eff_corr;
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

  TH1F *h_acc_LL_cos_eff_corr;
  TH1F *h_rej_LL_cos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_acc_LL_cos_eff = efficiencyCorrection(h_acc_LL_cos,LPFO_mode,file,category);
    TH1F* h_rej_LL_cos_eff = efficiencyCorrection(h_rej_LL_cos,LPFO_mode,file,category);
    h_acc_LL_cos_eff_corr = resolutionCorrection(h_acc_LL_cos_eff,LPFO_mode,file,category);
    h_rej_LL_cos_eff_corr = resolutionCorrection(h_rej_LL_cos_eff,LPFO_mode,file,category);
  }else{
    h_acc_LL_cos_eff_corr = (TH1F*) h_acc_LL_cos->Clone();
    h_rej_LL_cos_eff_corr = (TH1F*) h_rej_LL_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);

  const Int_t nbins = h_reco_LPFO_scos_eff_corr->GetNbinsX();

  // pq correction
  TString pValName = "p_" + LPFO_mode + LPFO_mode + "_" + category;
  TH1F *p_LL = new TH1F(pValName,pValName, 50,0,1);
  p_LL->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_LL_cos_eff_corr, h_rej_LL_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_LL->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_LL->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_LPFO_pq_cos;
  h_reco_LPFO_pq_cos = CorrectHist(category, h_reco_LPFO_qcos_eff_corr, p_vec);
  // if( category!="bg" ){
  //   h_reco_LPFO_pq_cos = CorrectHist(category, h_reco_LPFO_qcos_eff_corr, p_vec);
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

  Int_t luminosity = production.at({process,chiral}).second;
  h_gen_q_qcos->Scale(1.0 / luminosity);
  h_reco_LPFO_pq_cos->Scale(1.0 / luminosity);

  // output
  unordered_map<TString, TH1F*> hmap;
  hmap["gen"] = h_gen_q_qcos;
  hmap["reco"] = h_reco_LPFO_pq_cos;

  return hmap;


}

void pq_method_LPFO_total()
{
  try
  {
    unordered_map<TString, unordered_map<TString, TFile*> > file_map;
    unordered_map<TString, unordered_map<TString, unordered_map< TString, unordered_map<TString, TH1F*> > > > hmap;
    for( auto process : processes ){
      for( auto chiral : chirals ){
        Int_t processID = production.at({process,chiral}).first;
        cout << process << " " << chiral << " " << processID << endl;
        TString filename = inputDir + "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I" + processID + "." + process + "." + chiral + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root";
        TFile *file = new TFile(filename,"READ");
        if( !file->IsOpen() ) throw std::runtime_error("File not found");

        if( process=="P2f_z_h" ){
          for( auto category : qqbars ){
            cout << "======== read " << category << " ===========" << endl;
            hmap[process][chiral][category] = main_pq(file,process,chiral,category);
          }
        }else{
          hmap[process][chiral]["bg"] = main_pq(file,process,chiral,"bg");
        }
        file_map[process][chiral] = file;
      }
    }

    unordered_map<TString, THStack*> hs_reco;
    unordered_map<TString, TH1F*>    h_reco;
    unordered_map<TString, TLegend*> leg_reco;

    for( auto chiral : chirals ){
      hs_reco[chiral] = new THStack("hs_reco_" + chiral,";cos#theta;Entries / Int. Lumi.");
      h_reco[chiral]  = new TH1F("h_reco_" + chiral,";cos#theta;Entries / Int. Lumi.", 100,-1,1);
      leg_reco[chiral] = new TLegend(0.59,0.65,0.89,0.85);
      leg_reco[chiral]->SetMargin(0.4);
      leg_reco[chiral]->SetBorderSize(0);
      leg_reco[chiral]->SetFillStyle(0);
    }

    for( auto process : processes ){
      for( auto chiral : chirals ){
        if( process=="P2f_z_h" ){
          for( auto category : qqbars ){
            TH1F *h = hmap[process][chiral][category]["reco"];
            h->GetYaxis()->SetRangeUser(0,TopRange);
            h->SetFillStyle(0);
            TString histName = histLabel(process,category);
            h->SetTitle(histName);
            hs_reco.at(chiral)->Add(h);
            h_reco.at(chiral)->Add(h);

            cout << "====== " << category << " ======\n";
            cout << h->Integral(1,25) << endl;

          }
        }else{
          TH1F *h = hmap[process][chiral]["bg"]["reco"];
          h->GetYaxis()->SetRangeUser(0,TopRange);
          h->SetFillStyle(0);
          TString histName = histLabel(process,"bg");
          h->SetTitle(histName);
          h->SetLineStyle(7);
          hs_reco.at(chiral)->Add(h);
          h_reco.at(chiral)->Add(h);
        }
      }
    }
    for( auto process : leg_processes ){
      for( auto chiral : chirals ){
        if( process=="P2f_z_h" ){
          for( auto category : leg_qqbars ){
            TH1F *h = hmap[process][chiral][category]["reco"];
            TString histName = histLabel(process,category);
            leg_reco.at(chiral)->AddEntry(h,histName,"l");
          }
        }else{
          TH1F *h = hmap[process][chiral]["bg"]["reco"];
          TString histName = histLabel(process,"bg");
          leg_reco.at(chiral)->AddEntry(h,histName,"l");
        }
      }
    }


    cout << "========== Main Fit Results ==========\n";

    for( auto chiral : chirals ){
      TCanvas *c_hs_reco = new TCanvas("c_hs_reco_" + chiral,"c_hs_reco_" + chiral,900,900);
      TPad *pad_hs_reco = new TPad("pad_hs_reco_" + chiral, "pad_hs_reco_" + chiral,0,0,1,1);
      StylePad(pad_hs_reco,0,0.12,0,0.15);
      gStyle->SetHistTopMargin(0);
      gStyle->SetPalette(55);
      // hs_reco.at(chiral)->Draw("h plc nostack");
      hs_reco.at(chiral)->Draw("h plc");
      leg_reco.at(chiral)->Draw("same");
      hs_reco.at(chiral)->SetMaximum(TopRange);



      // Fit
      TCanvas *c_h_reco = new TCanvas("c_h_reco_" + chiral,"c_h_reco_" + chiral,900,900);
      TPad *pad_h_reco  = new TPad("pad_h_reco_" + chiral, "pad_h_reco_" + chiral,0,0,1,1);
      StylePad(pad_h_reco,0,0.12,0,0.15);
      h_reco.at(chiral)->SetLineWidth(3);
      h_reco.at(chiral)->Draw("");

      TF1 * f_total = new TF1("f_total","[0]*(1+x*x)+[1]*x",-0.6,0.6);
      f_total->SetParNames("S","A");
      h_reco.at(chiral)->Fit("f_total","MNRS");
      cout << "Gen Chi2 / ndf = " << f_total->GetChisquare() << " / " << f_total->GetNDF() << endl;
      f_total->SetLineWidth(3);
      f_total->SetLineColor(kRed);
      f_total->Draw("same");

      TLegend *leg = new TLegend(0.56,0.71,0.86,0.85);
      leg->SetMargin(0.2);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(h_reco.at(chiral),"Data","lep");
      leg->AddEntry(f_total,"#frac{d#sigma}{dcos#theta} = S(1+cos^{2}#theta) + Acos#theta","l");
      leg->Draw("same");

      // TCanvas *c_hs_reco_nostack = new TCanvas("c_hs_reco_nostack_" + chiral,"c_hs_reco_nostack_" + chiral,900,900);
      // TPad *pad_hs_reco_nostack = new TPad("pad_hs_reco_nostack_" + chiral, "pad_hs_reco_nostack_" + chiral,0,0,1,1);
      // StylePad(pad_hs_reco_nostack,0,0.12,0,0.15);
      // hs_reco.at(chiral)->Draw("h plc nostack");
      // leg_reco.at(chiral)->Draw();
      // hs_reco.at(chiral)->SetMaximum(TopRange);

      cout << "======================================\n";

    }




  /*
    unordered_map<TString, THStack*> hs_gen;
    for( auto chiral : chirals ){
      hs_gen[chiral] = new THStack("hs_gen_" + chiral,";cos#theta;Entries / Int. Lumi.");
    }

    for( auto process : processes ){
      for( auto chiral : chirals ){
        if( process=="P2f_z_h" ){
          for( auto category : qqbars ){
            if(category=="bb" || category=="cc" || category=="ss") continue;
            // if(category=="bb" || category=="cc" ) continue;
            TH1F *h = hmap[process][chiral][category]["gen"];
            // h->GetYaxis()->SetRangeUser(0,TopRange);
            h->SetFillStyle(0);
            h->SetTitle(histLabel(process,category));
            hs_gen.at(chiral)->Add(h);
          }
        }else{
          TH1F *h = hmap[process][chiral]["bg"]["gen"];
          // h->GetYaxis()->SetRangeUser(0,TopRange);
          h->SetFillStyle(0);
          h->SetTitle(histLabel(process,"bg"));
          h->SetLineStyle(7);
          hs_gen.at(chiral)->Add(h);
        }
      }
    }

    for( auto chiral : chirals ){
      TCanvas *c_hs_gen = new TCanvas("c_hs_gen_" + chiral,"c_hs_gen_" + chiral,900,900);
      TPad *pad_hs_gen = new TPad("pad_hs_gen_" + chiral, "pad_hs_gen_" + chiral,0,0,1,1);
      StylePad(pad_hs_gen,0,0.12,0,0.15);
      gStyle->SetHistTopMargin(0);
      gStyle->SetPalette(55);
      hs_gen.at(chiral)->Draw("h plc nostack");
      pad_hs_gen->BuildLegend(0.59,0.68,0.89,0.89);
    }
  */

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}