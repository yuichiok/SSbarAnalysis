#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString LPFO_mode = "Pi";
// TString ichiral = "eL.pR";
TString ichiral = "eR.pL";

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
array<TString,4> processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"};
array<TString,6> qqbars    = {"dd", "uu", "ss", "cc", "bb", "rr"};

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

Float_t fitRange = 0.78;

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

unordered_map<TString, TH1F*> main_pq(TFile* file, TString process, TString chiral, TString qq)
{
  gStyle->SetOptStat(0);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get(qq + "/gen/h_" + qq + "_qcos");
  TH1F *h_reco_LPFO_scos  = (TH1F*) file->Get(qq + "/cos/h_" + qq + "_" + LPFO_mode + "_scos");
  TH1F *h_reco_LPFO_qcos  = (TH1F*) file->Get(qq + "/cos/h_" + qq + "_" + LPFO_mode + "_qcos");
  h_gen_q_qcos->Sumw2();

  cout << qq << " reco entry = " << h_reco_LPFO_qcos->GetEntries() << endl;
  cout << qq << " gen entry  = " << h_gen_q_qcos->GetEntries() << endl;
  cout << qq << " reco eff   = " << (float)h_reco_LPFO_qcos->GetEntries() / (float)h_gen_q_qcos->GetEntries() << endl;
  
  // used for pq correction
  TH1F *h_acc_LL_cos  = (TH1F*) file->Get(qq + "/cos/h_" + qq + "_" + LPFO_mode + "_acc_cos");
  TH1F *h_rej_LL_cos  = (TH1F*) file->Get(qq + "/cos/h_" + qq + "_" + LPFO_mode + "_rej_cos");

  // efficiency correction
  Bool_t isEffCorr = true;
  TH1F *h_reco_LPFO_scos_eff_corr;
  TH1F *h_reco_LPFO_qcos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_reco_LPFO_scos_eff = efficiencyCorrection(h_reco_LPFO_scos,LPFO_mode,file,qq);
    TH1F* h_reco_LPFO_qcos_eff = efficiencyCorrection(h_reco_LPFO_qcos,LPFO_mode,file,qq);
    // h_reco_LPFO_scos_eff_corr = resolutionCorrection(h_reco_LPFO_scos_eff,LPFO_mode,file,qq);
    // h_reco_LPFO_qcos_eff_corr = resolutionCorrection(h_reco_LPFO_qcos_eff,LPFO_mode,file,qq);
    h_reco_LPFO_scos_eff_corr = (TH1F*) h_reco_LPFO_scos_eff->Clone();
    h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos_eff->Clone();
  }else{
    h_reco_LPFO_scos_eff_corr = (TH1F*) h_reco_LPFO_scos->Clone();
    h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos->Clone();
  }

  TH1F *h_acc_LL_cos_eff_corr;
  TH1F *h_rej_LL_cos_eff_corr;
  if (isEffCorr)
  {
    TH1F* h_acc_LL_cos_eff = efficiencyCorrection(h_acc_LL_cos,LPFO_mode,file,qq);
    TH1F* h_rej_LL_cos_eff = efficiencyCorrection(h_rej_LL_cos,LPFO_mode,file,qq);
    // h_acc_LL_cos_eff_corr = resolutionCorrection(h_acc_LL_cos_eff,LPFO_mode,file,qq);
    // h_rej_LL_cos_eff_corr = resolutionCorrection(h_rej_LL_cos_eff,LPFO_mode,file,qq);
    h_acc_LL_cos_eff_corr = (TH1F*) h_acc_LL_cos_eff->Clone();
    h_rej_LL_cos_eff_corr = (TH1F*) h_rej_LL_cos_eff->Clone();


  }else{
    h_acc_LL_cos_eff_corr = (TH1F*) h_acc_LL_cos->Clone();
    h_rej_LL_cos_eff_corr = (TH1F*) h_rej_LL_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);

  const Int_t nbins = h_reco_LPFO_scos_eff_corr->GetNbinsX();

  // pq correction
  TString pValName = "p_" + LPFO_mode + LPFO_mode + "_" + qq;
  TH1F *p_LL = new TH1F(pValName,pValName, 50,0,1);
  p_LL->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_LL_cos_eff_corr, h_rej_LL_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_LL->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_LL->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_LPFO_pq_cos;
  h_reco_LPFO_pq_cos = CorrectHist(qq, h_reco_LPFO_qcos_eff_corr, p_vec);
  // if( qq!="bg" ){
  //   h_reco_LPFO_pq_cos = CorrectHist(qq, h_reco_LPFO_qcos_eff_corr, p_vec);
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

void pq_method_PiLPFO_BG_sub()
{
  try
  {
    // unordered_map<TString, unordered_map<TString, TFile*> > file_map;
    unordered_map<TString, unordered_map<TString, unordered_map<TString, TH1F*>>> hmap; // [process][ud][reco/gen]
    TH1F * h_total_reco     = new TH1F("h_total_reco_" + ichiral,";cos#theta;Entries / Int. Lumi.", 100,-1,1);
    TH1F * h_total_gen      = new TH1F("h_total_gen_" + ichiral,";cos#theta;Entries / Int. Lumi.", 100,-1,1);
    TH1F * h_total_reco_bg  = new TH1F("h_total_reco_bg_" + ichiral,";cos#theta;Entries / Int. Lumi.", 100,-1,1);
    TH1F * h_total_gen_bg   = new TH1F("h_total_gen_bg_" + ichiral,";cos#theta;Entries / Int. Lumi.", 100,-1,1);
    h_total_reco->Sumw2();
    h_total_gen->Sumw2();
    h_total_reco_bg->Sumw2();
    for( auto process : processes ){

      Int_t processID = production.at({process,ichiral}).first;
      cout << process << " " << ichiral << " " << processID << endl;
      TString filename = inputDir + "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I" + processID + "." + process + "." + ichiral + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root";
      TFile *file = new TFile(filename,"READ");
      if( !file->IsOpen() ) throw std::runtime_error("File not found");
      // file_map[process][chiral] = file;

      if(process=="P2f_z_h"){

        for( auto qq : qqbars ){

          hmap[process][qq] = main_pq(file,process,ichiral,qq);
          if(qq=="uu"||qq=="dd"){
            Float_t eff = (Float_t) hmap[process][qq]["gen"]->Integral() / (Float_t) hmap[process][qq]["reco"]->Integral();
            hmap.at(process).at(qq).at("reco")->Scale(eff);
          }else{
            h_total_reco_bg->Add(hmap.at(process).at(qq).at("reco"));
            h_total_gen_bg->Add(hmap.at(process).at(qq).at("gen"));
          }

          h_total_reco->Add(hmap.at(process).at(qq).at("reco"));
          h_total_gen->Add(hmap.at(process).at(qq).at("gen"));

        }

      }else{

        TString qq = "bg";
        hmap[process][qq] = main_pq(file,process,ichiral,qq);
        h_total_reco->Add(hmap.at(process).at(qq).at("reco"));
        h_total_gen->Add(hmap.at(process).at(qq).at("gen"));
        h_total_reco_bg->Add(hmap.at(process).at(qq).at("reco"));
        h_total_gen_bg->Add(hmap.at(process).at(qq).at("gen"));

      }

    }


    TCanvas *c_total = new TCanvas("c_total","c_total",800,800);
    TPad *padTotal = new TPad("padTotal","padTotal",0,0,1,1);
    StylePad(padTotal,0,0.15,0,0.17);
    h_total_reco->SetTitle(";cos#theta;Entries");

    TH1F* h_total_reco_sig = (TH1F*) h_total_reco->Clone();
    h_total_reco_sig->Add(h_total_reco_bg,-1);
    TH1F* h_total_gen_sig = (TH1F*) h_total_gen->Clone();
    h_total_gen_sig->Add(h_total_gen_bg,-1);

    // Normalize(h_total_reco_sig);
    // Normalize(h_total_gen_sig);

    StyleHist(h_total_reco_sig,kBlue+2);
    StyleHist(h_total_gen_sig,kGreen+2);

    h_total_gen_sig->GetYaxis()->SetRangeUser(0,0.05);
    h_total_gen_sig->Draw("h");

    // Fitting
    TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
    f_gen->SetParNames("S","A");
    h_total_gen_sig->Fit("f_gen","MNRS");
    f_gen->SetLineColor(kGreen+2);
    f_gen->SetLineStyle(2);
    cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

    TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
    f_reco->SetParNames("S","A");
    h_total_reco_sig->Fit("f_reco","MNRS");
    f_reco->SetLineColor(kRed);
    cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;
    f_gen->Draw("same");

    h_total_reco_sig->Draw("same");
    f_reco->Draw("same");

    TLegend *legTotal = new TLegend(0.51,0.74,0.89,0.89);
    legTotal->SetMargin(0.4);
    legTotal->SetLineColor(0);
    legTotal->AddEntry(h_total_gen_sig,"Genrated #pi^{-}","f");
    legTotal->AddEntry(h_total_reco_sig,"Reconstructed #pi^{-}","lE");
    legTotal->AddEntry(f_gen,"Generated Fit","l");
    legTotal->AddEntry(f_reco,"Reconstructed Fit","l");
    legTotal->Draw();




  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}