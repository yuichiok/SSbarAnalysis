#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString LPFO_mode = "Pi";
TString chiral = "eL.pR";
// TString chiral = "eR.pL";

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

unordered_map<TString, TH1F*> main_pq(TFile* file, TString prodMode)
{
  gStyle->SetOptStat(0);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get(prodMode + "/gen/h_" + prodMode + "_qcos");
  TH1F *h_reco_LPFO_scos  = (TH1F*) file->Get(prodMode + "/cos/h_" + prodMode + "_" + LPFO_mode + "_scos");
  TH1F *h_reco_LPFO_qcos  = (TH1F*) file->Get(prodMode + "/cos/h_" + prodMode + "_" + LPFO_mode + "_qcos");

  cout << prodMode << " reco entry = " << h_reco_LPFO_qcos->GetEntries() << endl;
  cout << prodMode << " gen entry  = " << h_gen_q_qcos->GetEntries() << endl;
  cout << prodMode << " reco eff   = " << (float)h_reco_LPFO_qcos->GetEntries() / (float)h_gen_q_qcos->GetEntries() << endl;
  
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
    // h_reco_LPFO_scos_eff_corr = resolutionCorrection(h_reco_LPFO_scos_eff,LPFO_mode,file,prodMode);
    // h_reco_LPFO_qcos_eff_corr = resolutionCorrection(h_reco_LPFO_qcos_eff,LPFO_mode,file,prodMode);
    h_reco_LPFO_scos_eff_corr = (TH1F*) h_reco_LPFO_scos_eff->Clone();
    h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos_eff->Clone();

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
    // h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos_eff,LPFO_mode,file,prodMode);
    // h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos_eff,LPFO_mode,file,prodMode);
    h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos_eff->Clone();
    h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos_eff->Clone();

  }else{
    h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
    h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+2);

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

  TH1F *h_reco_LPFO_pq_cos = CorrectHist(prodMode, h_reco_LPFO_qcos_eff_corr, p_vec);
  // TH1F *h_reco_LPFO_pq_cos = (TH1F*) h_reco_LPFO_qcos_eff_corr->Clone();
  StyleHist(h_reco_LPFO_pq_cos,kBlue);

  // Fitting
  TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
  f_gen->SetParNames("S","A");
  f_gen->SetParameters(1.84809E-2+1.84846E-2,3.47499E-2-3.10189E-2);
  h_gen_q_qcos->Fit("f_gen","MNRS");
  cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

  TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
  f_reco->SetParNames("S","A");
  h_reco_LPFO_pq_cos->Fit("f_reco","MNRS");
  cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;


  // Int_t luminosity = production.at({"P2f_z_h",chiral}).second;
  // h_gen_q_qcos->Scale(1.0 / luminosity);
  // h_reco_LPFO_pq_cos->Scale(1.0 / luminosity);


  // output
  unordered_map<TString, TH1F*> hmap;
  hmap["gen"] = h_gen_q_qcos;
  hmap["reco"] = h_reco_LPFO_pq_cos;

  return hmap;


}

void pq_method_PiLPFO_add()
{
  TGaxis::SetMaxDigits(3);

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

    unordered_map<TString, unordered_map<TString, TH1F*>> hmap; // [ud][reco/gen]
    hmap["uu"] = main_pq(file_map["P2f_z_h"][chiral], "uu");
    hmap["dd"] = main_pq(file_map["P2f_z_h"][chiral], "dd");

    Int_t Nreco = hmap.at("uu").at("reco")->Integral() + hmap.at("dd").at("reco")->Integral();

    // cout << "uu ratio = " << uu_ratio << endl;
    // cout << "dd ratio = " << dd_ratio << endl;
    // cout << "weight   = " << dd_ratio / uu_ratio << endl;

    // hmap.at("dd").at("gen")->Scale( dd_ratio / uu_ratio );
    // hmap.at("dd").at("gen")->Scale( 0.85784 );
    // hmap.at("dd").at("gen")->Scale( 3618.54 / 5307.50 );

    Float_t eff_uu = (Float_t) hmap["uu"]["gen"]->Integral() / (Float_t) hmap["uu"]["reco"]->Integral();
    Float_t eff_dd = (Float_t) hmap["dd"]["gen"]->Integral() / (Float_t) hmap["dd"]["reco"]->Integral();

    cout << "=========\n";
    cout << "gen uu = " << hmap["uu"]["gen"]->Integral() << "\n";
    cout << "rec uu = " << hmap["uu"]["reco"]->Integral() << "\n";
    cout << "=========\n";
    cout << "gen dd = " << hmap["dd"]["gen"]->Integral() << "\n";
    cout << "rec dd = " << hmap["dd"]["reco"]->Integral() << "\n";
    cout << "=========\n";
    cout << "eff uu = " << eff_uu << ", eff dd = " << eff_dd << "\n";

    hmap.at("uu").at("reco")->Scale(eff_uu);
    hmap.at("dd").at("reco")->Scale(eff_dd);


    TH1F *h_reco = (TH1F*) hmap.at("uu").at("reco")->Clone();
    h_reco->Add(hmap.at("dd").at("reco"));
    StyleHist(h_reco,kBlue+2);
    
    TH1F *h_gen = (TH1F*) hmap.at("uu").at("gen")->Clone();
    h_gen->Add(hmap.at("dd").at("gen"));
    StyleHist(h_gen,kGreen+2);

    // Normalize(h_reco);
    // Normalize(h_gen);


    TCanvas *cts_reco = new TCanvas("cts_reco","cts_reco",800,800);
    TCanvas *cts_gen  = new TCanvas("cts_gen ","cts_gen ",800,800);

    // Total plots
      // Normalization
      // Normalize2Gen(h_gen,h_reco);

      // Fitting
      TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
      f_gen->SetParNames("S","A");
      h_gen->Fit("f_gen","MNRS");
      f_gen->SetLineColor(kGreen+2);
      f_gen->SetLineStyle(2);
      cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

      TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
      f_reco->SetParNames("S","A");
      h_reco->Fit("f_reco","MNRS");
      f_reco->SetLineColor(kRed);
      cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;

      // output AFB
      // Int_t Nreco = h_reco->Integral();
      cout << "Nreco = " << Nreco << endl;
      Float_t AFB_gen  = AFB_calculation(f_gen);
      Float_t AFB_reco = AFB_calculation(f_reco);
      Float_t AFB_reco_error = AFB_error(AFB_reco, Nreco);

      cout << "Gen  AFB = " << AFB_gen << endl;
      cout << "Reco AFB = " << AFB_reco << " +- " << AFB_reco_error << endl;
      cout << "Precision = " << AFB_reco_error / AFB_reco << endl;

      TCanvas *cTotal = new TCanvas("cTotal","cTotal",800,800);
      TPad *padTotal = new TPad("padTotal","padTotal",0,0,1,1);
      StylePad(padTotal,0,0.15,0,0.17);
      // h_gen->GetYaxis()->SetRangeUser(0,800E3);
      h_gen->SetTitle(";cos#theta_{#pi^{-}};Entries");
      h_gen->GetYaxis()->SetRangeUser(0,1E6);
      h_gen->Draw("h");
      h_reco->Draw("same");
      // f_reco->Draw("same");
      // f_gen->Draw("same");

/*
    Int_t bineval = 21;
    cout << h_reco->GetBinError(bineval) / h_reco->GetBinContent(bineval) << endl;


      // TLegend *legTotal = new TLegend(0.51,0.74,0.89,0.89);
      TLegend *legTotal = new TLegend(0.51,0.74,0.89,0.84);
      legTotal->SetMargin(0.4);
      legTotal->SetLineColor(0);
      legTotal->AddEntry(h_gen,"Parton level","f");
      legTotal->AddEntry(h_reco,"Reconstructed","lE");
      // legTotal->AddEntry(f_gen,"Generated Fit","l");
      // legTotal->AddEntry(f_reco,"Reconstructed Fit","l");
      legTotal->Draw();

    // Stack plots
      // Styling
      StyleHist(hmap.at("uu").at("reco"),kViolet);
      StyleHist(hmap.at("dd").at("reco"),kRed);

      cts_reco->cd();
      TH1F *totalReco = (TH1F*) hmap.at("uu").at("reco")->Clone();
      totalReco->Add(hmap.at("dd").at("reco"));
      StyleHist(totalReco,kBlue);
      totalReco->SetFillStyle(0);
      TPad *padts_reco = new TPad("padts_reco","padts_reco",0,0,1,1);
      StylePad(padts_reco,0,0.15,0,0.17);
      // totalReco->GetYaxis()->SetRangeUser(0,800E3);
      totalReco->SetTitle(";cos#theta;Entries");
      totalReco->Draw("h");
      hmap.at("uu").at("reco")->Draw("h same");
      hmap.at("dd").at("reco")->Draw("h same");

      TLegend *legStackReco = new TLegend(0.3,0.7,0.7,0.85);
      legStackReco->SetMargin(0.4);
      legStackReco->SetLineColor(0);
      legStackReco->AddEntry(totalReco,"Total Reconstructed #pi^{-}","l");
      legStackReco->AddEntry(hmap.at("uu").at("reco"),"u#bar{u}","l");
      legStackReco->AddEntry(hmap.at("dd").at("reco"),"d#bar{d}","l");
      legStackReco->Draw();


      cts_gen->cd();
      TH1F *h_gen_uu = (TH1F*) hmap.at("uu").at("gen")->Clone();
      TH1F *h_gen_dd = (TH1F*) hmap.at("dd").at("gen")->Clone();
      TH1F *totalGen = (TH1F*) h_gen_uu->Clone();
      totalGen->Add(h_gen_dd);

      double intCosReco = totalGen->Integral(20,80);
      double intCosGen  = h_reco->Integral(20,80);
      double recogenRatio = intCosGen / intCosReco;
      h_gen_uu->Scale( recogenRatio );
      h_gen_dd->Scale( recogenRatio );
      totalGen->Scale( recogenRatio );

      StyleHist(totalGen,kGreen+2);
      StyleHist(h_gen_uu,kViolet);
      StyleHist(h_gen_dd,kRed);
      totalGen->SetFillStyle(0);
      TPad *padts_gen = new TPad("padts_gen","padts_gen",0,0,1,1);
      StylePad(padts_gen,0,0.15,0,0.17);
      // totalGen->GetYaxis()->SetRangeUser(0,800E3);
      totalGen->SetTitle(";cos#theta;Entries");
      totalGen->Draw("h");
      h_gen_uu->Draw("h same");
      h_gen_dd->Draw("h same");

      TLegend *legStackGen = new TLegend(0.51,0.76,0.85,0.86);
      legStackGen->SetMargin(0.4);
      legStackGen->SetLineColor(0);
      legStackGen->AddEntry(totalGen,"Total Generated #pi^{-}","l");
      legStackGen->AddEntry(h_gen_uu,"u#bar{u}","l");
      legStackGen->AddEntry(h_gen_dd,"d#bar{d}","l");
      legStackGen->Draw();

    // Draw polar angle fit ratio
      TCanvas  *c_ratio_fit = new TCanvas("c_ratio_fit","c_ratio_fit",700,900);
      // Int_t binFitRange = fitRange * h_reco->GetNbinsX();
      TH1F *h_reco_LPFO_pq_cos_subhist = new TH1F("h_reco_LPFO_pq_cos_subhist",";LPFO Pion cos#theta; Entries",90,-fitRange,fitRange);
      for ( int ibin=1; ibin<=h_reco_LPFO_pq_cos_subhist->GetNbinsX(); ibin++ )
      {
        Int_t recobin = ibin + 5;
        h_reco_LPFO_pq_cos_subhist->SetBinContent(ibin,h_reco->GetBinContent(recobin));
        h_reco_LPFO_pq_cos_subhist->SetBinError(ibin,h_reco->GetBinError(recobin));
      }
      h_reco_LPFO_pq_cos_subhist->SetMarkerStyle(20);

      TF1 *f_reco_ratio = new TF1("f_reco_ratio","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
      f_reco_ratio->SetParNames("S","A");
      h_reco_LPFO_pq_cos_subhist->Fit("f_reco_ratio");
      cout << "Reco Chi2 / ndf = " << f_reco_ratio->GetChisquare() << " / " << f_reco_ratio->GetNDF() << endl;
      // h_reco_LPFO_pq_cos_subhist->GetYaxis()->SetRangeUser(0,800E3);
      c_ratio_fit->Clear();

      auto trp = new TRatioPlot(h_reco_LPFO_pq_cos_subhist);
      trp->SetGraphDrawOpt("P");
      trp->SetSeparationMargin(0.0);
      trp->Draw();
      trp->GetLowerRefYaxis()->SetTitle("Ratio");
      trp->GetLowerRefGraph()->SetMinimum(-4);
      trp->GetLowerRefGraph()->SetMaximum(4);
      trp->GetLowerRefYaxis()->SetLabelSize(0.02);

      trp->GetUpperPad()->cd();
      TLegend *leg_trp = new TLegend(0.3,0.7,0.7,0.85);
      leg_trp->SetMargin(0.4);
      leg_trp->SetLineColor(0);
      leg_trp->AddEntry(h_reco_LPFO_pq_cos_subhist,"Reconstructed #pi^{-} (corrected)","p");
      leg_trp->AddEntry(f_reco_ratio,"#frac{d#sigma}{dcos#theta} fit for LPFO","l");
      leg_trp->Draw();

    // Draw polar angle fit ratio
      TCanvas  *c_gen_ratio_fit = new TCanvas("c_gen_ratio_fit","c_gen_ratio_fit",700,900);
      // Int_t binFitRange = fitRange * h_reco->GetNbinsX();
      TH1F *h_gen_LPFO_pq_cos_subhist = new TH1F("h_gen_LPFO_pq_cos_subhist",";LPFO Pion cos#theta; Entries",90,-fitRange,fitRange);
      // h_gen_LPFO_pq_cos_subhist->Set
      for ( int ibin=1; ibin<=h_gen_LPFO_pq_cos_subhist->GetNbinsX(); ibin++ )
      {
        Int_t recobin = ibin + 5;
        h_gen_LPFO_pq_cos_subhist->SetBinContent(ibin,h_gen->GetBinContent(recobin));
        h_gen_LPFO_pq_cos_subhist->SetBinError(ibin,h_gen->GetBinError(recobin));
      }
      h_gen_LPFO_pq_cos_subhist->SetMarkerStyle(20);

      TF1 *f_gen_ratio = new TF1("f_gen_ratio","[0]*(1+x*x)+[1]*x",-fitRange,fitRange);
      f_gen_ratio->SetParNames("S","A");
      f_gen_ratio->SetParameters(1.84809E-2+1.84846E-2,3.47499E-2-3.10189E-2);
      h_gen_LPFO_pq_cos_subhist->Fit("f_gen_ratio");
      cout << "Reco Chi2 / ndf = " << f_reco_ratio->GetChisquare() << " / " << f_reco_ratio->GetNDF() << endl;
      // h_reco_LPFO_pq_cos_subhist->GetYaxis()->SetRangeUser(0,800E3);
      c_gen_ratio_fit->Clear();

      auto trp_gen = new TRatioPlot(h_gen_LPFO_pq_cos_subhist);
      trp_gen->SetGraphDrawOpt("P");
      trp_gen->SetSeparationMargin(0.0);
      trp_gen->Draw();
      trp_gen->GetLowerRefYaxis()->SetTitle("Ratio");
      trp_gen->GetLowerRefGraph()->SetMinimum(-4);
      trp_gen->GetLowerRefGraph()->SetMaximum(4);
      trp_gen->GetLowerRefYaxis()->SetLabelSize(0.02);

      trp_gen->GetUpperPad()->cd();
      TLegend *leg_trp_gen = new TLegend(0.3,0.7,0.7,0.85);
      leg_trp_gen->SetMargin(0.4);
      leg_trp_gen->SetLineColor(0);
      leg_trp_gen->AddEntry(h_gen_LPFO_pq_cos_subhist,"Parton level #pi^{-} (corrected)","p");
      leg_trp_gen->AddEntry(f_gen_ratio,"#frac{d#sigma}{dcos#theta} fit for LPFO","l");
      leg_trp_gen->Draw();


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
*/

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}