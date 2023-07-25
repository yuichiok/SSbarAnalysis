#include <iostream>
#include <vector>
#include <map>

#include "include/Styles.hh"
#include "include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector;

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
	double intCosReco = h->Integral(20,80);
	double intCosGen  = h_gen->Integral(20,80);

  cout << "intCosReco = " << intCosReco << ", intCosGen = " << intCosGen << "\n";

  h->Scale( intCosGen / intCosReco );
}

vector<TH1F*> GetHists(TFile *file, TString dir)
{
  vector<TH1F*> vec;
  file->cd(dir);
  for(auto k : *gDirectory->GetListOfKeys()) {
    TKey *key = static_cast<TKey*>(k);
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    vec.push_back( (TH1F*)key->ReadObj() );
  }
  return vec;
}

TCanvas * main_pq(TFile *file, TH1F *h_reco_LPFO_qcos, TString LPFO_mode)
{
  gStyle->SetOptStat(0);

  // reco and gen polar
  TH1F *h_gen_q_qcos_orig = (TH1F*) file->Get("h_gen_q_qcos");
  TH1F *h_gen_q_qcos      = (TH1F*) h_gen_q_qcos_orig->Clone();

  // TH1F *h_reco_LPFO_qcos = (TH1F*) file->Get("h_reco_" + LPFO_mode + "_qcos");

  // used for pq correction
  TH1F *h_acc_LPFO_cos  = (TH1F*) file->Get("pq/h_acc_" + LPFO_mode + LPFO_mode + "_cos");
  TH1F *h_rej_LPFO_cos  = (TH1F*) file->Get("pq/h_rej_" + LPFO_mode + LPFO_mode + "_cos");

  // efficiency correction
  Bool_t isEffCorr = true;
  TH1F *h_reco_LPFO_qcos_eff_corr;
  if (isEffCorr)
  {
    h_reco_LPFO_qcos_eff_corr = Efficiency_Correction(h_reco_LPFO_qcos, LPFO_mode, file);
  }else{
    h_reco_LPFO_qcos_eff_corr = (TH1F*) h_reco_LPFO_qcos->Clone();
  }

  TH1F *h_acc_LPFO_cos_eff_corr;
  TH1F *h_rej_LPFO_cos_eff_corr;
  if (isEffCorr)
  {
    h_acc_LPFO_cos_eff_corr = Efficiency_Correction(h_acc_LPFO_cos, LPFO_mode, file);
    h_rej_LPFO_cos_eff_corr = Efficiency_Correction(h_rej_LPFO_cos, LPFO_mode, file);
  }else{
    h_acc_LPFO_cos_eff_corr = (TH1F*) h_acc_LPFO_cos->Clone();
    h_rej_LPFO_cos_eff_corr = (TH1F*) h_rej_LPFO_cos->Clone();
  }

  StyleHist(h_gen_q_qcos,kGreen+1);
  h_gen_q_qcos->SetFillStyle(0);
  h_gen_q_qcos->SetLineStyle(2);

  StyleHist(h_reco_LPFO_qcos_eff_corr,kRed+2);
  StyleHist(h_acc_LPFO_cos_eff_corr,kRed+2);
  StyleHist(h_rej_LPFO_cos_eff_corr,kBlue+2);

  const Int_t nbins = h_reco_LPFO_qcos_eff_corr->GetNbinsX();

  TString p_name = "p_" + (TString)h_reco_LPFO_qcos->GetName();
  TH1F *p_LPFO = new TH1F(p_name, "p value", 50,0,1);
  p_LPFO->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_LPFO_cos_eff_corr, h_rej_LPFO_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_LPFO->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_LPFO->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_LPFO_pq_cos;
  TString hlast = "h_reco_" + LPFO_mode + "_cos_jet2_poff_pid_ud_spfo_chg";
  if ((TString) h_reco_LPFO_qcos_eff_corr->GetName() == hlast)
  {
    h_reco_LPFO_pq_cos = CorrectHist(h_reco_LPFO_qcos_eff_corr, p_vec);
  }else{
    h_reco_LPFO_pq_cos = (TH1F*) h_reco_LPFO_qcos_eff_corr->Clone();
  }  

  StyleHist(h_reco_LPFO_pq_cos,kBlue);

  // Normalize2Gen(h_gen_q_qcos,h_reco_LPFO_qcos_eff_corr);
  Normalize2Gen(h_gen_q_qcos,h_reco_LPFO_pq_cos);

  // Fitting
  TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_gen->SetParNames("S","A");
  h_gen_q_qcos->Fit("f_gen","MNRS");
  cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

  TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_reco->SetParNames("S","A");
  h_reco_LPFO_pq_cos->Fit("f_reco","MNRS");
  cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;
  
  // output AFB
  Float_t AFB_gen  = AFB_calculation(f_gen);
  Float_t AFB_reco = AFB_calculation(f_reco);
  cout << "Gen  AFB = " << AFB_gen << endl;
  cout << "Reco AFB = " << AFB_reco << endl;

  // Draw polar angle fit ratio
  TString c_name = "c_" + (TString)h_reco_LPFO_qcos->GetName();
  TCanvas  *c_dmc = new TCanvas(c_name,c_name,600,800);

  TPad *pad1 = new TPad("pad1", "pad1", 0.00, 0.3, 1.00, 1.00);
  TPad *pad2 = new TPad("pad2", "pad2", 0.00, 0.00, 1.00, 0.3);

  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad1->SetLeftMargin(0.14);
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGrid(1,1);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderMode(0);
  pad2->SetLeftMargin(0.14);
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetGrid(1,1);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();

  h_reco_LPFO_pq_cos->SetTitle(";cos#theta;Entries");
  h_reco_LPFO_pq_cos->Draw("h");
  h_gen_q_qcos->Draw("hsame");

  pad2->cd();

  TH1F *h_ratio = (TH1F*) h_reco_LPFO_qcos_eff_corr->Clone();
  cout << "reco integral = " << h_reco_LPFO_qcos->Integral() << ", gen integral = " << h_gen_q_qcos_orig->Integral() << endl;
  h_ratio->Divide(h_gen_q_qcos_orig);
  // h_ratio->GetYaxis()->SetRangeUser(0.8,1.2);

  h_ratio->SetTitle(";cos#theta;Data/MC");
  h_ratio->GetXaxis()->SetTitleSize(0.08);
  h_ratio->GetXaxis()->SetLabelSize(0.07);
  h_ratio->GetYaxis()->SetTitleSize(0.07);
  h_ratio->GetYaxis()->SetTitleOffset(0.75);
  h_ratio->GetYaxis()->SetLabelSize(0.07);
  h_ratio->Draw("P");

  TString iname = "h_reco_Pi_cos_jet2_poff_pid_ud_spfo_chg";
  TString pname = (TString) h_reco_LPFO_qcos->GetName();
  if(pname == iname){
    TCanvas *c_test = new TCanvas("c_test","c_test",800,800);
    h_ratio->Draw("h");
  }

  return c_dmc;

}

void SaveHists(TCanvas *c, TH1F *ih)
{
  TString printname = "c_" + prod_mode + "_" + (TString)ih->GetName() + ".png";
  c->Print("~/Desktop/" + printname);
}

void PrintEfficiency(TFile *file, vector<TH1F*> hvec)
{
  if (!file->IsOpen()) return;
  TH1F *h_gen_q_qcos     = (TH1F*) file->Get("h_gen_q_qcos");
  Int_t n_gen_events     = h_gen_q_qcos->GetEntries();
  cout << "name,nevents,efficiency\n";
  cout << "gen," << n_gen_events << ",-\n";
  for ( auto ih : hvec )
  {
    Int_t n_reco_events = ih->GetEntries();
    cout << ih->GetName() << ",";
    cout << ih->GetEntries() << "\n";
  }

}

void pq_method_PiLPFO_eff_DataMC()
{
  TGaxis::SetMaxDigits(3);

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR." + prod_mode + ".KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.spfox.eff.hists.all.root","READ");

  try
  {
    if (!file->IsOpen()) return;
    vector<TH1F*> hvec = GetHists(file, "cos_cut_eff");
    for ( auto ih : hvec )
    {
      TCanvas *c = main_pq(file, ih, LPFO_mode);
      // SaveHists(c, ih);
    }
    PrintEfficiency(file, hvec);

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}