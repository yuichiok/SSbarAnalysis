#include <iostream>
#include <vector>
#include <map>

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
  TH1F *h_gen_q_qcos     = (TH1F*) file->Get("h_gen_q_qcos");
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

  TH1F *h_reco_LPFO_pq_cos = CorrectHist(h_reco_LPFO_qcos_eff_corr, p_vec);
  StyleHist(h_reco_LPFO_pq_cos,kBlue);

  Normalize2Gen(h_gen_q_qcos,h_reco_LPFO_qcos_eff_corr);

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
  TString subh_name = "subh_" + (TString)h_reco_LPFO_qcos->GetName();
  TCanvas  *c_ratio = new TCanvas(c_name,c_name,800,800);
  TH1F *h_reco_LPFO_pq_cos_subhist = new TH1F(subh_name,";LPFO Pion cos#theta; Entries",80,-0.8,0.8);
  for ( int ibin=1; ibin<=h_reco_LPFO_pq_cos_subhist->GetNbinsX(); ibin++ )
  {
    Int_t recobin = ibin + 10;
    h_reco_LPFO_pq_cos_subhist->SetBinContent(ibin,h_reco_LPFO_pq_cos->GetBinContent(recobin));
    h_reco_LPFO_pq_cos_subhist->SetBinError(ibin,h_reco_LPFO_pq_cos->GetBinError(recobin));
  }
  h_reco_LPFO_pq_cos_subhist->SetMarkerStyle(20);

  TF1 *f_reco_ratio = new TF1("f_reco_ratio","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_reco_ratio->SetParNames("S","A");
  h_reco_LPFO_pq_cos_subhist->Fit("f_reco_ratio");
  c_ratio->Clear();

  auto trp = new TRatioPlot(h_reco_LPFO_pq_cos_subhist);
  trp->SetGraphDrawOpt("P");
  trp->SetSeparationMargin(0.0);
  trp->Draw();
  trp->GetLowerRefYaxis()->SetTitle("Ratio");
  trp->GetLowerRefGraph()->SetMinimum(-4);
  trp->GetLowerRefGraph()->SetMaximum(4);

  trp->GetUpperPad()->cd();
  TLegend *leg_trp = new TLegend(0.3,0.7,0.7,0.85);
  leg_trp->SetMargin(0.4);
  leg_trp->SetLineColor(0);
  leg_trp->AddEntry(h_reco_LPFO_pq_cos_subhist,"Reconstructed #pi^{-} (corrected)","p");
  leg_trp->AddEntry(f_reco_ratio,"#frac{d#sigma}{dcos#theta} fit for LPFO","l");
  leg_trp->Draw();

  return c_ratio;

}

void SaveHists(TFile *file, TString prod_mode, TString LPFO_mode, vector<TH1F*> hvec)
{
  if (!file->IsOpen()) return;

  for ( auto ih : hvec )
  {
    TCanvas *c = main_pq(file, ih, LPFO_mode);
    TString printname = "c_" + prod_mode + "_" + (TString)ih->GetName() + ".png";
    c->Print("~/Desktop/" + printname);
  }
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

void pq_method_PiLPFO_eff()
{
  TString prod_mode = "ud";
  TString LPFO_mode = "Pi";
  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR." + prod_mode + ".KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.eff.hists.all.root","READ");

  try
  {
    if (!file->IsOpen()) return;
    vector<TH1F*> hvec = GetHists(file, "cos_cut_eff");
    SaveHists(file, prod_mode, LPFO_mode, hvec);
    PrintEfficiency(file, hvec);

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}