#include <iostream>
#include <vector>
#include <map>

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
	double intCosReco = h->Integral(20,80);
	double intCosGen  = h_gen->Integral(20,80);

  cout << "intCosReco = " << intCosReco << ", intCosGen = " << intCosGen << "\n";

  h->Scale( intCosGen / intCosReco );
}

TH1F* plotEfficiency(TH1F *h_num, TH1F *h_denom)
{
  gStyle->SetOptStat(0);

  TH1F *h_eff = (TH1F*) h_num->Clone();
  h_eff->Divide(h_denom);

  return h_eff;
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

void efficiency_cos()
{
  TGaxis::SetMaxDigits(3);

  TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR." + prod_mode + ".KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.check1.hists.all.root","READ");

  vector<TString> gen_reco  = {"gen","reco"};
  vector<TString> PFO_mode  = {"K","Pi"};
  vector<TString> heff_name = {"momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
  unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_cos_eff;  // [GenReco][LPFO][hist]

  try
  {
    if (!file->IsOpen()) return;
    TString dir_name = "efficiency/";

    for ( auto igenreco : gen_reco ){
      for ( auto i_lmode : PFO_mode ){
        for ( auto ih : heff_name ){
          TString hname = "h_" + igenreco + "_" + i_lmode + "_" + ih;
          TH1F *h = (TH1F*) file->Get(dir_name + hname);
          h1_cos_eff[igenreco][i_lmode][ih] = h;
        }
      }
    }


    int count = 0;
    TH1F * h_denom;
    TCanvas *c_eff_gen_Pi = new TCanvas("c_eff_gen_Pi", "c_eff_gen_Pi", 1500,400);
    c_eff_gen_Pi->Divide(5,1);
    TCanvas *c_cos_gen_Pi = new TCanvas("c_cos_gen_Pi", "c_cos_gen_Pi", 1500,400);
    c_cos_gen_Pi->Divide(5,1);

    cout << "=== Gen Pi ===\n";
    
    for ( auto ih : heff_name ){

      TH1F * h_num = h1_cos_eff["gen"]["Pi"][ih];
      
      if (count) {
        
        TH1F *h_eff = plotEfficiency(h_num, h_denom);
        
        c_eff_gen_Pi->cd(count);
        h_eff->Draw("h");

      }
      
      h_denom = h_num;
      count++;

    }


  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}