#include <iostream>
#include "../include/Styles.hh"

TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.ISR.Kv35.all.root","READ");

void plot_ISR()
{
  vector<TString> h_category = {"npfos", "photon_Ecos"};
  vector<TString> sig_isr    = {"ISR","signal"};
  map<TString, map<TString, TH2F*>> h2d;

  for(auto& cat : h_category){
    for(auto& sig : sig_isr){
      TString dir   = "ISR_Analysis/";
      TString hname = "h2_" + cat + "_" + sig;
      h2d[cat][sig] = (TH2F*)file->Get(dir + hname);
    }
  }

}

void ISR_plots()
{
  gStyle->SetOptStat(0);
  plot_ISR();


}