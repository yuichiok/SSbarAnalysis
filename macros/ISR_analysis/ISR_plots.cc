#include <iostream>
#include "../include/Styles.hh"

TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ud.ISR.Kv35.all.root","READ");

TH2F* h2_ISR_signal(map<TString, map<TString, TH2F*>> hmap, TString category)
{
  TH2F* h2 = (TH2F*)hmap[category]["ISR"]->Clone();
  h2->Divide(hmap[category]["signal"]);
  return h2;
}

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
      if(cat == "npfos") h2d[cat][sig]->SetTitle(sig);
      if(cat == "photon_Ecos") h2d[cat][sig]->SetTitle(";|cos#theta_{#gamma_{clus}}|;E_{#gamma_{clus}} (GeV)");
    }
  }

  for(auto& cat : h_category){
    h2d[cat]["ratio"] = h2_ISR_signal(h2d, cat);
    // if(cat == "npfos") h2d[cat]["ratio"]->SetTitle("ISR/signal");
  }
  sig_isr.push_back("ratio");

  Int_t ncan = h_category.size() * sig_isr.size();
  TCanvas* cs[ncan];
  Int_t ccnt = 0;
  for( int ic = 0; ic < ncan; ic++ ) cs[ic] = new TCanvas("c" + TString::Format("%d",ic), "c" + TString::Format("%d",ic),800,800);

  Int_t ic = 0;
  for(auto& cat : h_category){
    for(auto& sig : sig_isr){
      cs[ic]->cd();
      StylePad(gPad,0,0,0.14,0.14);
      gPad->SetGrid(0,0);
      h2d[cat][sig]->SetTitle("");
      h2d[cat][sig]->Draw("colz");
      ic++;
    }
  }



  // TCanvas* c1 = new TCanvas("c1","c1",1000,700);
  // c1->Divide(3,2);
  // Int_t ic=1;
  // for(auto& cat : h_category){
  //   for(auto& sig : sig_isr){
  //     c1->cd(ic);
  //     StylePad(gPad,0,0,0.14,0.14);
  //     gPad->SetGrid(0,0);
  //     h2d[cat][sig]->Draw("colz");
  //     ic++;
  //   }
  // }

}

void ISR_plots()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  plot_ISR();


}