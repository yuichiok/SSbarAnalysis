#include <iostream>
#include <vector>
#include <map>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString prod_mode = "uu";
TString LPFO_mode = "Pi";

void dedx()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR." + prod_mode + ".KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.check2.hists.all.root","READ");

  vector<TString> PFO_mode  = {"K","Pi"};
  vector<TString> PFO_type  = {"K","Pi", "p", "e", "mu"};
  vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
  unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]

  try
  {
    if (!file->IsOpen()) return;
    TString dir_name = "dEdx/";

    for ( auto imode : PFO_mode ){
      for ( auto itype : PFO_type ){
        for ( auto ih : hdEdx_name ){
          TString hname = "h2_" + imode + "_" + itype + "_" + ih;
          TH2F *h = (TH2F*) file->Get(dir_name + hname);
          h2_dEdx[imode][itype][ih] = h;
        }
      }
    }


    TCanvas *c_dedx_dist_cos = new TCanvas("c_dedx_dist_cos", "c_dedx_dist_cos", 500,500);
    TCanvas *c_type_dedx_dist_cos = new TCanvas("c_type_dedx_dist_cos", "c_type_dedx_dist_cos", 1500,400);
    c_type_dedx_dist_cos->Divide(5,1);

    int count = 0;
    TH1F * h_dedx_dist_cos_proj_sum = new TH1F("h_dedx_dist_cos_proj_sum","h_dedx_dist_cos_proj_sum",100,-1,1);
    StyleHist(h_dedx_dist_cos_proj_sum,kBlue);
    for ( auto itype : PFO_type ){
      TH2F *h_dedx_dist_cos = h2_dEdx.at(LPFO_mode).at(itype).at("dEdx_dist_cos");
      TH1F *h_dedx_dist_cos_proj = (TH1F*) h_dedx_dist_cos->ProjectionX();
      StyleHist(h_dedx_dist_cos_proj,kBlue);

      h_dedx_dist_cos_proj_sum->Add(h_dedx_dist_cos_proj);

      c_type_dedx_dist_cos->cd(++count);
      h_dedx_dist_cos_proj->Draw("h");
    }

    c_dedx_dist_cos->cd();
    h_dedx_dist_cos_proj_sum->Draw("h");

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}