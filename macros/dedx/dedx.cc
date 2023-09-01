#include <iostream>
#include <vector>
#include <map>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

const TString prod_mode = "ud";
const TString LPFO_mode = "Pi";

const vector<TString> PFO_mode  = {"K","Pi"};
const vector<TString> PFO_type  = {"K","Pi", "p", "e", "mu"};

TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR." + prod_mode + ".KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root","READ");

TH1F* plotEfficiency(TH1F *h_num, TH1F *h_denom)
{
  gStyle->SetOptStat(0);

  TH1F *h_eff = (TH1F*) h_num->Clone();
  h_eff->Divide(h_denom);
  h_eff->GetYaxis()->SetRangeUser(0,1);

  return h_eff;
}

void dedxDistCosProjType()
{
  vector<TString> gen_reco  = {"gen","reco"};
  vector<TString> cut_name = {"nocut","momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
  vector<TString> heff_dedx_name = {"dEdx_p","dEdx_cos","dEdx_error_cos","dEdx_dist_cos"};
  unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > > h2_dEdx_eff;  // [GenReco][LPFO][TruthID][cut][hist]
  TString dir_name = "efficiency/dEdx/";

  for ( auto igen_reco : gen_reco ){
    for ( auto imode : PFO_mode ){
      for ( auto itype : PFO_type ){
        for ( auto icut : cut_name ){
          for ( auto ihist : heff_dedx_name ){
            TString hname = "h2_" + igen_reco + "_" + imode + "_" + itype + "_" + icut + "_" + ihist;
            TH2F *h = (TH2F*) file->Get(dir_name + hname);
            h2_dEdx_eff[igen_reco][imode][itype][icut][ihist] = h;
          }
        }
      }
    }
  }

  TCanvas *c_type_dedx_cos_type = new TCanvas("c_type_dedx_cos_type", "c_type_dedx_cos_type", 1500,1000);
  c_type_dedx_cos_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_dist_cos_type = new TCanvas("c_type_dedx_dist_cos_type", "c_type_dedx_dist_cos_type", 1500,1000);
  c_type_dedx_dist_cos_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_dist_cos_proj_type = new TCanvas("c_type_dedx_dist_cos_proj_type", "c_type_dedx_dist_cos_proj_type", 1500,1000);
  c_type_dedx_dist_cos_proj_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_dist_cos_proj_efficiency_type = new TCanvas("c_type_dedx_dist_cos_proj_efficiency_type", "c_type_dedx_dist_cos_proj_efficiency_type", 1500,1000);
  c_type_dedx_dist_cos_proj_efficiency_type->Divide(cut_name.size(),PFO_type.size());

  TH1F *PiID_eff;

  int count = 0;
  for ( auto itype : PFO_type ){
    TH1F *h_denom = nullptr;
    for ( auto iname : cut_name ){
      
      count++;

      TString hTitle = itype + " | " + iname;
      TH2F *h_dedx_cos           = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_cos");
      TH2F *h_dedx_dist_cos      = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_dist_cos");
      TH1F *h_dedx_dist_cos_proj = (TH1F*) h_dedx_dist_cos->ProjectionX();

      h_dedx_dist_cos->SetTitle(hTitle);
      h_dedx_dist_cos_proj->SetTitle(hTitle);

      c_type_dedx_cos_type->cd(count);
      h_dedx_cos->Draw("colz");
      h_dedx_cos->Draw("cont3 same");

      c_type_dedx_dist_cos_type->cd(count);
      h_dedx_dist_cos->Draw("colz");
      h_dedx_dist_cos->Draw("cont3 same");
      
      c_type_dedx_dist_cos_proj_type->cd(count);
      StyleHist(h_dedx_dist_cos_proj,kBlue);
      h_dedx_dist_cos_proj->Draw("h");

      TH1F *h_efficiency;
      if(h_denom){
        h_efficiency = plotEfficiency(h_dedx_dist_cos_proj, h_denom);
        h_denom = (TH1F*) h_dedx_dist_cos_proj->Clone();
        c_type_dedx_dist_cos_proj_efficiency_type->cd(count);
        StyleHist(h_efficiency,kBlue);
        h_efficiency->Draw("h");
        if(itype == "Pi" && iname == "PID"){
          PiID_eff = (TH1F*) h_efficiency->Clone();
        }
      }else{
        h_denom = (TH1F*) h_dedx_dist_cos_proj->Clone();        
      }

    }

  }

}

void dedxDistCosProj()
{
  vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
  unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]
  TString dir_name = "dEdx/";

  for ( auto imode : PFO_mode ){
    for ( auto itype : PFO_type ){
      for ( auto icut : hdEdx_name ){
        TString hname = "h2_" + imode + "_" + itype + "_" + icut;
        TH2F *h = (TH2F*) file->Get(dir_name + hname);
        h2_dEdx[imode][itype][icut] = h;
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

void dedxCos()
{
  vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
  unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]
  TString dir_name = "dEdx/";

  for ( auto imode : PFO_mode ){
    for ( auto itype : PFO_type ){
      for ( auto icut : hdEdx_name ){
        TString hname = "h2_" + imode + "_" + itype + "_" + icut;
        TH2F *h = (TH2F*) file->Get(dir_name + hname);
        h2_dEdx[imode][itype][icut] = h;
      }
    }
  }

  TCanvas *c_type_dedx_cos = new TCanvas("c_type_dedx_cos", "c_type_dedx_cos", 1500,400);
  c_type_dedx_cos->Divide(5,1);

  int count = 0;
  for ( auto itype : PFO_type ){

    count++;

    TH2F *h_dedx_cos = h2_dEdx.at(LPFO_mode).at(itype).at("dEdx_cos");

    c_type_dedx_cos->cd(count);
    h_dedx_cos->Draw("colz");

  }

}

void dedx()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  if (!file->IsOpen()) return;

  // dedxDistCosProj();
  dedxDistCosProjType();
  // dedxCos();


}