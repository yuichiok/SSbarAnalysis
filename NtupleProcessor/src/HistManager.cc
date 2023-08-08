#include <iostream>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h> 

#include "HistManager.hh"

using std::cout;   using std::endl;

namespace QQbarAnalysis
{
  HistManager::HistManager() {
    bins_dEdx[0]=0.1;
    for(int i=1;i<200;i++) bins_dEdx[i]=bins_dEdx[i-1]+0.1/100.;
  }

  void HistManager::InitializeHists()
  {
    Float_t bins_p[]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80,90,100};
    Int_t   nbins_p=sizeof(bins_p)/sizeof(Float_t) - 1;

    //////////////////
    //     TH1F     //
    //////////////////

    // Gen  histograms
    for( auto iname : hcos_gen_name){
      TString hname = "h_" + iname;
      h1_gen_cos[iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
    }

    // Reco histograms
    for( auto i_lmode : _pt.PFO_mode ){
      for( auto iname : hcos_name){
        TString hname = "h_" + i_lmode + "_" + iname;
        h1_cos[i_lmode][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
      }
    }

    for( auto i_lmode : _pt.PFO_mode ){
      for( auto iname : hres_name){
        TString hname = "h_" + i_lmode + "_" + iname;
        h1_resolution[i_lmode][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
      }
    }

    // effciency
    for( auto i_gen_reco : gen_reco ){
      for( auto i_lmode : _pt.PFO_mode ){
        for( auto iname : heff_name){
          TString hname = "h_" + i_gen_reco + "_" + i_lmode + "_" + iname;
          h1_cos_eff[i_gen_reco][i_lmode][iname]  = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
        }
      }
    }

    //////////////////
    //     TH2F     //
    //////////////////

    // dEdx
      for( auto i_lmode : _pt.PFO_mode ){
        for( auto i_type : _pt.PFO_type ){
          h2_dEdx[i_lmode][i_type]["dEdx_p"]        = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_p",";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
          h2_dEdx[i_lmode][i_type]["dEdx_dist_cos"] = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_dist_cos",";cos#theta;#frac{dE}{dx} distance",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,dEdx_dist_min,dEdx_dist_max);
        }
      }

      Hist2List();

  }

  void HistManager::Hist2List()
  {
    // gen
    for ( const auto &[iname, hist] : h1_gen_cos ){
      hList1_gen_cos->Add(hist);
    }

    // reco
    for ( const auto &[i_lmode, hists] : h1_cos ){
      for ( const auto &[icut, hist] : hists ){
        hList1_cos->Add(hist);
      }
    }

    for ( const auto &[i_lmode, hists] : h1_resolution ){
      for ( const auto &[icut, hist] : hists ){
        hList1_resolution->Add(hist);
      }
    }

    // effiency
    for ( const auto &[i_gen_reco, type_hists] : h1_cos_eff ){
      for ( const auto &[i_lmode, hists] : type_hists ){
        for ( const auto &[icut, hist] : hists ){
          hList1_efficiency->Add(hist);
        }
      }
    }

    // 2D hist
    for ( const auto &[i_lmode, type_hists] : h2_dEdx ){
      for ( const auto &[i_type, hists] : type_hists ){
        for ( const auto &[icut, hist] : hists ){
          hList2_dEdx->Add(hist);
        }
      }
    }

  }

  void HistManager::WriteLists( TFile * output)
  {
    // Focus to this file
    output->cd();

    // Use unordered_map
    TDirectory * d_gen = output->mkdir("gen");
      d_gen->cd();
      hList1_gen_cos->Write();
      output->cd();

    TDirectory * d_cos = output->mkdir("cos");
      d_cos->cd();
      hList1_cos->Write();
      output->cd();

    TDirectory * d_resolution = output->mkdir("resolution");
      d_resolution->cd();
      hList1_resolution->Write();
      output->cd();

    TDirectory * d_efficiency = output->mkdir("efficiency");
      d_efficiency->cd();
      hList1_efficiency->Write();
      output->cd();

    TDirectory * d_dEdx = output->mkdir("dEdx");
      d_dEdx->cd();
      hList2_dEdx->Write();
      output->cd();

  }
}