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
    bins_cos[0]=-1.0;
    for(int i=1;i<nbins_cos+1;i++) bins_cos[i]=bins_cos[i-1]+2./nbins_cos;
  }

  void HistManager::InitializeHists()
  {
    Float_t bins_p[]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80,90,100};
    Int_t   nbins_p=sizeof(bins_p)/sizeof(Float_t) - 1;

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
          for( auto i_type : _pt.PFO_type ){
            TString hname_dedx_dist = "h2_" + i_gen_reco + "_" + i_lmode + "_" + i_type + "_" + iname;
            h2_dEdx_eff[i_gen_reco][i_lmode][i_type][iname]["dEdx_dist_cos"]  = new TH2F(hname_dedx_dist + "_dEdx_dist_cos",";cos#theta;#frac{dE}{dx} distance",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,dEdx_dist_min,dEdx_dist_max);
            h2_dEdx_eff[i_gen_reco][i_lmode][i_type][iname]["dEdx_error_cos"] = new TH2F(hname_dedx_dist + "_dEdx_error_cos",";cos#theta;#frac{dE}{dx} error",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,0.0001,0.005);
            h2_dEdx_eff[i_gen_reco][i_lmode][i_type][iname]["dEdx_cos"]       = new TH2F(hname_dedx_dist + "_dEdx_cos",";cos#theta;#frac{dE}{dx}",nbins_cos,bins_cos,nbins_dEdx,bins_dEdx);
            h2_dEdx_eff[i_gen_reco][i_lmode][i_type][iname]["dEdx_p"]         = new TH2F(hname_dedx_dist + "_dEdx_p",";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
          }
        }
      }
    }

    // dEdx
    for( auto i_lmode : _pt.PFO_mode ){
      for( auto i_type : _pt.PFO_type ){
        h2_dEdx[i_lmode][i_type]["dEdx_p"]        = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_p",";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
        h2_dEdx[i_lmode][i_type]["dEdx_cos"]      = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_cos",";cos#theta;#frac{dE}{dx}",nbins_cos,bins_cos,nbins_dEdx,bins_dEdx);
        h2_dEdx[i_lmode][i_type]["dEdx_dist_cos"] = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_dist_cos",";cos#theta;#frac{dE}{dx} distance",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,dEdx_dist_min,dEdx_dist_max);
      }
    }

    CreateLists();

  }

  template <typename Key, typename Value>
  void HistManager::recursiveIterate(TList *list, const unordered_map<Key, Value>& map) {
    for (const auto &[key, val] : map) {
      if constexpr ( std::is_same_v<Value, TH1F*> || std::is_same_v<Value, TH2F*> ) {
        list->Add(val);
      } else if constexpr ( std::is_same_v<Value, unordered_map<Key, typename Value::mapped_type>> ) {
        recursiveIterate(list,val);
      }
    }
  }

  void HistManager::CreateLists()
  {
    // gen
    recursiveIterate( hList1_gen_cos, h1_gen_cos );
    // reco
    recursiveIterate( hList1_cos, h1_cos );
    recursiveIterate( hList1_resolution, h1_resolution );
    // efficiency
    recursiveIterate( hList1_efficiency, h1_cos_eff );
    recursiveIterate( hList2_efficiency, h2_dEdx_eff );
    // 2D hist
    recursiveIterate( hList2_dEdx, h2_dEdx );   

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
      TDirectory * d_efficiency_dEdx_dist_cos = d_efficiency->mkdir("dEdx");
        d_efficiency_dEdx_dist_cos->cd();
        hList2_efficiency->Write();
        output->cd();

    TDirectory * d_dEdx = output->mkdir("dEdx");
      d_dEdx->cd();
      hList2_dEdx->Write();
      output->cd();

  }
}