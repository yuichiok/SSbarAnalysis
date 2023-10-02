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

    for( auto iqq : QQ_mode ){

      InitLists( iqq );

      TString hname_prefix  = "h_" + iqq + "_";
      TString hname_prefix2 = "h2_" + iqq + "_";

      // Gen  histograms
      for( auto iname : hcos_gen_name){
        TString hname = hname_prefix + iname;
        h1_gen_cos[iqq][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
      }

      // Preselection  histograms
      h1_preselection[iqq]["sinacol"] = new TH1F(hname_prefix + "sinacol","sinacol;sin#Psi_{acol};Entries",100,0,1);
      h1_preselection[iqq]["invM"]    = new TH1F(hname_prefix + "invM","invM;m_{j_{1},j_{2}} [GeV];Entries",100,0,500);
      h1_preselection[iqq]["y23"]     = new TH1F(hname_prefix + "y23","y23;y_{23};Entries",50,0,0.25);
      h1_preselection[iqq]["LPFOacol"] = new TH1F(hname_prefix + "LPFOacol","LPFOacol;cos#theta_{L_{1},L_{2}};Entries",100,0,1);
      h1_preselection[iqq]["cosBF"]   = new TH1F(hname_prefix + "cosBF","cosBF;cos#theta;Entries",nbins_cos,cos_min,cos_max);
      h1_preselection[iqq]["cosAF"]   = new TH1F(hname_prefix + "cosAF","cosAF;cos#theta;Entries",nbins_cos,cos_min,cos_max);

      // Reco histograms
      for( auto i_lmode : _pt.PFO_mode ){
        for( auto iname : hcos_name){
          TString hname = hname_prefix + i_lmode + "_" + iname;
          h1_cos[iqq][i_lmode][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
        }
      }

      for( auto i_lmode : _pt.PFO_mode ){
        for( auto iname : hres_name){
          TString hname = hname_prefix + i_lmode + "_" + iname;
          h1_resolution[iqq][i_lmode][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
        }
      }

      // effciency
      for( auto i_gen_reco : gen_reco ){
        for( auto i_lmode : _pt.PFO_mode ){
          for( auto iname : heff_name){
            TString hname = hname_prefix + i_gen_reco + "_" + i_lmode + "_" + iname;
            h1_cos_eff[iqq][i_gen_reco][i_lmode][iname]  = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
            for( auto i_type : _pt.PFO_type ){
              TString hname_dedx_dist = hname_prefix2 + i_gen_reco + "_" + i_lmode + "_" + i_type + "_" + iname;
              h2_dEdx_eff[iqq][i_gen_reco][i_lmode][i_type][iname]["dEdx_dist_cos"]  = new TH2F(hname_dedx_dist + "_dEdx_dist_cos",";cos#theta;#frac{dE}{dx} distance",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,dEdx_dist_min,dEdx_dist_max);
              h2_dEdx_eff[iqq][i_gen_reco][i_lmode][i_type][iname]["dEdx_error_cos"] = new TH2F(hname_dedx_dist + "_dEdx_error_cos",";cos#theta;#frac{dE}{dx} error",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,0.0001,0.005);
              h2_dEdx_eff[iqq][i_gen_reco][i_lmode][i_type][iname]["dEdx_cos"]       = new TH2F(hname_dedx_dist + "_dEdx_cos",";cos#theta;#frac{dE}{dx}",nbins_cos,bins_cos,nbins_dEdx,bins_dEdx);
              h2_dEdx_eff[iqq][i_gen_reco][i_lmode][i_type][iname]["dEdx_p"]         = new TH2F(hname_dedx_dist + "_dEdx_p",";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
            }
          }
        }
      }

      // dEdx
      for( auto i_lmode : _pt.PFO_mode ){
        for( auto i_type : _pt.PFO_type ){
          h2_dEdx[iqq][i_lmode][i_type]["dEdx_p"]        = new TH2F(hname_prefix2 + i_lmode + "_" + i_type + "_dEdx_p",";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
          h2_dEdx[iqq][i_lmode][i_type]["dEdx_cos"]      = new TH2F(hname_prefix2 + i_lmode + "_" + i_type + "_dEdx_cos",";cos#theta;#frac{dE}{dx}",nbins_cos,bins_cos,nbins_dEdx,bins_dEdx);
          h2_dEdx[iqq][i_lmode][i_type]["dEdx_dist_cos"] = new TH2F(hname_prefix2 + i_lmode + "_" + i_type + "_dEdx_dist_cos",";cos#theta;#frac{dE}{dx} distance",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,dEdx_dist_min,dEdx_dist_max);
        }
      }

      CreateLists( iqq );

    } // end of QQ_mode loop


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

  void HistManager::InitLists( TString iqq )
  {
    // gen
    hList1[iqq]["gen_cos"] = new TList();
    // preselection
    hList1[iqq]["preselection"] = new TList();
    // reco
    hList1[iqq]["cos"] = new TList();
    hList1[iqq]["resolution"] = new TList();
    // efficiency
    hList1[iqq]["efficiency"] = new TList();
    hList2[iqq]["efficiency"] = new TList();
    // dEdx
    hList2[iqq]["dEdx"] = new TList();

  }

  void HistManager::CreateLists( TString iqq )
  {
    // gen
    recursiveIterate( hList1[iqq]["gen_cos"], h1_gen_cos[iqq] );
    // preselection
    recursiveIterate( hList1[iqq]["preselection"], h1_preselection[iqq] );
    // reco
    recursiveIterate( hList1[iqq]["cos"], h1_cos[iqq] );
    recursiveIterate( hList1[iqq]["resolution"], h1_resolution[iqq] );
    // efficiency
    recursiveIterate( hList1[iqq]["efficiency"], h1_cos_eff[iqq] );
    recursiveIterate( hList2[iqq]["efficiency"], h2_dEdx_eff[iqq] );
    // dEdx
    recursiveIterate( hList2[iqq]["dEdx"], h2_dEdx[iqq] );

  }

  void HistManager::CreateDirectories( TFile * output, const TString& iqq )
  {
    // Focus to this file
    TDirectory* baseDirectory = output->mkdir(iqq);
    baseDirectory->cd();

    // Use unordered_map
    TDirectory * d_gen = baseDirectory->mkdir("gen");
      d_gen->cd();
      hList1.at(iqq).at("gen_cos")->Write();
      baseDirectory->cd();

    TDirectory * d_preselection = baseDirectory->mkdir("preselection");
      d_preselection->cd();
      hList1.at(iqq).at("preselection")->Write();
      baseDirectory->cd();

    TDirectory * d_cos = baseDirectory->mkdir("cos");
      d_cos->cd();
      hList1.at(iqq).at("cos")->Write();
      baseDirectory->cd();

    TDirectory * d_resolution = baseDirectory->mkdir("resolution");
      d_resolution->cd();
      hList1.at(iqq).at("resolution")->Write();
      baseDirectory->cd();

    TDirectory * d_efficiency = baseDirectory->mkdir("efficiency");
      d_efficiency->cd();
      hList1.at(iqq).at("efficiency")->Write();
      TDirectory * d_efficiency_dEdx_dist_cos = d_efficiency->mkdir("dEdx");
        d_efficiency_dEdx_dist_cos->cd();
        hList2.at(iqq).at("efficiency")->Write();
        baseDirectory->cd();

    TDirectory * d_dEdx = baseDirectory->mkdir("dEdx");
      d_dEdx->cd();
      hList2.at(iqq).at("dEdx")->Write();
      baseDirectory->cd();

    output->cd();

  }

  void HistManager::WriteLists( TFile * output)
  {
    // Focus to this file
    output->cd();

    vector<TString> lists = {"gen_cos","cos","resolution","efficiency","dEdx"};

    // qq or bg selection
    for ( auto iqq : QQ_mode ){

      CreateDirectories( output, iqq );

    } // end of QQ_mode loop

  }

}