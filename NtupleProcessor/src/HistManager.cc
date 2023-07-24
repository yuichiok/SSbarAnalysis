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
  HistManager::HistManager() {}

  void HistManager::InitializeHists()
  {
    Int_t   nbins_cos = 100;
    Float_t cos_min = -1.0, cos_max = 1.0;
    Float_t bins_cos_fine[] = {-1.0,-0.98,-0.96,-0.94,-0.92,-0.90,-0.88,-0.86,-0.84,-0.82,-0.80,-0.75,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,
                              0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.0};
    Int_t   nbins_cos_fine = sizeof(bins_cos_fine) / sizeof(Float_t) - 1;

    Float_t bins_p[]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80,90,100};
    Int_t   nbins_p=sizeof(bins_p)/sizeof(Float_t) - 1;

    Float_t bins_dEdx[200];
    bins_dEdx[0]=0.1;
    for(int i=1;i<200;i++) bins_dEdx[i]=bins_dEdx[i-1]+0.1/100.;
    Int_t nbins_dEdx=199;

    Float_t nbins_dEdx_dist = 100;
    Float_t dEdx_dist_min = -20.0, dEdx_dist_max = 20.0;

    //////////////////
    //     TH1F     //
    //////////////////

    // Mapping
    vector<TString> hcos_name = {"cos","qcos","scos","acc_cos","rej_cos"};
    for( auto i_lmode : _pt.PFO_mode ){
      for( auto iname : hcos_name){
        TString hname = "h_" + i_lmode + "_" + iname;
        h1_cos[i_lmode][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
      }
    }

    vector<TString> hres_name = {"gen_N_cos","reco_N_cos","N_corr_cos"};
    for( auto i_lmode : _pt.PFO_mode ){
      for( auto iname : hres_name){
        TString hname = "h_" + i_lmode + "_" + iname;
        h1_resolution[i_lmode][iname] = new TH1F(hname,iname + ";cos#theta;Entries",nbins_cos,cos_min,cos_max);
      }
    }



    // Original

      h1[gen_q_cos]       = new TH1F("h_gen_q_cos","; Generated q#bar{q} cos#theta; Entries",nbins_cos,-1,1);
      h1[gen_q_qcos]      = new TH1F("h_gen_q_qcos","; Generated q#bar{q} qcos#theta; Entries",nbins_cos,-1,1);

      h1[gen_K_cos]       = new TH1F("h_gen_K_cos","; Generated Kaon cos#theta; Entries",nbins_cos,-1,1);
      h1[gen_K_qcos]      = new TH1F("h_gen_K_qcos","; Generated Kaon qcos#theta; Entries",nbins_cos,-1,1);
      
      h1[reco_K_cos]      = new TH1F("h_reco_K_cos",";LPFO Kaon cos#theta; Entries",nbins_cos,-1,1);
      h1[reco_K_qcos]     = new TH1F("h_reco_K_qcos",";LPFO Kaon cos#theta; Entries",nbins_cos,-1,1);
      h1[reco_K_scos]     = new TH1F("h_reco_K_scos",";LPFO Kaon cos#theta; Entries",nbins_cos,-1,1);
      h1[reco_K_mom]      = new TH1F("h_reco_K_mom",";LPFO Kaon momentum (GeV); Entries",140,10,150);
      h1[gen_reco_K_sep_cos] = new TH1F("h_gen_reco_K_sep_cos",";Kaon cos#theta_{gen-LPFOK}; Entries",nbins_cos,-1,1);

      h1[reco_Pi_cos]      = new TH1F("h_reco_Pi_cos",";LPFO Pion cos#theta; Entries",nbins_cos,-1,1);
      h1[reco_Pi_qcos]     = new TH1F("h_reco_Pi_qcos",";LPFO Pion cos#theta; Entries",nbins_cos,-1,1);
      h1[reco_Pi_scos]     = new TH1F("h_reco_Pi_scos",";LPFO Pion cos#theta; Entries",nbins_cos,-1,1);
      h1[reco_Pi_mom]      = new TH1F("h_reco_Pi_mom",";LPFO Pion momentum (GeV); Entries",140,10,150);
      h1[gen_reco_Pi_sep_cos] = new TH1F("h_gen_reco_Pi_sep_cos",";cos#theta_{gen-LPFOPi}; Entries",nbins_cos,-1,1);

    // Number of Gen Reco Kaons & Pions
    /*
    SP_count [gen,reco,correct][K,Pi]
    */

      h1[gen_N_K_cos]     = new TH1F("h_gen_N_K_cos",";cos#theta;Generated N Kaons", nbins_cos,-1,1);
      h1[reco_N_K_cos]    = new TH1F("h_reco_N_K_cos",";cos#theta;Reconstructed N Kaons", nbins_cos,-1,1);
      h1[N_K_corr_cos]    = new TH1F("h_N_K_corr_cos",";cos#theta;Correctly Reconstructed N Kaons", nbins_cos,-1,1);

      h1[gen_N_Pi_cos]     = new TH1F("h_gen_N_Pi_cos",";cos#theta;Generated N Pions", nbins_cos,-1,1);
      h1[reco_N_Pi_cos]    = new TH1F("h_reco_N_Pi_cos",";cos#theta;Reconstructed N Pions", nbins_cos,-1,1);
      h1[N_Pi_corr_cos]    = new TH1F("h_N_Pi_corr_cos",";cos#theta;Correctly Reconstructed N Pions", nbins_cos,-1,1);

    //////////////////
    //     TH2F     //
    //////////////////
    /*
    mom_cos [gen,reco][K,Pi]
    count_gen_reco [K,Pi]
    SP        [K,Pi][stability,purity]
    ratio_cos [K,pi,p,e,mu]
    dEdx_p        [K,Pi][K,Pi,p,e,mu]
    dEdx_dist_cos  [K,Pi][K,Pi][K,Pi,p,e,mu]
    */

      h2[gen_K_p_cos]  = new TH2F("h2_gen_K_p_cos" ,";cos#theta;p (GeV)",nbins_cos,-1,1,100,0,60);
      h2[reco_K_p_cos] = new TH2F("h2_reco_K_p_cos",";cos#theta;p (GeV)",nbins_cos,-1,1,100,0,60);
      
      h2[stable_K_cos]   = new TH2F("h2_stable_K_cos",";cos#theta;Stability",nbins_cos,-1,1,50,0,1);
      h2[purity_K_cos]   = new TH2F("h2_purity_K_cos",";cos#theta;Purity",   nbins_cos,-1,1,50,0,1);

      h2[stable_Pi_cos]   = new TH2F("h2_stable_Pi_cos",";cos#theta;Stability",nbins_cos,-1,1,50,0,1);
      h2[purity_Pi_cos]   = new TH2F("h2_purity_Pi_cos",";cos#theta;Purity",   nbins_cos,-1,1,50,0,1);

    // dEdx
      vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
      for( auto i_lmode : _pt.PFO_mode ){
        for( auto i_type : _pt.PFO_type ){
          h2_dEdx[i_lmode][i_type]["dEdx_p"]        = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_p",";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
          h2_dEdx[i_lmode][i_type]["dEdx_dist_cos"] = new TH2F("h2_" + i_lmode + "_" + i_type + "_dEdx_dist_cos",";cos#theta;#frac{dE}{dx}",nbins_cos,cos_min,cos_max,nbins_dEdx_dist,dEdx_dist_min,dEdx_dist_max);
        }
      }

      Hist2List();

  }

  void HistManager::Hist2List()
  {
    for (auto ih : h1) {
      hList1->Add(ih);
    }

    for (auto ih : h2) {
      hList2->Add(ih);
    }
    
    // Use unordered_map
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

    // Write histogram lists
    hList1->Write();
    hList2->Write();

    // Use unordered_map
    TDirectory * d_cos = output->mkdir("cos");
      d_cos->cd();
      hList1_cos->Write();
      output->cd();

    TDirectory * d_resolution = output->mkdir("resolution");
      d_resolution->cd();
      hList1_resolution->Write();
      output->cd();

    TDirectory * d_dEdx = output->mkdir("dEdx");
      d_dEdx->cd();
      hList2_dEdx->Write();
      output->cd();

  }
}