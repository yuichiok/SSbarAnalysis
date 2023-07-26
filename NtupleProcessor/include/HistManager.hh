#ifndef GUARD_HistManager_h
#define GUARD_HistManager_h

#include "PFOTools.hh"

#include <iostream>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h>

using std::unordered_map;

namespace QQbarAnalysis
{
  class HistManager
  {
    public:
      HistManager();
      ~HistManager(){}

    // Methods
      virtual void InitializeHists();
      virtual void Hist2List();
      virtual void WriteLists( TFile * output );

    // Variables
      Int_t   nbins_cos = 100;
      Float_t cos_min = -1.0, cos_max = 1.0;

      Float_t bins_dEdx[200];
      Int_t nbins_dEdx=199;

      Float_t nbins_dEdx_dist = 100;
      Float_t dEdx_dist_min = -20.0, dEdx_dist_max = 20.0;


    // Declear histograms
    // h1 hist
      enum h1_index {
        gen_q_cos,
        gen_q_qcos,

        gen_K_cos,
        gen_K_qcos,

        gen_N_K_cos,
        reco_N_K_cos,
        N_K_corr_cos,

        gen_N_Pi_cos,
        reco_N_Pi_cos,
        N_Pi_corr_cos,

        dummy_h1,
        Last_h1 = dummy_h1
      };
      TH1F * h1[Last_h1];

      vector<TString> hcos_name = {"cos","qcos","scos","acc_cos","rej_cos"};
      unordered_map< TString, unordered_map< TString, TH1F* > > h1_cos;        // [LPFO][hist]

      vector<TString> hres_name = {"gen_N_cos","reco_N_cos","N_corr_cos"};
      unordered_map< TString, unordered_map< TString, TH1F* > > h1_resolution; // [LPFO][hist]


    // h2 hist
      vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
      unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]

    private:
    // Lists
      TList* hList1                = new TList();

      TList* hList1_cos            = new TList();
      TList* hList1_resolution     = new TList();
      TList* hList2_dEdx           = new TList();

    // PFO Tools
      PFOTools _pt;

  };
}
#endif