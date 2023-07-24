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

      unordered_map< TString, unordered_map< TString, TH1F* > > h1_cos;        // [LPFO][hist]
      unordered_map< TString, unordered_map< TString, TH1F* > > h1_resolution; // [LPFO][hist]

    // h2 hist
      enum h2_index {

        gen_K_p_cos,
        reco_K_p_cos,

        stable_K_cos,
        purity_K_cos,

        stable_Pi_cos,
        purity_Pi_cos,

        dummy_h2,
        Last_h2 = dummy_h2
      };
      TH2F * h2[Last_h2];

      unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]

    private:
    // Lists
      TList* hList1                = new TList();
      TList* hList2                = new TList();

      TList* hList1_cos            = new TList();
      TList* hList1_resolution     = new TList();
      TList* hList2_dEdx           = new TList();

    // PFO Tools
      PFOTools _pt;

  };
}
#endif