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

        reco_K_cos,
        reco_K_qcos,
        reco_K_scos,
        reco_K_mom,
        gen_reco_K_sep_cos,
        jet_reco_K_sep_cos,

        reco_Pi_cos,
        reco_Pi_qcos,
        reco_Pi_scos,
        reco_Pi_mom,
        gen_reco_Pi_sep_cos,
        jet_reco_Pi_sep_cos,

        reco_sum_jetE,
        reco_jet_sep,

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

      unordered_map< TString, unordered_map< TString, TH1F* > > h1_cos; // [LPFO][hist]

      enum h1_pq_index {
        // KID
        acc_KK,
        rej_KK,

        // PiID
        acc_PiPi,
        rej_PiPi,

        dummy_pq,
        Last_h1_pq = dummy_pq
      };
      TH1F * h1_pq[Last_h1_pq];

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

      enum h2_jet_index {

        jet_mult_cos,
        jet_mult_cos_noISR,

        dummy_h2_jet,
        Last_dummy_h2_jet = dummy_h2_jet
      };
      TH2F * h2_jet[Last_dummy_h2_jet];

      unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]

    private:
    // Lists
      TList* hList1                = new TList();
      TList* hList1_pq             = new TList();
      TList* hList1_particle_ratio = new TList();
      TList* hList2                = new TList();
      TList* hList2_jet            = new TList();

      TList* hList1_cos            = new TList();
      TList* hList2_dEdx           = new TList();

    // PFO Tools
      PFOTools _pt;

  };
}
#endif