#ifndef GUARD_HistManager_h
#define GUARD_HistManager_h

#include <iostream>
#include <map>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h> 

const int NSLABS = 15;

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

      reco_sum_jetE,
      reco_jet_sep,

      lpfo_gen_K_mom,
      lpfo_reco_K_mom,

      gen_N_K_cos,
      reco_N_K_cos,
      N_K_corr_cos,

      dummy1,
      LastH1 = dummy1
    };
    TH1F * h1[LastH1];

    enum h1_pq_index {
      acc_KK,
      rej_KK,

      dummy2,
      LastH2 = dummy2
    };
    TH1F * h1_pq[LastH2];

  // h2 hist
    enum h2_index {

      stable_cos,
      purity_cos,

      dummy3,
      LastH3 = dummy3
    };
    TH2F * h2[LastH3];


  private:
    TList* hList1    = new TList();
    TList* hList1_pq = new TList();
    TList* hList2    = new TList();

};

#endif