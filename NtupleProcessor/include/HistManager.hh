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
    const Int_t hmax = 100;
    enum h_index {
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
      Last = lpfo_reco_K_mom
    };

    TH1F * h[Last];

  private:
    TList* hList = new TList();

};

#endif