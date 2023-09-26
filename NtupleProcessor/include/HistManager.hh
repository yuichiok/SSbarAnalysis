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
      template <typename Key, typename Value>
      void recursiveIterate(TList *list, const unordered_map<Key, Value>& map);
      virtual void InitLists( TString pmode );
      virtual void CreateLists( TString pmode );
      virtual void CreateDirectories( TFile * output, const TString& baseDir );
      virtual void WriteLists( TFile * output );

    // Variables
      const static Int_t   nbins_cos = 100;
      static constexpr Float_t cos_min = -1.0, cos_max = 1.0;

      Float_t bins_dEdx[200];
      const static Int_t nbins_dEdx=199;

      Float_t bins_cos[nbins_cos+1];

      static constexpr Float_t nbins_dEdx_dist = 100;
      static constexpr Float_t dEdx_dist_min = -20.0, dEdx_dist_max = 20.0;


    // Declear histograms

    // process declaration
      vector<TString> QQ_mode = {"dd","uu","ss","cc","bb","bg","rr"};

    // h1 hist
      vector<TString> hcos_gen_name = {"cos","qcos"};
      unordered_map<TString, unordered_map< TString, TH1F* > > h1_gen_cos;        // [QQ_mode][hist]

      vector<TString> hpreselection_name = {"sinacol","invM","y23","cosBF","cosAF"};
      unordered_map<TString, unordered_map< TString, TH1F* > > h1_preselection;        // [QQ_mode][hist]

      vector<TString> hcos_name = {"cos","qcos","scos","acc_cos","rej_cos"};
      unordered_map<TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_cos;        // [QQ_mode][LPFO][hist]

      vector<TString> hres_name = {"gen_N_cos","reco_N_cos","N_corr_cos"};
      unordered_map<TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_resolution; // [QQ_mode][LPFO][hist]

      // efficiency plots
      vector<TString> gen_reco  = {"gen","reco"};
      vector<TString> heff_name = {"nocut","momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
      vector<TString> heff_dedx_name = {"dEdx_p","dEdx_cos","dEdx_error_cos","dEdx_dist_cos"};
      unordered_map<TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > > h1_cos_eff;  // [QQ_mode][GenReco][LPFO][cut]
      unordered_map<TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > > > h2_dEdx_eff;  // [QQ_mode][GenReco][LPFO][TruthID][cut][hist]

    // h2 hist
      vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
      unordered_map<TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > h2_dEdx; // [QQ_mode][LPFO][TruthID][hist]

    private:
    // Lists
      unordered_map<TString, unordered_map< TString, TList* > > hList1; // [QQ_mode][category]
      unordered_map<TString, unordered_map< TString, TList* > > hList2; // [QQ_mode][category]

    // PFO Tools
      PFOTools _pt;

  };
}
#endif