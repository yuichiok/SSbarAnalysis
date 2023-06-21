#ifndef GUARD_EventAnalyzer_h
#define GUARD_EventAnalyzer_h

/*------------------------------------------------------------------------------
   NtupleProcessor
 Created : 2022-09-02  okugawa
 Main class of NtupleProcessor program. Created to handle input and running of
 ntuple processing.
 Takes input variables (datasets, configurations, etc) and sets up the
 appropriate classes to handle each portion of the analysis process.
------------------------------------------------------------------------------*/
#include "PFOTools.hh"
#include "TreeStructures.hh"
#include "TreeReader.hh"
#include "TreeWriter.hh"
#include "FileSelector.hh"
#include "HistManager.hh"
#include "AnalysisConfig.hh"

#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <map>
#include <vector>
#include <fstream>

namespace QQbarAnalysis
{
  class EventAnalyzer
  {
    public:
      EventAnalyzer(TString input, TString fnac, TString o="");
      virtual ~EventAnalyzer(){};

      enum Selector       { kQQ, kMC, kReco };
      enum MCParticlePair { FIRST_ENTRY, kDD, kUU, kSS, kCC, kBB, kTT };
      enum SelectID       { kKaon, kPion, kProton, kElectron, kMuon };
      enum PDGConfig      { noKPi, K_K, K_Pi, Pi_Pi };

    // methods
      Bool_t           InitReadTree( TTree* ); // Maps class variables to an input TTree.
      void             InitHists();
      
      void             CreateFile();
      void             WriteFile();

      void             AnalyzeGen();
      void             AnalyzeGenReco(PFOTools mct, PFOTools pfot);
      void             AnalyzeReco( Long64_t entry );
      
      Bool_t           Select( Selector s );          // Evaluates the class' list of event selection criteria
      Bool_t           GenPairPicker( Float_t mc_particle, std::vector<int> input_gen );

      Bool_t           Cut_ESum( VectorTools v[2] );
      Bool_t           Cut_ACol( VectorTools v[2] );
      Bool_t           Cut_ISR ( VectorTools v[2] );
      void             ClearStructs();
      virtual Bool_t   Notify();

    // Gadgets
      Int_t           *Gen_Reco_Stats_Stable( PFOTools mct, PFOTools pfot, SelectID pid , Float_t cos_min, Float_t cos_max );
      Int_t           *Gen_Reco_Stats_Cheat( PFOTools mct, PFOTools pfot, SelectID pid , Float_t cos_min, Float_t cos_max );
      Float_t         *Get_Stable_Purity( Int_t *N_Ks );

      void             Count_Particle( PFO_Info ipfo, Int_t pdg, TH1F *h_n_reco, TH1F *h_n_gen );
      Float_t         *Particle_Ratios( TH1F *h_n_particles[], Int_t mode );

    // Histogram extractor
      void             PolarAngleGen( PFOTools mct );
      void             Mom_Polar_Gen( PFOTools mct, PFOTools pfot );
      void             ProcessDoubleTag( PFOTools pfot, PFOTools mct, vector<Bool_t> cuts[3], PDGConfig double_tag[3]);
      void             ProcessDoubleTag( PFOTools pfot, PFOTools mct, std::map< TString, std::map<TString, Bool_t> > cuts );
      void             PolarAngle_acc_rej( PFOTools pfot, vector<Bool_t> cuts, Bool_t ss_config );

      void             Jet_sum_n_acol();

    // Running Variables
      TString  options;  // Options input with TreeIterator.
      TTree   *fChain;   //!pointer to the analyzed TTree or TChain
      Int_t    fCurrent; //!current Tree number in a TChain

    // Extra data
      long patEventsAnalyzed;     // Number of events that were processed to make the Ntuple.
      long entriesInNtuple  ;     // Number of events that were processed to make the Ntuple.

    // Fixed size dimensions of array or collections stored in the TTree if any.
      FileSelector  _fs;
      HistManager   _hm;

      AnalysisConfig _anCfg;
      TString       _config;

    // List of cut names
      vector<TString> PFO_mode  = {"K","Pi"};
      vector<TString> cut_names = {"jet_association","quality","SPFO","charge","PID"};

    private:

      Int_t ientry = -1;

      MC_QQbar      _mc      ;
      Jet_QQbar     _jet     ;
      VTX_QQbar     _vtx     ;
      PFO_QQbar     _pfo     ;
      Branch_QQbar  _branch  ;

      TFile * _hfile;

      TTree * _hTree;
      Tree_Data _data;

      TEvent           _eve;
      TreeVariables    _stats_lpfo;
      LPFO_Info        _data_lpfo;

      TString _hfilename;

  };
}
#endif