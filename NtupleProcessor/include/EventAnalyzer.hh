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
#include "MapTString.hh"

#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <vector>
#include <fstream>

using std::unordered_map;

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
      void             InitWeights();
      
      void             CreateFile();
      void             WriteFile();

      void             AnalyzeGen( Long64_t entry );
      void             AnalyzeReco( Long64_t entry );
      
      Bool_t           Select( Selector s );          // Evaluates the class' list of event selection criteria
      vector<Bool_t>   SelectRecoMC( TString recomc, VectorTools v_reco[2], VectorTools v_gen[2] );
      Bool_t           isSignal();                      // Checks if the event is qqbar
      Bool_t           GenPairPicker( Float_t mc_particle, std::vector<int> input_gen );

      unordered_map<TString, Bool_t> TriggerMap( PFOTools pfot, TString lmode, unordered_map< int, vector<PFO_Info> > subjet_pair, TString gen_reco );
      Bool_t           Cut_ESum( VectorTools v[2] );
      Bool_t           Cut_ACol( VectorTools v[2] );
      Bool_t           Cut_ISR ( VectorTools v[2] );

      Bool_t           Cut_PhotonJets( TString recomc, TString debug );
      Float_t          Cut_SinACol( VectorTools v[2] );
      Float_t          Cut_invM( VectorTools v[2] );
      Float_t          Cut_d23( TString recomc );

      /////
      bool PreSelection(int type=0, float acolcut=0.3);
      
      // calculation of photon related quantities
      double  a_npfo[2]={-1};
      double  a_npfo_photon[2]={-1};
      double  a_npfo_charge[2]={-1};
      double  a_photonjet_E[2]={0};
      double  a_photonjet_costheta[2]={-2};
      void PFOphotonQuantities();

      std::vector< float > getAngles(std::vector< float > & direction);
      float getModule(vector< float > & v);
      std::vector< float > getDirection(vector<float> & vectorPoint);
      float GetCostheta(std::vector<float> & vectorPoint);
      float GetSinacol(TVector3 v1, TVector3 v2);
      float AcolValue();

      /////



      void             ClearStructs();
      virtual Bool_t   Notify();

    // Gadgets
      void             ResolutionAnalysis( PFOTools pfot, PFOTools mct );
      unordered_map<TString, Int_t> Gen_Reco_Stats_Stable( PFOTools pfot, PFOTools mct, TString lmode , Float_t cos_min, Float_t cos_max );

    // Histogram extractor
      void             ProcessDoubleTag( PFOTools pfot, PFOTools mct, unordered_map< TString, unordered_map<TString, Bool_t> > cuts, unordered_map< int, vector<PFO_Info> > subjet_pair );
      void             ProcessDoubleTagEfficiency( PFOTools pfot, PFOTools mct, unordered_map< TString, unordered_map<TString, Bool_t> > cuts, unordered_map< int, vector<PFO_Info> > subjet_pair, TString gen_reco );

      void             Jet_sum_n_acol();

    // Running Variables
      TString  options;  // Options input with TreeIterator.
      TTree   *fChain;   //!pointer to the analyzed TTree or TChain
      Int_t    fCurrent; //!current Tree number in a TChain

    // Extra data
      long patEventsAnalyzed;     // Number of events that were processed to make the Ntuple.
      long entriesInNtuple  ;     // Number of events that were processed to make the Ntuple.

    private:

      Int_t ientry = -1;

      Int_t _check_pt = 0;

      Int_t _qmodePDG = 0;
      TString  _qmode = "";
      unordered_map< Int_t, TString > _qmode_map = {
        {0,"bg"}, {1,"dd"}, {2,"uu"}, {3,"ss"}, {4,"cc"}, {5,"bb"}
      };

      MC_QQbar      _mc      ;
      Jet_QQbar     _jet     ;
      VTX_QQbar     _vtx     ;
      PFO_QQbar     _pfo     ;
      Branch_QQbar  _branch  ;

      PFOTools      _pt;

      TFile * _hfile;

      TTree * _hTree;
      Tree_Data _data;

      TEvent           _eve;
      TreeVariables    _stats_lpfo;
      LPFO_Info        _data_lpfo;

      TString _hfilename;

    // Fixed size dimensions of array or collections stored in the TTree if any.
      FileSelector  _fs;
      HistManager   _hm;

      AnalysisConfig _anCfg;
      TString       _config;

      unordered_map< TString, TGraphErrors* > _gdedx;
      TGraphErrors * _gdedxPi;
      TGraphErrors * _gdedxK;

  };
}
#endif