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
#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <fstream>
#include "TreeStructures.hh"
#include "TreeReader.hh"

class EventAnalyzer
{
  public:
    EventAnalyzer(TString o="");
    virtual ~EventAnalyzer(){};

    enum       Selector { kMC, kLPFO };
    enum   ChargeConfig { kSame, kOpposite };
    enum MCParticlePair { FIRST_ENTRY, kDD, kUU, kSS, kCC, kBB, kTT };

  // methods
    Bool_t           MapTree( TTree* ); // Maps class variables to an input TTree.
    void             Analyze( Long64_t entry );
    Bool_t           Select( Selector s );          // Evaluates the class' list of event selection criteria
    Bool_t           GenPairPicker( Float_t mc_particle, MCParticlePair pair );
    Bool_t           ISRPicker( Float_t Kvcut );
    virtual Bool_t   Notify();

  // LPFO checks
    Bool_t           is_charge_config ( ChargeConfig cc );

    Bool_t           PFO_Quality_checks    ();
    Bool_t           is_momentum           ( PFO_Info iPFO, Float_t MINP, Float_t MAXP );
    Bool_t           is_tpc_hits           ( PFO_Info iPFO, Int_t MIN_TPC_HITS );
    Bool_t           is_offset_small       ( PFO_Info iPFO, Int_t MAX_OFFSET );


  // Running Variables
    TString  options;  // Options input with TreeIterator.
    TTree   *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t    fCurrent; //!current Tree number in a TChain

  // Extra data
    long patEventsAnalyzed;     // Number of events that were processed to make the Ntuple.
    long entriesInNtuple  ;     // Number of events that were processed to make the Ntuple.

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Leading PFO
    PFO_Info LPFO[2];
    std::vector<PFO_Info> SPFOs[2];

  private:

    MC_QQbar      _mc      ;
    Jet_QQbar     _jet     ;
    PFO_QQbar     _pfo     ;
    Branch_QQbar  _branch  ;

};

#endif