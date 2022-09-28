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
#include "PFOTools.hh"
#include "TreeStructures.hh"
#include "TreeReader.hh"
#include "TreeWriter.hh"
#include "FileSelector.hh"

#include "TreeVariables.hh"
ClassImp(TreeVariables)

class EventAnalyzer
{
  public:
    EventAnalyzer(TString o="");
    virtual ~EventAnalyzer(){};

    enum       Selector { kMC, kLPFO };
    enum MCParticlePair { FIRST_ENTRY, kDD, kUU, kSS, kCC, kBB, kTT };

  // methods
    Bool_t           InitReadTree( TTree* ); // Maps class variables to an input TTree.
    void             InitWriteTree(); // 
    void             WriteFile(); // 
    void             Analyze( Long64_t entry );
    Bool_t           Select( Selector s );          // Evaluates the class' list of event selection criteria
    Bool_t           GenPairPicker( Float_t mc_particle, MCParticlePair pair );
    Bool_t           ISRPicker( Float_t Kvcut );
    virtual Bool_t   Notify();



  // Running Variables
    TString  options;  // Options input with TreeIterator.
    TTree   *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t    fCurrent; //!current Tree number in a TChain

  // Extra data
    long patEventsAnalyzed;     // Number of events that were processed to make the Ntuple.
    long entriesInNtuple  ;     // Number of events that were processed to make the Ntuple.

  // Fixed size dimensions of array or collections stored in the TTree if any.
    FileSelector  _fs;


  private:

    MC_QQbar      _mc      ;
    Jet_QQbar     _jet     ;
    PFO_QQbar     _pfo     ;
    Branch_QQbar  _branch  ;

    TreeWriter writer;
    TFile * _hfile;
    TTree * _hTree;
    TTree * _hTree_LPFO;
    TTree * _hTree_LPFO_KK;
    TTree * _hTree_LPFO_KPi;

    TString       _hfilename;

    TreeVariables    _tree_lpfo;
    TreeVariables    _tree_lpfo_kk;
    TreeVariables    _tree_lpfo_kpi;

};

#endif