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
#include "HistManager.hh"

class EventAnalyzer
{
  public:
    EventAnalyzer(TString o="");
    virtual ~EventAnalyzer(){};

    enum       Selector { kQQ, kMC, kReco };
    enum MCParticlePair { FIRST_ENTRY, kDD, kUU, kSS, kCC, kBB, kTT };

  // methods
    Bool_t           InitReadTree( TTree* ); // Maps class variables to an input TTree.
    void             InitWriteTree();
    void             InitHists();
    
    void             CreateFile();
    void             WriteFile();

    void             AnalyzeGen();
    void             Analyze( Long64_t entry );
    
    Bool_t           Select( Selector s );          // Evaluates the class' list of event selection criteria
    Bool_t           GenPairPicker( Float_t mc_particle, MCParticlePair pair );

    Bool_t           Cut_ESum( VectorTools v[2] );
    Bool_t           Cut_ACol( VectorTools v[2] );
    Bool_t           Cut_ISR ( VectorTools v[2] );
    void             ClearStructs();
    virtual Bool_t   Notify();

  // Gadgets
    void             Gen_Reco_Stats( PFOTools mct, PFOTools pfot );

  // Histogram extractor
    void             PolarAngleGen(PFOTools mct);
    void             PolarAngle(PFOTools pfot, Bool_t b_reco);

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


  private:

    MC_QQbar      _mc      ;
    Jet_QQbar     _jet     ;
    PFO_QQbar     _pfo     ;
    Branch_QQbar  _branch  ;

    TreeWriter writer;
    TFile * _hfile;

    TTree * _hTree;
    Tree_Data _data;

    TEvent           _eve;
    TreeVariables    _stats_lpfo;
    LPFO_Info        _data_lpfo;

    TString _hfilename;

};

#endif