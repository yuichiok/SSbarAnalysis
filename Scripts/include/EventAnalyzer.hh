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
#include "FileSelector.hh"

class EventAnalyzer
{
  public:
    EventAnalyzer(TString o="");
    virtual ~EventAnalyzer(){};

  // methods
    Bool_t           InitReadTree( TTree* ); // Maps class variables to an input TTree.
    void             Analyze( Long64_t entry );
    virtual Bool_t   Notify();



  // Running Variables
    TString  options;  // Options input with TreeIterator.
    TTree   *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t    fCurrent; //!current Tree number in a TChain

  // Extra data
    long patEventsAnalyzed;     // Number of events that were processed to make the Ntuple.
    long entriesInNtuple  ;     // Number of events that were processed to make the Ntuple.

  // Fixed size dimensions of array or collections stored in the TTree if any.

  private:
    FileSelector _fs;

    MC_QQbar*      _mc      ;
    Jet_QQbar     _jet     ;
    PFO_QQbar     _lpfo     ;
    Branch_QQbar  _branch  ;

    // TreeWriter writer;
    // TFile * _hfile;
    // TTree * _hTree;
    // TTree * _hTree_LPFO;
    // TTree * _hTree_LPFO_KK;
    // TTree * _hTree_LPFO_KPi;

    // TString       _hfilename;

    TEvent*           _eve;
    TreeVariables*    _stats_lpfo;
    LPFO_Info*        _data_lpfo = nullptr;

    TBranch          *_b_data_lpfo;

};

#endif