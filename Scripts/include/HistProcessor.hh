#ifndef GUARD_HistProcessor_h
#define GUARD_HistProcessor_h

/*------------------------------------------------------------------------------
   HistProcessor
 Created : 2022-09-02  okugawa
 Main class of HistProcessor program. Created to handle input and running of
 ntuple processing.
 Takes input variables (datasets, configurations, etc) and sets up the
 appropriate classes to handle each portion of the analysis process.
------------------------------------------------------------------------------*/
#include <TString.h>
#include "EventAnalyzer.hh"
#include "TreeIterator.hh"

class HistProcessor
{
  public:
    HistProcessor(TString o="", int me = -1);
    ~HistProcessor(){}

    EventAnalyzer  eAnalyzer;
    TreeIterator   tIter;
    TString        options;        // Extra options for this processing.
    int            maxEvents;      // Max number of events to run over in this ntuple

    TFile          *ntupleFile;

  private: 

};

#endif