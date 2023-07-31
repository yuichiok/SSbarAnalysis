
/*------------------------------------------------------------------------------
TreeIterator.cpp
 Created : 2022-09-05  okugawa
------------------------------------------------------------------------------*/

#include <iomanip>
#include <iostream>
#include <string>
#include <TString.h>
#include "TreeIterator.hh"
#include "HistManager.hh"

using std::cout;  using std::endl;   using std::setw;   using std::string;

namespace QQbarAnalysis
{
  void TreeIterator::Begin(TTree * /*tree*/){}

  void TreeIterator::SlaveBegin(TTree * /*tree*/)
  {
    // Initialize log, counters
      nEntries          = 0;
      nEntriesProcessed = 0;

      eAnalyzer.CreateFile();
      eAnalyzer.InitHists();

  }


  void TreeIterator::Init(TTree *tree)
  {
    // Initialize the current tree.
      if (!tree) return;
      fChain = tree;              // Commented, don't think I need to access tree information outside of this function.
      fChain->SetMakeClass(1);    // ??? Can't seem to get a good answer about what this is.
                                  //    I think it recreates classes stored in the root tree as opposed to just mapping variables from the object.
      nEntries = tree->GetEntries();
      finalEntry = nEntries-1;

    // Initialize Event Handler, adding the criteria of each HistoMaker to it's list of criteria.
      eAnalyzer.InitReadTree(fChain);  

  }


  Bool_t TreeIterator::Notify()
  {
      Init(fChain);
      return kTRUE;
  }


  Bool_t TreeIterator::Process(Long64_t entry)
  {
    // Load current entry
      fChain->GetEntry(entry);
      nEntriesProcessed++;

      if ( entry % 1000 == 0 ) cout << "    [TreeIterator] Event " << entry << endl;
    // Evaluate the criteria for this entry
      if ( eAnalyzer.Select( EventAnalyzer::kQQ ) ){

        eAnalyzer.Jet_sum_n_acol();

        if( eAnalyzer.Select( EventAnalyzer::kMC ) ){
          eAnalyzer.AnalyzeGen(entry);
        }
        if( eAnalyzer.Select( EventAnalyzer::kReco ) ){
          eAnalyzer.AnalyzeReco(entry);
        }
        
      }

      return true;
  }

  void TreeIterator::SlaveTerminate(){}
  void TreeIterator::Terminate()
  {
    // Write Trees
      eAnalyzer.WriteFile();
  }
}