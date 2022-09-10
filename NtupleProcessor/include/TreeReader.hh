#ifndef GUARD_TreeReader_h
#define GUARD_TreeReader_h

/*------------------------------------------------------------------------------
   TreeReader
 Created : 2022-09-02  okugawa
 Main class of TreeReader program.
------------------------------------------------------------------------------*/

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "TreeStructures.hh"

class TreeReader
{
  public:
    TreeReader();
    virtual ~TreeReader() {};

    void InitializeMCReadTree (TTree *tree, MC_QQbar & _data,  Branch_QQbar & _branch);
    void InitializeJetReadTree(TTree *tree, Jet_QQbar & _data, Branch_QQbar & _branch);
    void InitializePFOReadTree(TTree *tree, PFO_QQbar & _data, Branch_QQbar & _branch);

  private: 

};

#endif