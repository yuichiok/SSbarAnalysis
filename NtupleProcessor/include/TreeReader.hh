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

class TreeReader
{
  public:
    TreeReader();
    virtual ~TreeReader() {};

    void InitializeReadTree(TTree *tree, StatsData_QQbar & data);

  private: 

};

#endif