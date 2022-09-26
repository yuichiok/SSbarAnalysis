#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "TreeStructures.hh"
#ifndef _TreeWriter_hh
#define _TreeWriter_hh

class TreeWriter 
{
  public:
    //
    //	Constants
    //
    
    //
    //	Constructors
    //
    TreeWriter ();
    virtual ~TreeWriter () {};
    //
    //	Methods
    //
    
    void InitializeLPFOTree(TTree * tree, Tree_SSbar & data);
  private:
    //
    //	Data
    //
    
    //
    //	Private methods
    //
};

#endif
