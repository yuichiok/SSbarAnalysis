#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "PFOTools.hh"
#include "TreeStructures.hh"
#include "TreeVariables.hh"
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
    
    void InitializeLPFOTree( TTree * tree, Tree_SSbar & data);
    void WriteLPFOVariables( PFOTools pt, PFO_QQbar *pqq, TreeVariables *data );

  private:
    //
    //	Data
    //
    
    //
    //	Private methods
    //
};

#endif
