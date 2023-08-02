#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "PFOTools.hh"
#include "TreeStructures.hh"
#ifndef _TreeWriter_hh
#define _TreeWriter_hh

namespace QQbarAnalysis
{
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
      
      void InitializeDataTree( TTree * tree, Tree_Data & data);

    private:
      //
      //	Data
      //
      
      //
      //	Private methods
      //
  };
}
#endif
