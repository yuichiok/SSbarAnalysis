#ifndef GUARD_PFOTools_h
#define GUARD_PFOTools_h

/*------------------------------------------------------------------------------
   PFOTools
 Created : 2022-09-11  okugawa
 Main class of PFOTool program.
------------------------------------------------------------------------------*/

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "../../SSbarLibrary/TreeStructures.hh"

class PFOTools
{
  public:
    PFOTools();
    virtual ~PFOTools() {};

  private:
    PFO_QQbar _pfo;

};

#endif