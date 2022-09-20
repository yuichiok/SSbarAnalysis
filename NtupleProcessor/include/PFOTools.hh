#ifndef GUARD_PFOTools_h
#define GUARD_PFOTools_h

/*------------------------------------------------------------------------------
   PFOTools
 Created : 2022-09-11  okugawa
 Main class of PFOTool program.
------------------------------------------------------------------------------*/

#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "TreeStructures.hh"

using std::vector;

class PFOTools
{
  public:
    PFOTools( PFO_QQbar *ptr );
    virtual ~PFOTools() {};

    virtual void              SortJet  ( vector<PFO_Info> jet );
    virtual bool              ValidPFO ();

    virtual vector<PFO_Info>  GetJet        ( int ijet );
    virtual vector<PFO_Info>  GetSortedJet  ( int ijet );

    // List of PFOs in jets
    vector<PFO_Info> PFO_jet[2];

  private:
    PFO_QQbar *data     ;
    PFO_Info   PFO      ;

};

#endif