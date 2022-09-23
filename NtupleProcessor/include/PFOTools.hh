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

    virtual vector<PFO_Info>  SortJet  ( vector<PFO_Info> jet );
    virtual Bool_t            ValidPFO ();

    virtual vector<PFO_Info>  GetJet            ( int ijet );
    virtual vector<PFO_Info>  GetSortedJet      ( int ijet );
    virtual Int_t             Get_dEdx_dist_PID ( Float_t kdEdx_dist, Float_t pidEdx_dist, Float_t pdEdx_dist );

    // List of PFOs in jets
    vector<PFO_Info> PFO_jet[2];

  private:
    PFO_QQbar *data     ;
    PFO_Info   PFO      ;

};

#endif