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
    PFOTools();
    PFOTools( PFO_QQbar *ptr );
    virtual ~PFOTools() {};

    enum   ChargeConfig { kSame, kOpposite };

    virtual vector<PFO_Info>  SortJet  ( vector<PFO_Info> jet );
    virtual Bool_t            ValidPFO ();

    virtual vector<PFO_Info>  GetJet            ( int ijet );
    virtual vector<PFO_Info>  GetSortedJet      ( int ijet );
    virtual Int_t             Get_dEdx_dist_PID ( Float_t kdEdx_dist, Float_t pidEdx_dist, Float_t pdEdx_dist );
    virtual Bool_t            isKaon            ( PFO_Info iPFO );
    virtual Bool_t            isPion            ( PFO_Info iPFO );
    virtual Bool_t            isProton          ( PFO_Info iPFO );

  // LPFO checks
    Bool_t           is_charge_config ( ChargeConfig cc );

    Bool_t           PFO_Quality_checks    ( PFO_Info iPFO );
    Bool_t           is_momentum           ( PFO_Info iPFO, Float_t MINP, Float_t MAXP );
    Bool_t           is_tpc_hits           ( PFO_Info iPFO, Int_t MIN_TPC_HITS );
    Bool_t           is_offset_small       ( PFO_Info iPFO, Int_t MAX_OFFSET );


  // List of PFOs in jets
    vector<PFO_Info> PFO_jet[2];

  // Leading/Sub-Leading PFOs
    PFO_Info LPFO[2];
    std::vector<PFO_Info> SPFOs[2];
    std::vector<PFO_Info> SPFOs_K[2];


  private:
    PFO_QQbar *data     ;
    PFO_Info   PFO      ;

};

#endif