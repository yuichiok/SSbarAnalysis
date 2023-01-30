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
#include "AnalysisConfig.hh"

using std::vector;

class PFOTools
{
  public:
    PFOTools();
    PFOTools( MC_QQbar *ptr, TString fnac );
    PFOTools( MC_QQbar *ptr_mc, PFO_QQbar *ptr, TString fnac );
    virtual ~PFOTools() {};
    virtual void InitializeMCTools( MC_QQbar *ptr );
    virtual void InitializePFOTools( MC_QQbar *ptr_mc, PFO_QQbar *ptr );

    enum   ChargeConfig { kSame, kOpposite };

    virtual vector<PFO_Info>  SortJet  ( vector<PFO_Info> jet );
    virtual Bool_t            ValidPFO ();

    virtual vector<PFO_Info>  Get_Valid_PFOs    ();
    virtual vector<PFO_Info>  GetJet            ( int ijet );
    virtual vector<PFO_Info>  GetSortedJet      ( int ijet );
    virtual PFO_Info          Get_KLPFO         ( int ijet );
    virtual Int_t             Get_dEdx_dist_PID ( Float_t kdEdx_dist, Float_t pidEdx_dist, Float_t pdEdx_dist );
    static  Bool_t            isKaon            ( PFO_Info iPFO );
    static  Bool_t            isPion            ( PFO_Info iPFO );
    static  Bool_t            isProton          ( PFO_Info iPFO );

  // LPFO checks
    virtual Bool_t           is_charge_config ( ChargeConfig cc, Int_t charge0, Int_t charge1 );

    virtual Bool_t           LPFO_Quality_checks    ( PFO_Info iPFO );
    virtual Bool_t           is_momentum           ( PFO_Info iPFO, Float_t MINP, Float_t MAXP );
    virtual Bool_t           is_tpc_hits           ( PFO_Info iPFO, Int_t MIN_TPC_HITS );
    virtual Bool_t           is_offset_small       ( PFO_Info iPFO, Int_t MAX_OFFSET );
    virtual Bool_t           is_dEdxdist_bad       ( Float_t e_dist, Float_t mu_dist, Float_t pi_dist, Float_t k_dist, Float_t p_dist );

    AnalysisConfig _anCfg;

  // MC gen info
    MC_Info    mc_quark[2];
    MC_Info    mc_stable[1000];

  // List of PFOs in jets
    vector<PFO_Info> PFO_jet[2];

  // Leading/Sub-Leading PFOs
    PFO_Info LPFO[2];
    PFO_Info KLPFO[2];
    std::vector<PFO_Info> SPFOs[2];
    std::vector<PFO_Info> SPFOs_K[2];
    std::vector<PFO_Info> SPFOs_Pi[2];


  private:
    MC_QQbar  *mc_data  ;
    PFO_QQbar *data     ;
    PFO_Info   PFO      ;

};

#endif