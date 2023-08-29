#ifndef GUARD_PFOTools_h
#define GUARD_PFOTools_h

/*------------------------------------------------------------------------------
   PFOTools
 Created : 2022-09-11  okugawa
 Main class of PFOTool program.
------------------------------------------------------------------------------*/
#include "TreeStructures.hh"
#include "AnalysisConfig.hh"
#include "MapTString.hh"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>

#include <vector>

using std::vector; using std::unordered_map;

namespace QQbarAnalysis
{
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
      enum   ParticleID   { kKaon, kPion, kProton };

      virtual vector<PFO_Info>  SortJet  ( vector<PFO_Info> jet );

      virtual vector<PFO_Info>  GetJet            ( int ijet );
      virtual vector<PFO_Info>  GetSortedJet      ( int ijet );
      virtual vector<PFO_Info>  GetSubjet         ( int ijet, TString lmode );
      virtual Int_t             Get_dEdx_dist_PID ( Float_t kdEdx_dist, Float_t pidEdx_dist, Float_t pdEdx_dist );
      virtual Float_t           Get_dEdx_dist     ( PFO_Info iPFO, TString particle );
      virtual Bool_t            isKaon            ( PFO_Info iPFO );
      virtual Bool_t            isPion            ( PFO_Info iPFO );
      virtual Bool_t            isProton          ( PFO_Info iPFO );
      virtual Bool_t            is_cheatNoOthers  ( PFO_Info iPFO );
      virtual Bool_t            is_PID            ( TString lmode, PFO_Info iPFO );

    // LPFO checks
      virtual Bool_t           is_jet_mult_non0 ();
      virtual Bool_t           is_charge_config ( ChargeConfig cc, Int_t charge0, Int_t charge1 );

      virtual Bool_t           is_momentum           ( PFO_Info iPFO, Float_t MINP, Float_t MAXP );
      virtual Bool_t           is_tpc_hits           ( PFO_Info iPFO, Int_t MIN_TPC_HITS );
      virtual Bool_t           is_offset_small       ( PFO_Info iPFO, Int_t MAX_OFFSET );
      virtual Bool_t           is_dEdxdist_bad       ( Float_t e_dist, Float_t mu_dist, Float_t pi_dist, Float_t k_dist, Float_t p_dist );

    // MC gen info
      MC_Info    mc_quark[2];
      MC_Info    mc_stable[1000];

    // List of PFOs in jets
      vector<PFO_Info> Valid_PFOs;

      vector<PFO_Info> PFO_jet[2];
      unordered_map<int, vector<PFO_Info> > PFO_sorted_jet;
      unordered_map< TString, unordered_map<int, vector<PFO_Info> > > PFO_subjet;
      unordered_map< TString, unordered_map<int, vector<PFO_Info> > > PFO_unsorted_subjet_cheat;
      unordered_map< TString, unordered_map<int, vector<PFO_Info> > > PFO_subjet_cheat;

      vector<PFO_Info> LPFO;      

    // PFO modes and types
      vector<TString> PFO_mode  = {"K","Pi"};
      vector<TString> PFO_type  = {"K","Pi", "p", "e", "mu"};
      unordered_map<Int_t,TString> PFO_type_map = { {321,"K"}, {211,"Pi"}, {2212,"p"}, {11,"e"}, {13,"mu"} };
      unordered_map<TString, Int_t> PFO_type_map_rev = { {"K",321}, {"Pi",211}, {"p",2212}, {"e",11}, {"mu",13} };


    private:
      MC_QQbar  *mc_data  ;
      PFO_QQbar *data     ;
      PFO_Info   PFO      ;

      TF1 *fCorrection = new TF1("fCorrection","-2*(1.6 - pol8)",-1,1);
      Double_t pars[9] = {0.599455,-0.000377453,-0.120892,0.00118809,0.82,-0.000671831,-1.72613,-0.000781913,1.04};


      AnalysisConfig _anCfg;

  };
}
#endif