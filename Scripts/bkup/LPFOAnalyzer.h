//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 29 16:42:47 2022 by ROOT version 6.26/02
// from TTree event/tree
// found on file: ../rootfiles/tmp_root/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n002.d_dstm_15162_000.ss.tmp.root
//////////////////////////////////////////////////////////

#ifndef LPFOAnalyzer_h
#define LPFOAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class LPFOAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //TEvent          *Event;
   Int_t           eve_valid_lpfo;
 //MC_QQbar        *MC;
   Float_t         mc_quark_E[2];
   Float_t         mc_quark_px[2];
   Float_t         mc_quark_py[2];
   Float_t         mc_quark_pz[2];
   Float_t         mc_quark_m[2];
   Float_t         mc_quark_pdg[2];
   Float_t         mc_quark_charge[2];
   Float_t         mc_ISR_E[2];
   Float_t         mc_ISR_px[2];
   Float_t         mc_ISR_py[2];
   Float_t         mc_ISR_pz[2];
   Float_t         mc_ISR_m[2];
   Float_t         mc_ISR_pdg[2];
   Float_t         mc_ISR_charge[2];
   Int_t           mc_quark_ps_n;
   Float_t         mc_quark_ps_E[1000];
   Float_t         mc_quark_ps_px[1000];
   Float_t         mc_quark_ps_py[1000];
   Float_t         mc_quark_ps_pz[1000];
   Float_t         mc_quark_ps_m[1000];
   Float_t         mc_quark_ps_pdg[1000];
   Float_t         mc_quark_ps_charge[1000];
   Float_t         mc_quark_ps_y12;
   Float_t         mc_quark_ps_y23;
   Float_t         mc_quark_ps_d12;
   Float_t         mc_quark_ps_d23;
   Float_t         mc_quark_ps_jet_E[2];
   Float_t         mc_quark_ps_jet_px[2];
   Float_t         mc_quark_ps_jet_py[2];
   Float_t         mc_quark_ps_jet_pz[2];
   Int_t           mc_stable_n;
   Float_t         mc_stable_E[1000];
   Float_t         mc_stable_px[1000];
   Float_t         mc_stable_py[1000];
   Float_t         mc_stable_pz[1000];
   Float_t         mc_stable_m[1000];
   Int_t           mc_stable_pdg[1000];
   Float_t         mc_stable_charge[1000];
   Int_t           mc_stable_isoverlay[1000];
   Int_t           mc_stable_isisr[1000];
   Float_t         mc_stable_y12;
   Float_t         mc_stable_y23;
   Float_t         mc_stable_d12;
   Float_t         mc_stable_d23;
   Float_t         mc_stable_jet_E[2];
   Float_t         mc_stable_jet_px[2];
   Float_t         mc_stable_jet_py[2];
   Float_t         mc_stable_jet_pz[2];
 //TreeVariables   *Stats_LPFO;
   Int_t           lpfo_match[2];
   Int_t           lpfo_truejet_pdg[2];
   Int_t           lpfo_truejet_type[2];
   Int_t           lpfo_pdgcheat[2];
   Int_t           lpfo_nparents[2];
   Int_t           lpfo_pdgcheat_parent[2][1000];
   Float_t         lpfo_E[2];
   Float_t         lpfo_px[2];
   Float_t         lpfo_py[2];
   Float_t         lpfo_pz[2];
   Float_t         lpfo_m[2];
   Int_t           lpfo_type[2];
   Int_t           lpfo_isoverlay[2];
   Int_t           lpfo_isisr[2];
   Int_t           lpfo_vtx[2];
   Int_t           lpfo_charge[2];
   Int_t           lpfo_ntracks[2];
   Int_t           lpfo_tpc_hits[2];
   Float_t         lpfo_dedx[2];
   Float_t         lpfo_dedxerror[2];
   Float_t         lpfo_d0[2];
   Float_t         lpfo_d0error[2];
   Float_t         lpfo_z0[2];
   Float_t         lpfo_z0error[2];
   Float_t         lpfo_phi[2];
   Float_t         lpfo_phierror[2];
   Float_t         lpfo_omega[2];
   Float_t         lpfo_omegaerror[2];
   Float_t         lpfo_tanlambda[2];
   Float_t         lpfo_tanlambdaerror[2];
   Float_t         lpfo_chi2[2];
   Float_t         lpfo_ndf[2];
   Float_t         lpfo_vtxpt[2][3];
   Float_t         lpfo_endpt[2][3];
   Int_t           lpfo_pid[2];
   Float_t         lpfo_pid_likelihood[2];
   Float_t         lpfo_pid_eprob[2];
   Float_t         lpfo_pid_muprob[2];
   Float_t         lpfo_pid_piprob[2];
   Float_t         lpfo_pid_kprob[2];
   Float_t         lpfo_pid_pprob[2];
   Float_t         lpfo_pid_hprob[2];
   Int_t           lpfo_piddedx[2];
   Float_t         lpfo_piddedx_likelihood[2];
   Float_t         lpfo_piddedx_eprob[2];
   Float_t         lpfo_piddedx_muprob[2];
   Float_t         lpfo_piddedx_piprob[2];
   Float_t         lpfo_piddedx_kprob[2];
   Float_t         lpfo_piddedx_pprob[2];
   Float_t         lpfo_piddedx_hprob[2];
   Float_t         lpfo_piddedx_e_dedxdist[2];
   Float_t         lpfo_piddedx_mu_dedxdist[2];
   Float_t         lpfo_piddedx_pi_dedxdist[2];
   Float_t         lpfo_piddedx_k_dedxdist[2];
   Float_t         lpfo_piddedx_p_dedxdist[2];
   Float_t         lpfo_piddedx_e_lkhood[2];
   Float_t         lpfo_piddedx_mu_lkhood[2];
   Float_t         lpfo_piddedx_pi_lkhood[2];
   Float_t         lpfo_piddedx_k_lkhood[2];
   Float_t         lpfo_piddedx_p_lkhood[2];
   Float_t         lpfo_pidtof_p_at_calo[2];
   Float_t         lpfo_pidtof_closest_beta_0ps[2];
   Float_t         lpfo_pidtof_closest_beta_10ps[2];
   Float_t         lpfo_pidtof_closest_beta_50ps[2];
   Float_t         lpfo_pidtof_closest_beta_100ps[2];
   Float_t         lpfo_pidtof_fastest_beta_0ps[2];
   Float_t         lpfo_pidtof_fastest_beta_10ps[2];
   Float_t         lpfo_pidtof_fastest_beta_50ps[2];
   Float_t         lpfo_pidtof_fastest_beta_100ps[2];
   Float_t         lpfo_pidtof_cylfit_beta_0ps[2];
   Float_t         lpfo_pidtof_cylfit_beta_10ps[2];
   Float_t         lpfo_pidtof_cylfit_beta_50ps[2];
   Float_t         lpfo_pidtof_cylfit_beta_100ps[2];
   Float_t         lpfo_pidtof_closestfit_beta_0ps[2];
   Float_t         lpfo_pidtof_closestfit_beta_10ps[2];
   Float_t         lpfo_pidtof_closestfit_beta_50ps[2];
   Float_t         lpfo_pidtof_closestfit_beta_100ps[2];
 //LPFO_Info       *Data_LPFO;
   Int_t           lpfo_config;
   Float_t         lpfo_p_mag[2];
   Int_t           lpfo_dEdx_dist_pdg[2];
   Float_t         lpfo_cos[2];
   Float_t         lpfo_qcos[2];

   // List of branches
   TBranch        *b_Event_eve_valid_lpfo;   //!
   TBranch        *b_MC_mc_quark_E;   //!
   TBranch        *b_MC_mc_quark_px;   //!
   TBranch        *b_MC_mc_quark_py;   //!
   TBranch        *b_MC_mc_quark_pz;   //!
   TBranch        *b_MC_mc_quark_m;   //!
   TBranch        *b_MC_mc_quark_pdg;   //!
   TBranch        *b_MC_mc_quark_charge;   //!
   TBranch        *b_MC_mc_ISR_E;   //!
   TBranch        *b_MC_mc_ISR_px;   //!
   TBranch        *b_MC_mc_ISR_py;   //!
   TBranch        *b_MC_mc_ISR_pz;   //!
   TBranch        *b_MC_mc_ISR_m;   //!
   TBranch        *b_MC_mc_ISR_pdg;   //!
   TBranch        *b_MC_mc_ISR_charge;   //!
   TBranch        *b_MC_mc_quark_ps_n;   //!
   TBranch        *b_MC_mc_quark_ps_E;   //!
   TBranch        *b_MC_mc_quark_ps_px;   //!
   TBranch        *b_MC_mc_quark_ps_py;   //!
   TBranch        *b_MC_mc_quark_ps_pz;   //!
   TBranch        *b_MC_mc_quark_ps_m;   //!
   TBranch        *b_MC_mc_quark_ps_pdg;   //!
   TBranch        *b_MC_mc_quark_ps_charge;   //!
   TBranch        *b_MC_mc_quark_ps_y12;   //!
   TBranch        *b_MC_mc_quark_ps_y23;   //!
   TBranch        *b_MC_mc_quark_ps_d12;   //!
   TBranch        *b_MC_mc_quark_ps_d23;   //!
   TBranch        *b_MC_mc_quark_ps_jet_E;   //!
   TBranch        *b_MC_mc_quark_ps_jet_px;   //!
   TBranch        *b_MC_mc_quark_ps_jet_py;   //!
   TBranch        *b_MC_mc_quark_ps_jet_pz;   //!
   TBranch        *b_MC_mc_stable_n;   //!
   TBranch        *b_MC_mc_stable_E;   //!
   TBranch        *b_MC_mc_stable_px;   //!
   TBranch        *b_MC_mc_stable_py;   //!
   TBranch        *b_MC_mc_stable_pz;   //!
   TBranch        *b_MC_mc_stable_m;   //!
   TBranch        *b_MC_mc_stable_pdg;   //!
   TBranch        *b_MC_mc_stable_charge;   //!
   TBranch        *b_MC_mc_stable_isoverlay;   //!
   TBranch        *b_MC_mc_stable_isisr;   //!
   TBranch        *b_MC_mc_stable_y12;   //!
   TBranch        *b_MC_mc_stable_y23;   //!
   TBranch        *b_MC_mc_stable_d12;   //!
   TBranch        *b_MC_mc_stable_d23;   //!
   TBranch        *b_MC_mc_stable_jet_E;   //!
   TBranch        *b_MC_mc_stable_jet_px;   //!
   TBranch        *b_MC_mc_stable_jet_py;   //!
   TBranch        *b_MC_mc_stable_jet_pz;   //!
   TBranch        *b_Stats_LPFO_lpfo_match;   //!
   TBranch        *b_Stats_LPFO_lpfo_truejet_pdg;   //!
   TBranch        *b_Stats_LPFO_lpfo_truejet_type;   //!
   TBranch        *b_Stats_LPFO_lpfo_pdgcheat;   //!
   TBranch        *b_Stats_LPFO_lpfo_nparents;   //!
   TBranch        *b_Stats_LPFO_lpfo_pdgcheat_parent;   //!
   TBranch        *b_Stats_LPFO_lpfo_E;   //!
   TBranch        *b_Stats_LPFO_lpfo_px;   //!
   TBranch        *b_Stats_LPFO_lpfo_py;   //!
   TBranch        *b_Stats_LPFO_lpfo_pz;   //!
   TBranch        *b_Stats_LPFO_lpfo_m;   //!
   TBranch        *b_Stats_LPFO_lpfo_type;   //!
   TBranch        *b_Stats_LPFO_lpfo_isoverlay;   //!
   TBranch        *b_Stats_LPFO_lpfo_isisr;   //!
   TBranch        *b_Stats_LPFO_lpfo_vtx;   //!
   TBranch        *b_Stats_LPFO_lpfo_charge;   //!
   TBranch        *b_Stats_LPFO_lpfo_ntracks;   //!
   TBranch        *b_Stats_LPFO_lpfo_tpc_hits;   //!
   TBranch        *b_Stats_LPFO_lpfo_dedx;   //!
   TBranch        *b_Stats_LPFO_lpfo_dedxerror;   //!
   TBranch        *b_Stats_LPFO_lpfo_d0;   //!
   TBranch        *b_Stats_LPFO_lpfo_d0error;   //!
   TBranch        *b_Stats_LPFO_lpfo_z0;   //!
   TBranch        *b_Stats_LPFO_lpfo_z0error;   //!
   TBranch        *b_Stats_LPFO_lpfo_phi;   //!
   TBranch        *b_Stats_LPFO_lpfo_phierror;   //!
   TBranch        *b_Stats_LPFO_lpfo_omega;   //!
   TBranch        *b_Stats_LPFO_lpfo_omegaerror;   //!
   TBranch        *b_Stats_LPFO_lpfo_tanlambda;   //!
   TBranch        *b_Stats_LPFO_lpfo_tanlambdaerror;   //!
   TBranch        *b_Stats_LPFO_lpfo_chi2;   //!
   TBranch        *b_Stats_LPFO_lpfo_ndf;   //!
   TBranch        *b_Stats_LPFO_lpfo_vtxpt;   //!
   TBranch        *b_Stats_LPFO_lpfo_endpt;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_likelihood;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_eprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_muprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_piprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_kprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_pprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_pid_hprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_likelihood;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_eprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_muprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_piprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_kprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_pprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_hprob;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_e_dedxdist;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_mu_dedxdist;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_pi_dedxdist;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_k_dedxdist;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_p_dedxdist;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_e_lkhood;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_mu_lkhood;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_pi_lkhood;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_k_lkhood;   //!
   TBranch        *b_Stats_LPFO_lpfo_piddedx_p_lkhood;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_p_at_calo;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closest_beta_0ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closest_beta_10ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closest_beta_50ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closest_beta_100ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_fastest_beta_0ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_fastest_beta_10ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_fastest_beta_50ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_fastest_beta_100ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_cylfit_beta_0ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_cylfit_beta_10ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_cylfit_beta_50ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_cylfit_beta_100ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closestfit_beta_0ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closestfit_beta_10ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closestfit_beta_50ps;   //!
   TBranch        *b_Stats_LPFO_lpfo_pidtof_closestfit_beta_100ps;   //!
   TBranch        *b_Data_LPFO_lpfo_config;   //!
   TBranch        *b_Data_LPFO_lpfo_p_mag;   //!
   TBranch        *b_Data_LPFO_lpfo_dEdx_dist_pdg;   //!
   TBranch        *b_Data_LPFO_lpfo_cos;   //!
   TBranch        *b_Data_LPFO_lpfo_qcos;   //!

   LPFOAnalyzer(TTree *tree=0);
   virtual ~LPFOAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef LPFOAnalyzer_cxx
LPFOAnalyzer::LPFOAnalyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../rootfiles/tmp_root/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n002.d_dstm_15162_000.ss.tmp.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../rootfiles/tmp_root/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n002.d_dstm_15162_000.ss.tmp.root");
      }
      f->GetObject("event",tree);

   }
   Init(tree);
}

LPFOAnalyzer::~LPFOAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LPFOAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LPFOAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void LPFOAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eve_valid_lpfo", &eve_valid_lpfo, &b_Event_eve_valid_lpfo);
   fChain->SetBranchAddress("mc_quark_E[2]", mc_quark_E, &b_MC_mc_quark_E);
   fChain->SetBranchAddress("mc_quark_px[2]", mc_quark_px, &b_MC_mc_quark_px);
   fChain->SetBranchAddress("mc_quark_py[2]", mc_quark_py, &b_MC_mc_quark_py);
   fChain->SetBranchAddress("mc_quark_pz[2]", mc_quark_pz, &b_MC_mc_quark_pz);
   fChain->SetBranchAddress("mc_quark_m[2]", mc_quark_m, &b_MC_mc_quark_m);
   fChain->SetBranchAddress("mc_quark_pdg[2]", mc_quark_pdg, &b_MC_mc_quark_pdg);
   fChain->SetBranchAddress("mc_quark_charge[2]", mc_quark_charge, &b_MC_mc_quark_charge);
   fChain->SetBranchAddress("mc_ISR_E[2]", mc_ISR_E, &b_MC_mc_ISR_E);
   fChain->SetBranchAddress("mc_ISR_px[2]", mc_ISR_px, &b_MC_mc_ISR_px);
   fChain->SetBranchAddress("mc_ISR_py[2]", mc_ISR_py, &b_MC_mc_ISR_py);
   fChain->SetBranchAddress("mc_ISR_pz[2]", mc_ISR_pz, &b_MC_mc_ISR_pz);
   fChain->SetBranchAddress("mc_ISR_m[2]", mc_ISR_m, &b_MC_mc_ISR_m);
   fChain->SetBranchAddress("mc_ISR_pdg[2]", mc_ISR_pdg, &b_MC_mc_ISR_pdg);
   fChain->SetBranchAddress("mc_ISR_charge[2]", mc_ISR_charge, &b_MC_mc_ISR_charge);
   fChain->SetBranchAddress("mc_quark_ps_n", &mc_quark_ps_n, &b_MC_mc_quark_ps_n);
   fChain->SetBranchAddress("mc_quark_ps_E[1000]", mc_quark_ps_E, &b_MC_mc_quark_ps_E);
   fChain->SetBranchAddress("mc_quark_ps_px[1000]", mc_quark_ps_px, &b_MC_mc_quark_ps_px);
   fChain->SetBranchAddress("mc_quark_ps_py[1000]", mc_quark_ps_py, &b_MC_mc_quark_ps_py);
   fChain->SetBranchAddress("mc_quark_ps_pz[1000]", mc_quark_ps_pz, &b_MC_mc_quark_ps_pz);
   fChain->SetBranchAddress("mc_quark_ps_m[1000]", mc_quark_ps_m, &b_MC_mc_quark_ps_m);
   fChain->SetBranchAddress("mc_quark_ps_pdg[1000]", mc_quark_ps_pdg, &b_MC_mc_quark_ps_pdg);
   fChain->SetBranchAddress("mc_quark_ps_charge[1000]", mc_quark_ps_charge, &b_MC_mc_quark_ps_charge);
   fChain->SetBranchAddress("mc_quark_ps_y12", &mc_quark_ps_y12, &b_MC_mc_quark_ps_y12);
   fChain->SetBranchAddress("mc_quark_ps_y23", &mc_quark_ps_y23, &b_MC_mc_quark_ps_y23);
   fChain->SetBranchAddress("mc_quark_ps_d12", &mc_quark_ps_d12, &b_MC_mc_quark_ps_d12);
   fChain->SetBranchAddress("mc_quark_ps_d23", &mc_quark_ps_d23, &b_MC_mc_quark_ps_d23);
   fChain->SetBranchAddress("mc_quark_ps_jet_E[2]", mc_quark_ps_jet_E, &b_MC_mc_quark_ps_jet_E);
   fChain->SetBranchAddress("mc_quark_ps_jet_px[2]", mc_quark_ps_jet_px, &b_MC_mc_quark_ps_jet_px);
   fChain->SetBranchAddress("mc_quark_ps_jet_py[2]", mc_quark_ps_jet_py, &b_MC_mc_quark_ps_jet_py);
   fChain->SetBranchAddress("mc_quark_ps_jet_pz[2]", mc_quark_ps_jet_pz, &b_MC_mc_quark_ps_jet_pz);
   fChain->SetBranchAddress("mc_stable_n", &mc_stable_n, &b_MC_mc_stable_n);
   fChain->SetBranchAddress("mc_stable_E[1000]", mc_stable_E, &b_MC_mc_stable_E);
   fChain->SetBranchAddress("mc_stable_px[1000]", mc_stable_px, &b_MC_mc_stable_px);
   fChain->SetBranchAddress("mc_stable_py[1000]", mc_stable_py, &b_MC_mc_stable_py);
   fChain->SetBranchAddress("mc_stable_pz[1000]", mc_stable_pz, &b_MC_mc_stable_pz);
   fChain->SetBranchAddress("mc_stable_m[1000]", mc_stable_m, &b_MC_mc_stable_m);
   fChain->SetBranchAddress("mc_stable_pdg[1000]", mc_stable_pdg, &b_MC_mc_stable_pdg);
   fChain->SetBranchAddress("mc_stable_charge[1000]", mc_stable_charge, &b_MC_mc_stable_charge);
   fChain->SetBranchAddress("mc_stable_isoverlay[1000]", mc_stable_isoverlay, &b_MC_mc_stable_isoverlay);
   fChain->SetBranchAddress("mc_stable_isisr[1000]", mc_stable_isisr, &b_MC_mc_stable_isisr);
   fChain->SetBranchAddress("mc_stable_y12", &mc_stable_y12, &b_MC_mc_stable_y12);
   fChain->SetBranchAddress("mc_stable_y23", &mc_stable_y23, &b_MC_mc_stable_y23);
   fChain->SetBranchAddress("mc_stable_d12", &mc_stable_d12, &b_MC_mc_stable_d12);
   fChain->SetBranchAddress("mc_stable_d23", &mc_stable_d23, &b_MC_mc_stable_d23);
   fChain->SetBranchAddress("mc_stable_jet_E[2]", mc_stable_jet_E, &b_MC_mc_stable_jet_E);
   fChain->SetBranchAddress("mc_stable_jet_px[2]", mc_stable_jet_px, &b_MC_mc_stable_jet_px);
   fChain->SetBranchAddress("mc_stable_jet_py[2]", mc_stable_jet_py, &b_MC_mc_stable_jet_py);
   fChain->SetBranchAddress("mc_stable_jet_pz[2]", mc_stable_jet_pz, &b_MC_mc_stable_jet_pz);
   fChain->SetBranchAddress("lpfo_match[2]", lpfo_match, &b_Stats_LPFO_lpfo_match);
   fChain->SetBranchAddress("lpfo_truejet_pdg[2]", lpfo_truejet_pdg, &b_Stats_LPFO_lpfo_truejet_pdg);
   fChain->SetBranchAddress("lpfo_truejet_type[2]", lpfo_truejet_type, &b_Stats_LPFO_lpfo_truejet_type);
   fChain->SetBranchAddress("lpfo_pdgcheat[2]", lpfo_pdgcheat, &b_Stats_LPFO_lpfo_pdgcheat);
   fChain->SetBranchAddress("lpfo_nparents[2]", lpfo_nparents, &b_Stats_LPFO_lpfo_nparents);
   fChain->SetBranchAddress("lpfo_pdgcheat_parent[2][1000]", lpfo_pdgcheat_parent, &b_Stats_LPFO_lpfo_pdgcheat_parent);
   fChain->SetBranchAddress("lpfo_E[2]", lpfo_E, &b_Stats_LPFO_lpfo_E);
   fChain->SetBranchAddress("lpfo_px[2]", lpfo_px, &b_Stats_LPFO_lpfo_px);
   fChain->SetBranchAddress("lpfo_py[2]", lpfo_py, &b_Stats_LPFO_lpfo_py);
   fChain->SetBranchAddress("lpfo_pz[2]", lpfo_pz, &b_Stats_LPFO_lpfo_pz);
   fChain->SetBranchAddress("lpfo_m[2]", lpfo_m, &b_Stats_LPFO_lpfo_m);
   fChain->SetBranchAddress("lpfo_type[2]", lpfo_type, &b_Stats_LPFO_lpfo_type);
   fChain->SetBranchAddress("lpfo_isoverlay[2]", lpfo_isoverlay, &b_Stats_LPFO_lpfo_isoverlay);
   fChain->SetBranchAddress("lpfo_isisr[2]", lpfo_isisr, &b_Stats_LPFO_lpfo_isisr);
   fChain->SetBranchAddress("lpfo_vtx[2]", lpfo_vtx, &b_Stats_LPFO_lpfo_vtx);
   fChain->SetBranchAddress("lpfo_charge[2]", lpfo_charge, &b_Stats_LPFO_lpfo_charge);
   fChain->SetBranchAddress("lpfo_ntracks[2]", lpfo_ntracks, &b_Stats_LPFO_lpfo_ntracks);
   fChain->SetBranchAddress("lpfo_tpc_hits[2]", lpfo_tpc_hits, &b_Stats_LPFO_lpfo_tpc_hits);
   fChain->SetBranchAddress("lpfo_dedx[2]", lpfo_dedx, &b_Stats_LPFO_lpfo_dedx);
   fChain->SetBranchAddress("lpfo_dedxerror[2]", lpfo_dedxerror, &b_Stats_LPFO_lpfo_dedxerror);
   fChain->SetBranchAddress("lpfo_d0[2]", lpfo_d0, &b_Stats_LPFO_lpfo_d0);
   fChain->SetBranchAddress("lpfo_d0error[2]", lpfo_d0error, &b_Stats_LPFO_lpfo_d0error);
   fChain->SetBranchAddress("lpfo_z0[2]", lpfo_z0, &b_Stats_LPFO_lpfo_z0);
   fChain->SetBranchAddress("lpfo_z0error[2]", lpfo_z0error, &b_Stats_LPFO_lpfo_z0error);
   fChain->SetBranchAddress("lpfo_phi[2]", lpfo_phi, &b_Stats_LPFO_lpfo_phi);
   fChain->SetBranchAddress("lpfo_phierror[2]", lpfo_phierror, &b_Stats_LPFO_lpfo_phierror);
   fChain->SetBranchAddress("lpfo_omega[2]", lpfo_omega, &b_Stats_LPFO_lpfo_omega);
   fChain->SetBranchAddress("lpfo_omegaerror[2]", lpfo_omegaerror, &b_Stats_LPFO_lpfo_omegaerror);
   fChain->SetBranchAddress("lpfo_tanlambda[2]", lpfo_tanlambda, &b_Stats_LPFO_lpfo_tanlambda);
   fChain->SetBranchAddress("lpfo_tanlambdaerror[2]", lpfo_tanlambdaerror, &b_Stats_LPFO_lpfo_tanlambdaerror);
   fChain->SetBranchAddress("lpfo_chi2[2]", lpfo_chi2, &b_Stats_LPFO_lpfo_chi2);
   fChain->SetBranchAddress("lpfo_ndf[2]", lpfo_ndf, &b_Stats_LPFO_lpfo_ndf);
   fChain->SetBranchAddress("lpfo_vtxpt[2][3]", lpfo_vtxpt, &b_Stats_LPFO_lpfo_vtxpt);
   fChain->SetBranchAddress("lpfo_endpt[2][3]", lpfo_endpt, &b_Stats_LPFO_lpfo_endpt);
   fChain->SetBranchAddress("lpfo_pid[2]", lpfo_pid, &b_Stats_LPFO_lpfo_pid);
   fChain->SetBranchAddress("lpfo_pid_likelihood[2]", lpfo_pid_likelihood, &b_Stats_LPFO_lpfo_pid_likelihood);
   fChain->SetBranchAddress("lpfo_pid_eprob[2]", lpfo_pid_eprob, &b_Stats_LPFO_lpfo_pid_eprob);
   fChain->SetBranchAddress("lpfo_pid_muprob[2]", lpfo_pid_muprob, &b_Stats_LPFO_lpfo_pid_muprob);
   fChain->SetBranchAddress("lpfo_pid_piprob[2]", lpfo_pid_piprob, &b_Stats_LPFO_lpfo_pid_piprob);
   fChain->SetBranchAddress("lpfo_pid_kprob[2]", lpfo_pid_kprob, &b_Stats_LPFO_lpfo_pid_kprob);
   fChain->SetBranchAddress("lpfo_pid_pprob[2]", lpfo_pid_pprob, &b_Stats_LPFO_lpfo_pid_pprob);
   fChain->SetBranchAddress("lpfo_pid_hprob[2]", lpfo_pid_hprob, &b_Stats_LPFO_lpfo_pid_hprob);
   fChain->SetBranchAddress("lpfo_piddedx[2]", lpfo_piddedx, &b_Stats_LPFO_lpfo_piddedx);
   fChain->SetBranchAddress("lpfo_piddedx_likelihood[2]", lpfo_piddedx_likelihood, &b_Stats_LPFO_lpfo_piddedx_likelihood);
   fChain->SetBranchAddress("lpfo_piddedx_eprob[2]", lpfo_piddedx_eprob, &b_Stats_LPFO_lpfo_piddedx_eprob);
   fChain->SetBranchAddress("lpfo_piddedx_muprob[2]", lpfo_piddedx_muprob, &b_Stats_LPFO_lpfo_piddedx_muprob);
   fChain->SetBranchAddress("lpfo_piddedx_piprob[2]", lpfo_piddedx_piprob, &b_Stats_LPFO_lpfo_piddedx_piprob);
   fChain->SetBranchAddress("lpfo_piddedx_kprob[2]", lpfo_piddedx_kprob, &b_Stats_LPFO_lpfo_piddedx_kprob);
   fChain->SetBranchAddress("lpfo_piddedx_pprob[2]", lpfo_piddedx_pprob, &b_Stats_LPFO_lpfo_piddedx_pprob);
   fChain->SetBranchAddress("lpfo_piddedx_hprob[2]", lpfo_piddedx_hprob, &b_Stats_LPFO_lpfo_piddedx_hprob);
   fChain->SetBranchAddress("lpfo_piddedx_e_dedxdist[2]", lpfo_piddedx_e_dedxdist, &b_Stats_LPFO_lpfo_piddedx_e_dedxdist);
   fChain->SetBranchAddress("lpfo_piddedx_mu_dedxdist[2]", lpfo_piddedx_mu_dedxdist, &b_Stats_LPFO_lpfo_piddedx_mu_dedxdist);
   fChain->SetBranchAddress("lpfo_piddedx_pi_dedxdist[2]", lpfo_piddedx_pi_dedxdist, &b_Stats_LPFO_lpfo_piddedx_pi_dedxdist);
   fChain->SetBranchAddress("lpfo_piddedx_k_dedxdist[2]", lpfo_piddedx_k_dedxdist, &b_Stats_LPFO_lpfo_piddedx_k_dedxdist);
   fChain->SetBranchAddress("lpfo_piddedx_p_dedxdist[2]", lpfo_piddedx_p_dedxdist, &b_Stats_LPFO_lpfo_piddedx_p_dedxdist);
   fChain->SetBranchAddress("lpfo_piddedx_e_lkhood[2]", lpfo_piddedx_e_lkhood, &b_Stats_LPFO_lpfo_piddedx_e_lkhood);
   fChain->SetBranchAddress("lpfo_piddedx_mu_lkhood[2]", lpfo_piddedx_mu_lkhood, &b_Stats_LPFO_lpfo_piddedx_mu_lkhood);
   fChain->SetBranchAddress("lpfo_piddedx_pi_lkhood[2]", lpfo_piddedx_pi_lkhood, &b_Stats_LPFO_lpfo_piddedx_pi_lkhood);
   fChain->SetBranchAddress("lpfo_piddedx_k_lkhood[2]", lpfo_piddedx_k_lkhood, &b_Stats_LPFO_lpfo_piddedx_k_lkhood);
   fChain->SetBranchAddress("lpfo_piddedx_p_lkhood[2]", lpfo_piddedx_p_lkhood, &b_Stats_LPFO_lpfo_piddedx_p_lkhood);
   fChain->SetBranchAddress("lpfo_pidtof_p_at_calo[2]", lpfo_pidtof_p_at_calo, &b_Stats_LPFO_lpfo_pidtof_p_at_calo);
   fChain->SetBranchAddress("lpfo_pidtof_closest_beta_0ps[2]", lpfo_pidtof_closest_beta_0ps, &b_Stats_LPFO_lpfo_pidtof_closest_beta_0ps);
   fChain->SetBranchAddress("lpfo_pidtof_closest_beta_10ps[2]", lpfo_pidtof_closest_beta_10ps, &b_Stats_LPFO_lpfo_pidtof_closest_beta_10ps);
   fChain->SetBranchAddress("lpfo_pidtof_closest_beta_50ps[2]", lpfo_pidtof_closest_beta_50ps, &b_Stats_LPFO_lpfo_pidtof_closest_beta_50ps);
   fChain->SetBranchAddress("lpfo_pidtof_closest_beta_100ps[2]", lpfo_pidtof_closest_beta_100ps, &b_Stats_LPFO_lpfo_pidtof_closest_beta_100ps);
   fChain->SetBranchAddress("lpfo_pidtof_fastest_beta_0ps[2]", lpfo_pidtof_fastest_beta_0ps, &b_Stats_LPFO_lpfo_pidtof_fastest_beta_0ps);
   fChain->SetBranchAddress("lpfo_pidtof_fastest_beta_10ps[2]", lpfo_pidtof_fastest_beta_10ps, &b_Stats_LPFO_lpfo_pidtof_fastest_beta_10ps);
   fChain->SetBranchAddress("lpfo_pidtof_fastest_beta_50ps[2]", lpfo_pidtof_fastest_beta_50ps, &b_Stats_LPFO_lpfo_pidtof_fastest_beta_50ps);
   fChain->SetBranchAddress("lpfo_pidtof_fastest_beta_100ps[2]", lpfo_pidtof_fastest_beta_100ps, &b_Stats_LPFO_lpfo_pidtof_fastest_beta_100ps);
   fChain->SetBranchAddress("lpfo_pidtof_cylfit_beta_0ps[2]", lpfo_pidtof_cylfit_beta_0ps, &b_Stats_LPFO_lpfo_pidtof_cylfit_beta_0ps);
   fChain->SetBranchAddress("lpfo_pidtof_cylfit_beta_10ps[2]", lpfo_pidtof_cylfit_beta_10ps, &b_Stats_LPFO_lpfo_pidtof_cylfit_beta_10ps);
   fChain->SetBranchAddress("lpfo_pidtof_cylfit_beta_50ps[2]", lpfo_pidtof_cylfit_beta_50ps, &b_Stats_LPFO_lpfo_pidtof_cylfit_beta_50ps);
   fChain->SetBranchAddress("lpfo_pidtof_cylfit_beta_100ps[2]", lpfo_pidtof_cylfit_beta_100ps, &b_Stats_LPFO_lpfo_pidtof_cylfit_beta_100ps);
   fChain->SetBranchAddress("lpfo_pidtof_closestfit_beta_0ps[2]", lpfo_pidtof_closestfit_beta_0ps, &b_Stats_LPFO_lpfo_pidtof_closestfit_beta_0ps);
   fChain->SetBranchAddress("lpfo_pidtof_closestfit_beta_10ps[2]", lpfo_pidtof_closestfit_beta_10ps, &b_Stats_LPFO_lpfo_pidtof_closestfit_beta_10ps);
   fChain->SetBranchAddress("lpfo_pidtof_closestfit_beta_50ps[2]", lpfo_pidtof_closestfit_beta_50ps, &b_Stats_LPFO_lpfo_pidtof_closestfit_beta_50ps);
   fChain->SetBranchAddress("lpfo_pidtof_closestfit_beta_100ps[2]", lpfo_pidtof_closestfit_beta_100ps, &b_Stats_LPFO_lpfo_pidtof_closestfit_beta_100ps);
   fChain->SetBranchAddress("lpfo_config", &lpfo_config, &b_Data_LPFO_lpfo_config);
   fChain->SetBranchAddress("lpfo_p_mag[2]", lpfo_p_mag, &b_Data_LPFO_lpfo_p_mag);
   fChain->SetBranchAddress("lpfo_dEdx_dist_pdg[2]", lpfo_dEdx_dist_pdg, &b_Data_LPFO_lpfo_dEdx_dist_pdg);
   fChain->SetBranchAddress("lpfo_cos[2]", lpfo_cos, &b_Data_LPFO_lpfo_cos);
   fChain->SetBranchAddress("lpfo_qcos[2]", lpfo_qcos, &b_Data_LPFO_lpfo_qcos);
   Notify();
}

Bool_t LPFOAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LPFOAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LPFOAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LPFOAnalyzer_cxx
