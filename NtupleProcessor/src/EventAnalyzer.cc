
/*------------------------------------------------------------------------------
EventAnalyzer.cpp
 Created : 2022-09-05  okugawa
------------------------------------------------------------------------------*/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>

#include "EventAnalyzer.hh"
#include "TreeReader.hh"
#include "PFOTools.hh"
#include "VectorTools.hh"

using std::cout;   using std::endl;
typedef unsigned int Index;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

EventAnalyzer::EventAnalyzer(TString input, TString fnac, TString o)
: options(o), _config(fnac)
{
    _fs.SetNames(input.Data());
    _anCfg.SetConfig(_config);
    patEventsAnalyzed = 0;
    entriesInNtuple   = 0;
}

Bool_t EventAnalyzer::InitReadTree(TTree* tree)
{
  // Maps TTree to class' variables.
  // TO DO: Implement check for correct mapping, return result?
  //    - Set up exception handling for negative result.

    entriesInNtuple = tree->GetEntries();

  // Set branch addresses and branch pointers
    if (!tree) return false;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

  // Read Tree
    TreeReader reader;
    reader.InitializeMCReadTree(fChain, _mc, _branch);
    reader.InitializeJetReadTree(fChain, _jet, _branch);
    reader.InitializePFOReadTree(fChain, _pfo, _branch);
    reader.InitializeVTXReadTree(fChain,_vtx, _branch);

    Notify();

    return true;

}

void EventAnalyzer::CreateFile()
{
  // Write Tree
    _hfilename = TString(_fs.GetOutName_withPath());
    _hfile = new TFile( _hfilename, "RECREATE", _hfilename ) ;

}

void EventAnalyzer::InitWriteTree()
{
  // Initialize Write Tree
    _hTree     = new TTree( "data", "tree" );
    writer.InitializeDataTree(_hTree,_data);
}

void EventAnalyzer::InitHists()
{
    _hm.InitializeHists();
}

void EventAnalyzer::WriteFile()
{
  // Write Tree
    _hm.WriteLists(_hfile);
    _hTree->Write();
    // _hfile->Write();
    _hfile->Close();

}

void EventAnalyzer::AnalyzeReco(Long64_t entry)
{
  // MC, PFO Analysis
  PFOTools mct( &_mc, _config );
  PFOTools pfot( &_mc, &_pfo, _config );

  // cout << "===== " << entry + 1 << " =====" << endl;

  Bool_t low_ctag = true;
  for(Int_t ijet=0; ijet<2; ijet++){
    _hm.h1_K_reco[_hm.ctag]->Fill(_jet.jet_ctag[ijet]);
    if(_jet.jet_ctag[ijet] > 0.8) {
        low_ctag = false;
    }
  }

  cout<<_jet.jet_nvtx[0]<<" "<<_jet.jet_nvtx[1]<<endl;
  CCbarAnalysis(pfot, low_ctag);


  _hm.h_tagging[_hm.jets_info]->Fill(0);

  if(_jet.jet_npfo[0] >= 1 and _jet.jet_npfo[1] >= 1) {
    _hm.h_tagging[_hm.jets_info]->Fill(6);
  }
  if(_jet.jet_npfo[0] == 0 and _jet.jet_npfo[1] == 0) {
    _hm.h_tagging[_hm.jets_info]->Fill(7);
  }
  // cout<<_jet.jet_nvtx[0]<<"  "<<_jet.jet_nvtx[1]<<endl;
  // with SV
  if (_jet.jet_nvtx[0] >= 1) {
    _hm.h_tagging[_hm.s_ctag]->Fill(_jet.jet_ctag[0]);
    _hm.h_tagging[_hm.s_btag]->Fill(_jet.jet_btag[0]);
    _hm.h_tagging[_hm.jets_info]->Fill(1);
    _hm.h_tagging[_hm.jets_info]->Fill(4);
  }
  if (_jet.jet_nvtx[1] >= 1) {
    _hm.h_tagging[_hm.s_ctag]->Fill(_jet.jet_ctag[1]);
    _hm.h_tagging[_hm.s_btag]->Fill(_jet.jet_btag[1]);
    _hm.h_tagging[_hm.jets_info]->Fill(2);
    _hm.h_tagging[_hm.jets_info]->Fill(4);
  }

  // w/o SV
  if (_jet.jet_nvtx[0] == 0) {
    _hm.h_tagging[_hm.s_ctag]->Fill(_jet.jet_ctag[0]);
    _hm.h_tagging[_hm.s_btag]->Fill(_jet.jet_btag[0]);
    _hm.h_tagging[_hm.jets_info]->Fill(1);
    _hm.h_tagging[_hm.jets_info]->Fill(5);
  }
  if (_jet.jet_nvtx[1] == 0) {
    _hm.h_tagging[_hm.p_ctag]->Fill(_jet.jet_ctag[1]);
    _hm.h_tagging[_hm.p_btag]->Fill(_jet.jet_btag[1]);
    _hm.h_tagging[_hm.jets_info]->Fill(2);
    _hm.h_tagging[_hm.jets_info]->Fill(5);
  }

  // w/o Vertex
  if (_jet.jet_nvtx[0] == -1) {
    _hm.h_tagging[_hm.t_ctag]->Fill(_jet.jet_ctag[0]);
    _hm.h_tagging[_hm.t_btag]->Fill(_jet.jet_btag[0]);
    _hm.h_tagging[_hm.jets_info]->Fill(1);
    _hm.h_tagging[_hm.jets_info]->Fill(3);
  }
  if (_jet.jet_nvtx[1] == -1) {
    _hm.h_tagging[_hm.t_ctag]->Fill(_jet.jet_ctag[1]);
    _hm.h_tagging[_hm.t_btag]->Fill(_jet.jet_btag[1]);
    _hm.h_tagging[_hm.jets_info]->Fill(2);
    _hm.h_tagging[_hm.jets_info]->Fill(3);
  }
  
  // _hTree->Fill();
  ClearStructs();

}

void EventAnalyzer::ClearStructs()
{
  _eve        = {};
  _stats_lpfo = {};
  _data_lpfo  = {};

  _data       = {};
}

Bool_t EventAnalyzer::Select(Selector sel)
{ // Evaluates the class' list of event selection criteria

  /*
  Must initialize
    - Float_t MINP_CUT (= 20 GeV)
    - TString PROCESS  (= "SS")
    - TString FILE_OUT (?)
  */

    VectorTools mcvt[2];
    VectorTools jetvt[2];
    for (int i=0; i<2; i++){
      mcvt[i].SetCoordinates(_mc.mc_quark_px[i],_mc.mc_quark_py[i],_mc.mc_quark_pz[i],_mc.mc_quark_E[i]);
      jetvt[i].SetCoordinates(_jet.jet_px[i],_jet.jet_py[i],_jet.jet_pz[i],_jet.jet_E[i]);
    }

    vector<Bool_t> CutTrigger;

    switch (sel){
      case kQQ:
        CutTrigger.push_back( GenPairPicker( _mc.mc_quark_pdg[0], _anCfg.gen_quarks ) );
        break;
      case kMC:
        CutTrigger.push_back( Cut_ESum( mcvt ) );
        CutTrigger.push_back( Cut_ACol( mcvt ) );
        CutTrigger.push_back( Cut_ISR ( mcvt ) );
        break;
      case kReco:
        CutTrigger.push_back( Cut_ESum( jetvt ) );
        CutTrigger.push_back( Cut_ACol( jetvt ) );
        CutTrigger.push_back( Cut_ISR ( jetvt ) );
        break;

      default:
        break;
    }

  // Options

    for (auto trigger : CutTrigger ){
      if (!trigger) { return false; }
    }

    return true;

}

Bool_t EventAnalyzer::GenPairPicker ( Float_t mc_particle, std::vector<int> input_gen )
{
  for ( auto igen : input_gen ){
    if( fabs(mc_particle) == igen ) return true;
  }
  
  return false;
}

Bool_t EventAnalyzer::Cut_ESum ( VectorTools v[2] )
{
  Float_t SumE = v[0].v4().E() + v[1].v4().E();

  return (SumE > 220);
}

Bool_t EventAnalyzer::Cut_ACol ( VectorTools v[2] )
{
  Float_t cosacol = std::cos( VectorTools::GetThetaBetween( v[0].v3(), v[1].v3() ) );

  return (cosacol > 0.95);
}

Bool_t EventAnalyzer::Cut_ISR ( VectorTools v[2] )
{
  using namespace ROOT::Math;

	if (v[0].v4().E() < 0.5 || v[1].v4().E() < 0.5)
		return false;

  Float_t abscos[2] = {0};
  for (int i=0; i < 2; i++)
  {
    abscos[i] = fabs( std::cos( v[i].v3().Theta() ) );
  }

	Float_t ssmass  = (v[0].v4() + v[1].v4()).M();
  Float_t sinacol = std::sin( VectorTools::GetThetaBetween( v[0].v3(), v[1].v3() ) );

	Float_t Kv = 250. * sinacol / (sinacol + sqrt(1 - abscos[0] * abscos[0]) + sqrt(1 - abscos[1] * abscos[1]));
  // Float_t K[2] = {0};
  //         K[0] = jet_vec[0].v3().R() * sinacol / sqrt(1 - abscos[1] * abscos[1]);
  //         K[1] = jet_vec[1].v3().R() * sinacol / sqrt(1 - abscos[0] * abscos[0]);

  return (Kv < 35 && ssmass > 130);

}

Bool_t EventAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventAnalyzer::CCbarAnalysis(PFOTools pfot, Bool_t low_ctag)
{
  std::vector<PFO_Info> PFO_Collection = pfot.Get_Valid_PFOs();

  Int_t nprong[2][2] = {{0,0},{0,0}};

  for ( long unsigned int i=0; i < PFO_Collection.size(); i++ )
  {
    PFO_Info ipfo = PFO_Collection.at(i);
    Int_t imatch = ipfo.pfo_match;

    if(ipfo.pfo_vtx==0){
      nprong[imatch][0]++;
    }else if(ipfo.pfo_vtx==1){
      nprong[imatch][1]++;
    }

  }

  for ( long unsigned int i=0; i < PFO_Collection.size(); i++ )
  {
    PFO_Info ipfo = PFO_Collection.at(i);
    Int_t  imatch = ipfo.pfo_match;

    //---------------Own-----------------//
    if(ipfo.p_mag>5 and ipfo.pfo_d0!=0 and ipfo.pfo_z0!=0)
    {
      if( nprong[imatch][0]==1 && ipfo.pfo_vtx==0 ){
        _hm.h_PS[_hm.d0_P_single]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h_PS[_hm.z0_P_single]->Fill(TMath::Abs(ipfo.pfo_z0));
      }else if( nprong[imatch][1]==1 && ipfo.pfo_vtx==0 ){
        _hm.h_PS[_hm.d0_S_single]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h_PS[_hm.z0_S_single]->Fill(TMath::Abs(ipfo.pfo_z0));
      }else if ( nprong[imatch][0]>1 && ipfo.pfo_vtx==1 ){
        _hm.h_PS[_hm.d0_P_mult]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h_PS[_hm.z0_P_mult]->Fill(TMath::Abs(ipfo.pfo_z0));
      }else if ( nprong[imatch][1]>1 && ipfo.pfo_vtx==1 ){
        _hm.h_PS[_hm.d0_S_mult]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h_PS[_hm.z0_S_mult]->Fill(TMath::Abs(ipfo.pfo_z0));
      }

      if(abs(ipfo.pfo_pdgcheat)==321)
      {
        _hm.h_general[_hm.n_K_ecal]->Fill(ipfo.pfo_charge);
      }
      // < new >
      // charged kaons <=> abs(ipfo.pfo_pdgcheat)==321
      // charged kaons of A category <=> abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==0
      // charged kaons of A category which have a multiple prong <=> abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==0 and nprong[imatch][0]>1
      // => have only primary vertex      
      if(abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==0 and nprong[imatch][0]>1)
      {
        _hm.h1_K_reco[_hm.d0_K_reco_primary]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h1_K_reco[_hm.d0_sigma_K_reco_primary]->Fill(TMath::Sqrt(ipfo.pfo_d0error));
        _hm.h1_K_reco[_hm.d0_sigma_d0_K_reco_primary]->Fill(TMath::Abs(ipfo.pfo_d0)/TMath::Sqrt(ipfo.pfo_d0error));
        
        _hm.h1_K_reco[_hm.z0_K_reco_primary]->Fill(TMath::Abs(ipfo.pfo_z0));
        _hm.h1_K_reco[_hm.z0_sigma_K_reco_primary]->Fill(TMath::Sqrt(ipfo.pfo_z0error));
        _hm.h1_K_reco[_hm.z0_sigma_z0_K_reco_primary]->Fill(TMath::Abs(ipfo.pfo_z0)/TMath::Sqrt(ipfo.pfo_z0error));
      }

      // charged kaons <=> abs(ipfo.pfo_pdgcheat)==321
      // charged kaons of B or C category <=> abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==1
      // nprong[imatch][1]>1 => not a single prong => not a category C => charged kaons of B category
      // charged kaons of B category <=> abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==1 and nprong[imatch][1]>1
      // => have primary and secondary vertices
      if(abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==1 and nprong[imatch][1]>1)
      {
        _hm.h1_K_reco[_hm.d0_K_reco_secondary]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h1_K_reco[_hm.d0_sigma_K_reco_secondary]->Fill(TMath::Sqrt(ipfo.pfo_d0error));
        _hm.h1_K_reco[_hm.d0_sigma_d0_K_reco_secondary]->Fill(TMath::Abs(ipfo.pfo_d0)/TMath::Sqrt(ipfo.pfo_d0error));
        
        _hm.h1_K_reco[_hm.z0_K_reco_secondary]->Fill(TMath::Abs(ipfo.pfo_z0));
        _hm.h1_K_reco[_hm.z0_sigma_K_reco_secondary]->Fill(TMath::Sqrt(ipfo.pfo_z0error));
        _hm.h1_K_reco[_hm.z0_sigma_z0_K_reco_secondary]->Fill(TMath::Abs(ipfo.pfo_z0)/TMath::Sqrt(ipfo.pfo_z0error));
      }

      if(abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==1 and nprong[imatch][1]==1)
      {
        _hm.h1_K_reco[_hm.d0_K_reco_pseudo]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h1_K_reco[_hm.d0_sigma_K_reco_pseudo]->Fill(TMath::Sqrt(ipfo.pfo_d0error));
        _hm.h1_K_reco[_hm.d0_sigma_d0_K_reco_pseudo]->Fill(TMath::Abs(ipfo.pfo_d0)/TMath::Sqrt(ipfo.pfo_d0error));
        
        _hm.h1_K_reco[_hm.z0_K_reco_pseudo]->Fill(TMath::Abs(ipfo.pfo_z0));
        _hm.h1_K_reco[_hm.z0_sigma_K_reco_pseudo]->Fill(TMath::Sqrt(ipfo.pfo_z0error));
        _hm.h1_K_reco[_hm.z0_sigma_z0_K_reco_pseudo]->Fill(TMath::Abs(ipfo.pfo_z0)/TMath::Sqrt(ipfo.pfo_z0error));
      }



      if(low_ctag and abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==0 and nprong[imatch][0]>1)
      {

        _hm.h1_K_reco[_hm.d0_K_reco_primary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h1_K_reco[_hm.d0_sigma_K_reco_primary_ctag_cut]->Fill(TMath::Sqrt(ipfo.pfo_d0error));
        _hm.h1_K_reco[_hm.d0_sigma_d0_K_reco_primary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_d0)/TMath::Sqrt(ipfo.pfo_d0error));
        
        _hm.h1_K_reco[_hm.z0_K_reco_primary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_z0));
        _hm.h1_K_reco[_hm.z0_sigma_K_reco_primary_ctag_cut]->Fill(TMath::Sqrt(ipfo.pfo_z0error));
        _hm.h1_K_reco[_hm.z0_sigma_z0_K_reco_primary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_z0)/TMath::Sqrt(ipfo.pfo_z0error));
      }

      if(low_ctag and abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==1 and nprong[imatch][1]>1)
      {
        _hm.h1_K_reco[_hm.d0_K_reco_secondary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h1_K_reco[_hm.d0_sigma_K_reco_secondary_ctag_cut]->Fill(TMath::Sqrt(ipfo.pfo_d0error));
        _hm.h1_K_reco[_hm.d0_sigma_d0_K_reco_secondary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_d0)/TMath::Sqrt(ipfo.pfo_d0error));
        
        _hm.h1_K_reco[_hm.z0_K_reco_secondary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_z0));
        _hm.h1_K_reco[_hm.z0_sigma_K_reco_secondary_ctag_cut]->Fill(TMath::Sqrt(ipfo.pfo_z0error));
        _hm.h1_K_reco[_hm.z0_sigma_z0_K_reco_secondary_ctag_cut]->Fill(TMath::Abs(ipfo.pfo_z0)/TMath::Sqrt(ipfo.pfo_z0error));
      }
      
      if(abs(ipfo.pfo_pdgcheat)==321 and ipfo.pfo_vtx==-1)
      {
        _hm.h1_K_reco[_hm.d0_K_reco_garbage]->Fill(TMath::Abs(ipfo.pfo_d0));
        _hm.h1_K_reco[_hm.d0_sigma_K_reco_garbage]->Fill(TMath::Sqrt(ipfo.pfo_d0error));
        _hm.h1_K_reco[_hm.d0_sigma_d0_K_reco_garbage]->Fill(TMath::Abs(ipfo.pfo_d0)/TMath::Sqrt(ipfo.pfo_d0error));
        
        _hm.h1_K_reco[_hm.z0_K_reco_garbage]->Fill(TMath::Abs(ipfo.pfo_z0));
        _hm.h1_K_reco[_hm.z0_sigma_K_reco_garbage]->Fill(TMath::Sqrt(ipfo.pfo_z0error));
        _hm.h1_K_reco[_hm.z0_sigma_z0_K_reco_garbage]->Fill(TMath::Abs(ipfo.pfo_z0)/TMath::Sqrt(ipfo.pfo_z0error));
      }
      // </ new >
      if(abs(ipfo.pfo_pdgcheat)==321)
      {
        _hm.h1_K_reco[_hm.pmag_K_reco]->Fill(ipfo.p_mag);
        _hm.h1_K_reco[_hm.cos_theta_K_reco]->Fill(ipfo.qcos);
      }

    }

  }
    
}