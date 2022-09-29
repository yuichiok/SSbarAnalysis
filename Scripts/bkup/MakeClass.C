#include "TROOT.h"
#include "TFile.h"

int MakeClass(){

  TString filename = "../rootfiles/tmp_root/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n002.d_dstm_15162_000.ss.tmp.root";
  TFile *file = new TFile(filename);

  TString treename = "event";
  TTree* tree = (TTree*)file->Get( treename ) ;

  file->ls();

  tree->MakeClass("LPFOAnalyzer");

  std::cout << "processed..." << std::endl;

  return 0;


}