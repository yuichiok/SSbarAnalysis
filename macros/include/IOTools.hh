#include <unordered_map>
#include "../../SSbarLibrary/include/MapTString.hh"

using std::unordered_map;

class IOTools
{
  public:
    IOTools( TFile * infile );
    ~IOTools();
  private:
    TFile * file;

    // Declear histograms
    // h1 hist
      vector<TString> hcos_gen_name = {"cos","qcos"};
      unordered_map< TString, TH1F* > h1_gen_cos;        // [hist]

      vector<TString> hcos_name = {"cos","qcos","scos","acc_cos","rej_cos"};
      unordered_map< TString, unordered_map< TString, TH1F* > > h1_cos;        // [LPFO][hist]

      vector<TString> hres_name = {"gen_N_cos","reco_N_cos","N_corr_cos"};
      unordered_map< TString, unordered_map< TString, TH1F* > > h1_resolution; // [LPFO][hist]

      // efficiency plots
      vector<TString> gen_reco  = {"gen","reco"};
      vector<TString> heff_name = {"momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
      vector<TString> heff_dedx_name = {"dEdx_p","dEdx_error_cos","dEdx_dist_cos"};
      unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_cos_eff;  // [GenReco][LPFO][cut]
      unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > > h2_dEdx_eff;  // [GenReco][LPFO][TruthID][cut][hist]

    // h2 hist
      vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};
      unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]


};

IOTools::IOTools( TFile * infile )
: file(infile)
{
}

IOTools::~IOTools()
{
}
