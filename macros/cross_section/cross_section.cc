#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::array; using std::unordered_map;

TString LPFO_mode = "K";
TString ichiral = "eL.pR";
// TString ichiral = "eR.pL";

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
array<TString,4> processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"};
array<TString,6> qqbars    = {"dd", "uu", "ss", "cc", "bb", "rr"};

unordered_map<pair<TString,TString>,pair<Int_t,Int_t>, hash_pair> production = {
    {{"P2f_z_h", "eL.pR"}, {500010,4994}},
    {{"P2f_z_h", "eR.pL"}, {500012,4994}},
    {{"P4f_ww_h", "eL.pR"}, {500066,4996}},
    {{"P4f_ww_h", "eR.pL"}, {500068,5116}},
    {{"P4f_zz_h", "eL.pR"}, {500062,5052}},
    {{"P4f_zz_h", "eR.pL"}, {500064,5109}},
    {{"Pqqh", "eL.pR"}, {402011,1457}},
    {{"Pqqh", "eR.pL"}, {402012,2278}},
};

void latexTransformation(unordered_map<TString, unordered_map<TString, unordered_map< TString, Float_t > > > cross_sections)
{
  for( auto process : processes ){
    if (process!="P2f_z_h") continue;

    for ( auto iqq : qqbars ){
      cout << iqq;
      for ( auto chiral : chirals ){
        cout << " & " << cross_sections[process][chiral][iqq];
      }
      cout << "  \\\\" << endl;

    }

  }

}

void cross_section()
{
  try
  {
    unordered_map<TString, unordered_map<TString, TFile*> > file_map;
    unordered_map<TString, unordered_map<TString, unordered_map< TString, unordered_map<TString, TH1F*> > > > hmap;
    unordered_map<TString, unordered_map<TString, unordered_map< TString, Float_t > > > cross_sections;
    unordered_map<TString, unordered_map<TString, unordered_map< TString, Int_t > > > Nreco;
    for( auto process : processes ){
      for( auto chiral : chirals ){
        Int_t processID = production.at({process,chiral}).first;
        cout << process << " " << chiral << " " << processID << endl;
        TString filename = inputDir + "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I" + processID + "." + process + "." + chiral + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root";
        TFile *file = new TFile(filename,"READ");
        if( !file->IsOpen() ) throw std::runtime_error("File not found");

        if( process=="P2f_z_h" ){
          for( auto iqq : qqbars ){
            TH1F *h = (TH1F*) file->Get( iqq + "/gen/h_" + iqq + "_qcos" );
            TH1F *hreco = (TH1F*) file->Get( iqq + "/cos/h_" + iqq + "_" + LPFO_mode + "_qcos" );
            Float_t luminosity = production.at({process,chiral}).second;
            Float_t cross_section = h->Integral() / luminosity;
            cout << "qq = " << iqq << ", cross_section = " << cross_section << endl;
            cross_sections[process][chiral][iqq] = cross_section;
            Nreco[process][chiral][iqq] = hreco->GetEntries();
          }
        }else{
          TString iqq = "bg";
          TH1F *hreco = (TH1F*) file->Get( iqq + "/cos/h_" + iqq + "_" + LPFO_mode + "_qcos" );
          Nreco[process][chiral][iqq] = hreco->GetEntries();
        }
      }
    }

    latexTransformation(cross_sections);

    cout << Nreco["P2f_z_h"][ichiral]["ss"] << endl;
    
    Int_t Nsig = 0;
    Int_t Nbg  = 0;

    for( auto process : processes ){
      for( auto iqq : qqbars ){
        if(process=="P2f_z_h"){
          if(iqq=="ss"){
            Nsig += Nreco[process][ichiral][iqq];
          }else{
            Nbg += Nreco[process][ichiral][iqq];
          }
        }else{
            Nbg += Nreco[process][ichiral]["bg"];
        }
      }
    }
    cout << "N signal = " << Nsig << ", sigma = " << 1.0 / sqrt(Nsig) << "\n";
    cout << "N bg     = " << Nbg  << ", sigma = " << 1.0 / sqrt(Nbg) << "\n";

    Float_t sigma_sig  = 1.0 / sqrt(Nsig);
    Float_t sigma_stat = sqrt( sigma_sig*sigma_sig + 0.001*0.001 );
    cout << "sigma_sys = " << sigma_stat << "\n";

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}