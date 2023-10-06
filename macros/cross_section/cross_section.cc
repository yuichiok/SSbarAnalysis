#include <iostream>
#include <vector>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::array; using std::unordered_map;

TString inputDir = "../../rootfiles/merged/";
array<TString,2> chirals   = {"eL.pR", "eR.pL"};
array<TString,4> processes = {"P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"};
array<TString,6> qqbars    = {"dd", "uu", "ss", "cc", "bb", "rr"};

unordered_map<pair<TString,TString>,pair<Int_t,Int_t>, hash_pair> production = {
    {{"P2f_z_h", "eL.pR"}, {500010,4994}},
    {{"P2f_z_h", "eR.pL"}, {500012,4994}},
    {{"P4f_ww_h", "eL.pR"}, {500066,4996}},
    {{"P4f_ww_h", "eR.pL"}, {500068,5116}},
    {{"P4f_zz_h", "eL.pR"}, {500062,5052}},
    {{"P4f_zz_h", "eR.pL"}, {500064,5109}},
    {{"Pe1e1h", "eL.pL"}, {402013,801943}},
    {{"Pe1e1h", "eL.pR"}, {402001,28294}},
    {{"Pe1e1h", "eR.pL"}, {402002,44887}},
    {{"Pe1e1h", "eR.pR"}, {402014,800018}}
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
            Float_t luminosity = production.at({process,chiral}).second;
            Float_t cross_section = h->Integral() / luminosity;
            cout << "qq = " << iqq << ", cross_section = " << cross_section << endl;
            cross_sections[process][chiral][iqq] = cross_section;
          }
        }
      }
    }

    latexTransformation(cross_sections);

    




  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}