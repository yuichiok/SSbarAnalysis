#include <iostream>
#include <vector>
#include <map>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString prod_mode = "uu";
TString chiral    = "eL.pR";
TString LPFO_mode = "Pi";
// TString LPFO_mode = "K";
// TString qq[3] = {"uu","dd","ss"};
TString qq[5] = {"dd","uu","ss","cc","bb"};
// TString qq[2] = {"dd","uu"};
// TString prod_modes[1] = {"ss"};

vector<TString> gen_reco  = {"gen","reco"};
vector<TString> PFO_mode  = {"K","Pi"};
// vector<TString> heff_name = {"nocut","momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
vector<TString> heff_name = {"nocut", "btag", "ctag", "nvtx", "momentum", "LPFOacol", "offset", "PID", "SPFO", "charge"};

array<TString,2> chirals   = {"eL.pR", "eR.pL"};
array<TString,1> processes = {"P2f_z_h"};


TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root","READ");
// TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500012.P2f_z_h.eR.pL.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root","READ");

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

void BinNormal(TH1F *h)
{
  const Int_t nbins = h->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++){
    Float_t binc = h->GetBinContent(ibin);
    Float_t binw = h->GetBinWidth(ibin);
    Float_t bin_normal = binc / binw;
    h->SetBinContent(ibin,bin_normal);
  }

}

void Normalize(TH1F *h)
{
  // h->Scale( 1.0 / h->GetEntries() );
  const Int_t nbins = h->GetNbinsX();
  Int_t nbins4 = nbins / 4;
  Int_t int_high = (nbins / 2) + nbins4;
  Int_t int_low  = (nbins / 2 + 1) - nbins4;
  h->Scale( 1.0 / h->Integral(int_low,int_high) );
}

void Normalize2Gen(TH1F *h, TH1F *h_gen)
{
	double intCosReco = h->Integral(20,80);
	double intCosGen  = h_gen->Integral(20,80);

  cout << "intCosReco = " << intCosReco << ", intCosGen = " << intCosGen << "\n";

  h->Scale( intCosGen / intCosReco );
}

TH1F* plotEfficiency(TH1F *h_num, TH1F *h_denom)
{
  gStyle->SetOptStat(0);

  cout << "h_num->GetName() = " << h_num->GetName() << "\n";

  Int_t Nnum = h_num->GetEntries();
  Int_t Ndenom = h_denom->GetEntries();
  cout << "Efficiency = " << (Float_t) Nnum / (Float_t) Ndenom  << "\t Nnum = " << Nnum << ", Ndenom = " << Ndenom << "\n";

  TH1F *h_eff = (TH1F*) h_num->Clone();
  h_eff->Divide(h_denom);
  h_eff->GetYaxis()->SetRangeUser(0,1);

  return h_eff;
}

void SaveHists(TCanvas *c, TH1F *ih)
{
  TString printname = "c_" + prod_mode + "_" + (TString)ih->GetName() + ".png";
  c->Print("~/Desktop/" + printname);
}

void PrintEfficiency(TFile *file, vector<TH1F*> hvec)
{
  if (!file->IsOpen()) return;
  TH1F *h_gen_q_qcos     = (TH1F*) file->Get(prod_mode + "/gen/h_" + prod_mode + "_qcos");
  Int_t n_gen_events     = h_gen_q_qcos->GetEntries();
  cout << "name,nevents,efficiency\n";
  cout << "gen," << n_gen_events << ",-\n";
  for ( auto ih : hvec )
  {
    Int_t n_reco_events = ih->GetEntries();
    cout << ih->GetName() << ",";
    cout << ih->GetEntries() << "\n";
  }

}

void calcEff()
{
  unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > > h1_cos_eff;  // [qq][GenReco][LPFO][hist]

  for ( auto igenreco : gen_reco ){
    for ( auto i_lmode : PFO_mode ){
      for ( auto ih : heff_name ){
        for( auto iqq : qq ){
          TString dir_name = "/efficiency/";
          TString hname = "h_" + iqq + "_" + igenreco + "_" + i_lmode + "_" + ih;
          TH1F *h = (TH1F*) file->Get(iqq + dir_name + hname);
          h1_cos_eff[iqq][igenreco][i_lmode][ih] = h;
        }
      }
    }
  }

  unordered_map< TString, unordered_map< TString, Float_t > > eff_map; // [qq][cut]
  for ( auto iqq : qq ){
    Int_t total = 0;
    for ( auto ih : heff_name ){

      Int_t nevt = h1_cos_eff[iqq]["reco"][LPFO_mode][ih]->GetEntries();
      if(ih=="nocut"){
        eff_map[iqq][ih] = nevt;
        total = nevt;
      }else{
        Float_t eff = (Float_t) nevt / (Float_t) total;
        eff_map[iqq][ih] = eff;
      }

    } // end of heff

  } // end of iqq

  // output efficiencies

  int cutno = 0;
  for (auto ih : heff_name ){
    if(!cutno){
      cout << "       & ";
    }else{
      cout << "+ Cut " << cutno << " & ";
    }
    for (auto iqq : qq){
      if(ih=="nocut"){
        // cout << std::setprecision (3) << iqq << " 100\\% (" << eff_map[iqq][ih] << ") & ";
        cout <<  iqq << " 100\\% (" << eff_map[iqq][ih] / 2.0 << ") & ";
      }else{
        cout << std::setprecision (3) << eff_map[iqq][ih] * 100. << "\\% & "; 
      }
    }
    cout << "\\\\" << endl;
    cutno++;
  }

}

// DO IT WHEN THERE IS TIME.
// unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > getCosEff(TFile *ifile)
// {
//   unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_cos_eff;  // [GenReco][LPFO][hist]

//   for ( auto igenreco : gen_reco ){
//     for ( auto i_lmode : PFO_mode ){
//       for ( auto ih : heff_name ){
//         Int_t tmp = 0;
//         for( auto iqq : qq ){
//           TString dir_name = "/efficiency/";
//           TString hname = "h_" + iqq + "_" + igenreco + "_" + i_lmode + "_" + ih;
//           TH1F *h = (TH1F*) ifile->Get(iqq + dir_name + hname);
//           if(tmp==0) {
//             TH1F *htmp = (TH1F*) h->Clone();
//             htmp->SetName("h_" + igenreco + "_" + i_lmode + "_" + ih);
//             h1_cos_eff[igenreco][i_lmode][ih] = htmp;
//           }else{
//             h1_cos_eff[igenreco][i_lmode][ih]->Add(h);
//           }
//           tmp++;
//         }
//       }
//     }
//   }

// }

// void eff_chiral()
// {
//   TString inputDir = "../../rootfiles/merged/";
//   unordered_map<TString, unordered_map<TString, TFile*> > file_map;
//   for( auto process : processes ){
//     for( auto chiral : chirals ){
//       Int_t processID = production.at({process,chiral}).first;
//       cout << process << " " << chiral << " " << processID << endl;
//       TString filename = inputDir + "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I" + processID + "." + process + "." + chiral + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root";
//       TFile *file = new TFile(filename,"READ");
//       if( !file->IsOpen() ) throw std::runtime_error("File not found");
//       file_map[process][chiral] = file;
//     }
//   }





// }




void total()
{
  TGaxis::SetMaxDigits(3);

  unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_cos_eff;  // [GenReco][LPFO][hist]

  for ( auto igenreco : gen_reco ){
    for ( auto i_lmode : PFO_mode ){
      for ( auto ih : heff_name ){
        Int_t tmp = 0;
        for( auto iqq : qq ){
          TString dir_name = "/efficiency/";
          TString hname = "h_" + iqq + "_" + igenreco + "_" + i_lmode + "_" + ih;
          TH1F *h = (TH1F*) file->Get(iqq + dir_name + hname);
          if(tmp==0) {
            TH1F *htmp = (TH1F*) h->Clone();
            htmp->SetName("h_" + igenreco + "_" + i_lmode + "_" + ih);
            h1_cos_eff[igenreco][i_lmode][ih] = htmp;
          }else{
            h1_cos_eff[igenreco][i_lmode][ih]->Add(h);
          }
          tmp++;
        }
      }
    }
  }

  int count = 0;
  TH1F * h_denom;
  TCanvas *c_eff_Pi = new TCanvas("c_eff_Pi", "c_eff_Pi", 1200,1200);
  // c_eff_Pi->Divide(heff_name.size()-1,1);
  c_eff_Pi->Divide(3,3);
  TCanvas *c_cos_Pi = new TCanvas("c_cos_Pi", "c_cos_Pi", 1200,1200);
  // c_cos_Pi->Divide(heff_name.size()-1,1);
  c_cos_Pi->Divide(3,3);
  
  for ( auto ih : heff_name ){

    TH1F * h_num = h1_cos_eff["reco"][LPFO_mode][ih];
    
    if (count) {
      
      TH1F *h_eff = plotEfficiency(h_num, h_denom);
      TString hname = "After " + ih + " selection (Cut " + count + ")";
      if(count==7) {
        // cutid = count+1;
        hname = TString::Format("After K^{-} PID selection (Cut %db)",count);
      }

      TString axis = TString::Format(";cos#theta;#epsilon_{%d}",count);
      
      TH1F *h_num_norm   = (TH1F*) h_num->Clone();
      TH1F *h_denom_norm = (TH1F*) h_denom->Clone();
      Normalize(h_num_norm);
      Normalize(h_denom_norm);

      c_eff_Pi->cd(count);
      StylePad(gPad,0,0,0.1,0.1);
      StyleHist(h_eff,kBlue);
      h_eff->SetTitle(hname+axis);
      h_eff->Draw("h");

      c_cos_Pi->cd(count);
      h_denom_norm->SetTitle(hname+axis);
      h_denom_norm->Draw("h");
      h_num_norm->Draw("hsame");

    }
    
    h_denom = h_num;
    count++;

  }

  TFile *file_eff_weight = new TFile("eff_weight.root","RECREATE");
  for ( auto igenreco : gen_reco ){
    for ( auto i_lmode : PFO_mode ){
      for ( auto ih : heff_name ){
        h1_cos_eff[igenreco][i_lmode][ih]->Write();
      }
    }
  }
  file_eff_weight->Close();



}

void partial()
{
  TGaxis::SetMaxDigits(3);

  unordered_map< TString, unordered_map< TString, unordered_map< TString, TH1F* > > > h1_cos_eff;  // [GenReco][LPFO][hist]

  try
  {
    if (!file->IsOpen()) return;
    TString dir_name = "/efficiency/";

    for ( auto igenreco : gen_reco ){
      for ( auto i_lmode : PFO_mode ){
        for ( auto ih : heff_name ){
          TString hname = "h_" + prod_mode + "_" + igenreco + "_" + i_lmode + "_" + ih;
          TH1F *h = (TH1F*) file->Get(prod_mode + dir_name + hname);
          h1_cos_eff[igenreco][i_lmode][ih] = h;
        }
      }
    }


    int count = 0;
    TH1F * h_denom;
    TCanvas *c_eff_Pi = new TCanvas("c_eff_Pi", "c_eff_Pi", 1500,400);
    c_eff_Pi->Divide(heff_name.size()-1,1);
    TCanvas *c_cos_Pi = new TCanvas("c_cos_Pi", "c_cos_Pi", 1500,400);
    c_cos_Pi->Divide(heff_name.size()-1,1);
    
    Int_t NbfCut = h1_cos_eff["reco"][LPFO_mode]["nocut"]->GetEntries();
    for ( auto ih : heff_name ){

      TH1F * h_num = h1_cos_eff["reco"][LPFO_mode][ih];

      Float_t total_ratio = (Float_t) h_num->GetEntries()/(Float_t) NbfCut;
      if(!count) cout << std::setprecision (3) << h_num->GetEntries() << endl;
      cout << std::setprecision (3) << total_ratio * 100 << "\\%" << endl;
      
      if (count) {
        
        TH1F *h_eff = plotEfficiency(h_num, h_denom);
        TString hname = "After " + ih + " selection";
        
        TH1F *h_num_norm   = (TH1F*) h_num->Clone();
        TH1F *h_denom_norm = (TH1F*) h_denom->Clone();
        Normalize(h_num_norm);
        Normalize(h_denom_norm);

        c_eff_Pi->cd(count);
        StylePad(gPad,0,0,0.17,0.1);
        StyleHist(h_eff,kBlue);
        h_eff->SetTitle(hname);
        h_eff->Draw("h");

        c_cos_Pi->cd(count);
        h_denom_norm->SetTitle(hname);
        h_denom_norm->Draw("h");
        h_num_norm->Draw("hsame");

      }
      
      h_denom = h_num;
      count++;

    }


  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}

void efficiency_cos()
{
  // partial();
  // total();

  calcEff();
  // eff_chiral();

}