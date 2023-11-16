#include <iostream>
#include <vector>
#include <map>

#include "../../SSbarLibrary/include/MapTString.hh"
#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector;
using std::unordered_map;
using std::pair;

TString prod_mode = "ss";
TString chiral    = "eL.pR";
TString LPFO_mode = "Pi";

// TString prod_modes[5] = {"dd","uu","ss","cc","bb"};
TString prod_modes[3] = {"dd","uu","ss"};

// TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h." + chiral + "." + prod_mode + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root","READ");
TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root","READ");

const vector<TString> PFO_mode  = {"K","Pi"};
const vector<TString> PFO_type  = {"K","Pi", "p", "e", "mu"};
vector< pair< TString, Color_t > > type_color = {{"Pi",kBlue},{"K",kRed},{"p",kGreen+1},{"e",kMagenta},{"mu",kCyan}};
unordered_map<TString, Color_t> type_color_map = {{"Pi",kBlue},{"K",kRed},{"p",kGreen+1},{"e",kMagenta},{"mu",kCyan}};

TH1F* plotEfficiency(TH1F *h_num, TH1F *h_denom)
{
  gStyle->SetOptStat(0);

  TH1F *h_eff = (TH1F*) h_num->Clone();
  h_eff->Divide(h_denom);
  h_eff->GetYaxis()->SetRangeUser(0,1);

  return h_eff;
}

// unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > > > h2_dEdx_eff;  // [prod_mode][GenReco][LPFO][TruthID][cut][hist]
unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TString > > > > > > h2_dEdx_eff_str;  // [prod_mode][GenReco][LPFO][TruthID][cut][hist]
vector<TString> gen_reco  = {"gen","reco"};
vector<TString> cut_name = {"nocut","momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
vector<TString> heff_dedx_name = {"dEdx_p","dEdx_cos","dEdx_error_cos","dEdx_dist_cos"};

// unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > h2_dEdx; // [prod_mode][LPFO][TruthID][hist]
unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TString > > > > h2_dEdx_str; // [prod_mode][LPFO][TruthID][hist]
vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};

void get_dEdx_eff()
{
  for (auto iprod_mode : prod_modes ){
    TString dir_eff_name = iprod_mode + "/efficiency/dEdx/";
    for ( auto igen_reco : gen_reco ){
      for ( auto imode : PFO_mode ){
        for ( auto itype : PFO_type ){
          for ( auto icut : cut_name ){
            for ( auto ihist : heff_dedx_name ){
              TString hname = "h2_" + iprod_mode + "_" + igen_reco + "_" + imode + "_" + itype + "_" + icut + "_" + ihist;
              h2_dEdx_eff_str[iprod_mode][igen_reco][imode][itype][icut][ihist] = dir_eff_name + hname;
              // TH2F *h = (TH2F*) file->Get(dir_eff_name + hname);
              // h2_dEdx_eff[iprod_mode][igen_reco][imode][itype][icut][ihist] = h;
            }
          }
        }
      }
    }
  }
}

void dedxPurity()
{
  get_dEdx_eff();

  vector<TH2F*> h_dedx_p_vec;
  vector<TH2F*> h_dedx_cos_vec;

  for ( const auto ipair : type_color ){

    TString itype  = ipair.first;
    Color_t icolor = ipair.second;

    TH2F *h_dedx_p;
    TH2F *h_dedx_cos;
    Int_t countAdd = 0;
    for ( const auto iprod_mode : prod_modes ){

      TString key_dedx_p   = h2_dEdx_eff_str.at(iprod_mode).at("reco").at(LPFO_mode).at(itype).at("offset").at("dEdx_p");
      TString key_dedx_cos = h2_dEdx_eff_str.at(iprod_mode).at("reco").at(LPFO_mode).at(itype).at("offset").at("dEdx_cos");
      TH2F *h_dedx_p_tmp   = (TH2F*) file->Get(key_dedx_p);
      TH2F *h_dedx_cos_tmp = (TH2F*) file->Get(key_dedx_cos);
      if(countAdd){
        h_dedx_p->Add(h_dedx_p_tmp);
        h_dedx_cos->Add(h_dedx_cos_tmp);
      }else{
        h_dedx_p   = (TH2F*) h_dedx_p_tmp->Clone();
        h_dedx_cos = (TH2F*) h_dedx_cos_tmp->Clone();
      }
      countAdd++;
    }

    leg_dedx_p->AddEntry(h_dedx_p,itype,"f");
    leg_dedx_cos->AddEntry(h_dedx_cos,itype,"f");

    h_dedx_p->SetTitle(";p [GeV];dE/dx #times 10^{-6} [GeV/mm]");
    h_dedx_p->SetLineColor(icolor);
    h_dedx_p->SetFillColor(icolor);
    h_dedx_p->SetMarkerColor(icolor);
    h_dedx_p->SetMarkerStyle(20);
    h_dedx_p->SetMarkerSize(0.5);

    h_dedx_cos->SetTitle(";cos#theta;dE/dx #times 10^{-6} [GeV/mm]");
    h_dedx_cos->SetLineColor(icolor);
    h_dedx_cos->SetFillColor(icolor);
    h_dedx_cos->SetMarkerColor(icolor);
    h_dedx_cos->SetMarkerStyle(20);
    h_dedx_cos->SetMarkerSize(0.5);

    h_dedx_p_vec.push_back(h_dedx_p);
    h_dedx_cos_vec.push_back(h_dedx_cos);

  }

  // Projection plot
  TCanvas *c_type_dedx_p_proj = new TCanvas("c_type_dedx_p_proj", "c_type_dedx_p_proj", 800,800);
  TPad *pad_type_dedx_p_proj  = new TPad("pad_type_dedx_p_proj", "pad_type_dedx_p_proj",0,0,1,1);
  StylePad(pad_type_dedx_p_proj,0,0.15,0,0.17);

  TLegend *leg_dedx_p_proj = new TLegend(0.62,0.67,0.80,0.83);
  leg_dedx_p_proj->SetLineColor(0);
  leg_dedx_p_proj->SetFillStyle(0);
  leg_dedx_p_proj->SetMargin(0.8);

  count_draw=0;
  for (int i=0; i < h_dedx_p_vec.size(); i++){
    TH1F *h_dedx_p_proj   = (TH1F*) h_dedx_p_vec.at(i)->ProjectionY();

    StyleHist(h_dedx_p_proj,type_color.at(i).second);

    leg_dedx_p_proj->AddEntry(h_dedx_p_proj,PFO_type.at(i),"l");
    pad_type_dedx_p_proj->cd();

    if(count_draw){
      h_dedx_p_proj->Draw("box same");
    }else{
      h_dedx_p_proj->Draw("box");
    }
    count_draw++;

  }

  pad_type_dedx_p_proj->cd();
  leg_dedx_p_proj->Draw();

}

void dedxSeparation()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  if (!file->IsOpen()) return;

  dedxPurity();

}