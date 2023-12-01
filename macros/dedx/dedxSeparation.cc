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

const vector<TString> PFO_mode  = {"Pi","K"};
const vector<TString> PFO_type  = {"Pi","K", "p", "e", "mu"};
const int nPFO_type = PFO_type.size();
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

Float_t calculateSepPow(Float_t mean1, Float_t ey1, Float_t mean2, Float_t ey2)
{
  Float_t sigma1 = ey1;
  Float_t sigma2 = ey2;
  Float_t sigma12 = sqrt( (pow(sigma1,2) + pow(sigma2,2)) / 2 );
  return fabs(mean1 - mean2) / sigma12;
}

void registerSepPow(Int_t nbins, TH2F *h1, TH2F *h2, Float_t *spwr, Float_t *p)
{
  for (int i = 0; i < nbins; i++) {
    TH1D* proj1 = h1->ProjectionY("proj1", i, i+1);
    TH1D* proj2 = h2->ProjectionY("proj2", i, i+1);

    p[i] = h1->GetXaxis()->GetBinCenter(i+1);

    double mean1 = proj1->GetMean();
    double mean2 = proj2->GetMean();

    double ey1 = proj1->GetRMS();
    double ey2 = proj2->GetRMS();

    spwr[i] = calculateSepPow(mean1, ey1, mean2, ey2);

  }

}

void registerPurity(TH2F *hPi, TH2F *hK, TH2F *hp, Float_t *purity, Float_t *efficiency, Float_t *x)
{

  float momentum_min = 5;
  double ea=0.178;
  Int_t ea_bin = 0;
  Int_t dedx_Nbins = hPi->GetNbinsY();
  Int_t x_Nbins = hPi->GetNbinsX();
  
  // search threshold bin
  for(int ibin=1; ibin<=dedx_Nbins; ibin++){
    if(hPi->GetYaxis()->GetBinLowEdge(ibin)>ea){
      ea_bin=ibin;
      break;
    }
  }

  for (int i = 0; i < x_Nbins; i++) {

    double int_pion       = 0;
    double int_kaon       = 0;
    double int_proton     = 0;
    double int_pion_total = 0;

    TH1D * proj_pion = hPi->ProjectionY("proj_pion",i,i+1);
    int_pion_total = proj_pion->GetEntries();

    TH1D * proj_kaon   = hK->ProjectionY("proj_kaon",i,i+1);
    TH1D * proj_proton = hp->ProjectionY("proj_proton",i,i+1);

    for(int j1=0; j1<proj_pion->GetNbinsX(); j1++) {
      if(proj_pion->GetXaxis()->GetBinCenter(j1)> ea) int_pion+=proj_pion->GetBinContent(j1);
    }
    for(int j1=0; j1<proj_kaon->GetNbinsX(); j1++) {
      if(proj_kaon->GetXaxis()->GetBinCenter(j1)> ea) int_kaon+=proj_kaon->GetBinContent(j1);
    }
    for(int j1=0; j1<proj_proton->GetNbinsX(); j1++) {
      if(proj_proton->GetXaxis()->GetBinCenter(j1)> ea) int_proton+=proj_proton->GetBinContent(j1);
    }

    x[i]          = hPi->GetXaxis()->GetBinCenter(i+1);
    purity[i]     = (int_pion/(int_proton+int_pion+int_kaon))*100.;
    efficiency[i] = (int_pion/int_pion_total)*100.;

  }

}

void StyleGraph(TGraph *g, Color_t col, Int_t style)
{
  g->GetXaxis()->SetRangeUser(5,100);
  g->SetTitle(";p [GeV];Separation Power");
  g->SetLineColor(col);
  g->SetMarkerColor(col);
  g->SetLineWidth(3);
  g->SetMarkerStyle(style);
  g->SetMarkerSize(1.5);
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

      TString key_dedx_p   = h2_dEdx_eff_str.at(iprod_mode).at("reco").at(LPFO_mode).at(itype).at("nocut").at("dEdx_p");
      TString key_dedx_cos = h2_dEdx_eff_str.at(iprod_mode).at("reco").at(LPFO_mode).at(itype).at("nocut").at("dEdx_cos");
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


  TH2F *h_dedx_p_Pi = (TH2F*) h_dedx_p_vec.at(0)->Clone();
  TH2F *h_dedx_p_K  = (TH2F*) h_dedx_p_vec.at(1)->Clone();
  TH2F *h_dedx_p_p  = (TH2F*) h_dedx_p_vec.at(2)->Clone();

  TH2F *h_dedx_cos_Pi = (TH2F*) h_dedx_cos_vec.at(0)->Clone();
  TH2F *h_dedx_cos_K  = (TH2F*) h_dedx_cos_vec.at(1)->Clone();
  TH2F *h_dedx_cos_p  = (TH2F*) h_dedx_cos_vec.at(2)->Clone();



  Int_t NBinsP = h_dedx_p_Pi->GetNbinsX();
  Float_t p[NBinsP];

  // Projection plot
  TCanvas *c_type_dedx_p_proj = new TCanvas("c_type_dedx_p_proj", "c_type_dedx_p_proj", 800,800);
  TPad *pad_type_dedx_p_proj  = new TPad("pad_type_dedx_p_proj", "pad_type_dedx_p_proj",0,0,1,1);
  StylePad(pad_type_dedx_p_proj,0,0.15,0,0.17);
  c_type_dedx_p_proj->SetLogx();
  pad_type_dedx_p_proj->SetLogx();

  TLegend *leg_dedx_p_proj = new TLegend(0.62,0.67,0.80,0.83);
  leg_dedx_p_proj->SetTextSize(0.04);
  leg_dedx_p_proj->SetLineColor(0);
  leg_dedx_p_proj->SetFillStyle(0);
  leg_dedx_p_proj->SetMargin(0.8);  

  Float_t spwr_PiK[NBinsP];
  Float_t spwr_Pip[NBinsP];
  Float_t spwr_Kp[NBinsP];


  registerSepPow(NBinsP, h_dedx_p_Pi, h_dedx_p_K, spwr_PiK, p);
  registerSepPow(NBinsP, h_dedx_p_Pi, h_dedx_p_p, spwr_Pip, p);
  registerSepPow(NBinsP, h_dedx_p_K, h_dedx_p_p, spwr_Kp, p);

  TGraph *g_sep_pow_PiK = new TGraph(NBinsP, p, spwr_PiK);
  StyleGraph(g_sep_pow_PiK, kRed+1, 20);
  leg_dedx_p_proj->AddEntry(g_sep_pow_PiK, "#pi / K", "p");

  TGraph *g_sep_pow_Pip = new TGraph(NBinsP, p, spwr_Pip);
  StyleGraph(g_sep_pow_Pip, kBlue+1, 21);
  leg_dedx_p_proj->AddEntry(g_sep_pow_Pip, "#pi / p", "p");

  TGraph *g_sep_pow_Kp = new TGraph(NBinsP, p, spwr_Kp);
  StyleGraph(g_sep_pow_Kp, kGreen+1, 22);
  leg_dedx_p_proj->AddEntry(g_sep_pow_Kp, "K / p", "p");


  g_sep_pow_PiK->Draw("AP");
  g_sep_pow_Pip->Draw("Psame");
  g_sep_pow_Kp->Draw("Psame");

  leg_dedx_p_proj->Draw("same");


  // efficiency and purity

  TCanvas *c_pureff_dedx_p = new TCanvas("c_pureff_dedx_p", "c_pureff_dedx_p", 800,800);
  TPad  *pad_pureff_dedx_p = new TPad("pad_pureff_dedx_p", "pad_pureff_dedx_p",0,0,1,1);
  StylePad(pad_pureff_dedx_p,0,0.15,0,0.17);
  pad_pureff_dedx_p->SetGrid(1,1);

  TCanvas *c_pureff_dedx_cos = new TCanvas("c_pureff_dedx_cos", "c_pureff_dedx_cos", 800,800);
  TPad  *pad_pureff_dedx_cos = new TPad("pad_pureff_dedx_cos", "pad_pureff_dedx_cos",0,0,1,1);
  StylePad(pad_pureff_dedx_cos,0,0.15,0,0.17);
  pad_pureff_dedx_cos->SetGrid(1,1);

  TLegend *leg_dedx_p_purity_efficiency = new TLegend(0.57,0.18,0.75,0.34);
  leg_dedx_p_purity_efficiency->SetTextSize(0.035);
  leg_dedx_p_purity_efficiency->SetLineColor(0);
  leg_dedx_p_purity_efficiency->SetFillStyle(0);
  leg_dedx_p_purity_efficiency->SetMargin(0.8);

  TLegend *leg_dedx_cos_purity_efficiency = new TLegend(0.57,0.18,0.75,0.34);
  leg_dedx_cos_purity_efficiency->SetTextSize(0.035);
  leg_dedx_cos_purity_efficiency->SetLineColor(0);
  leg_dedx_cos_purity_efficiency->SetFillStyle(0);
  leg_dedx_cos_purity_efficiency->SetMargin(0.8);

  Int_t NBinsCos = h_dedx_cos_Pi->GetNbinsX();
  Float_t cos[NBinsCos];
  Float_t purity_p[NBinsP];
  Float_t purity_cos[NBinsCos];
  Float_t efficiency_p[NBinsP];
  Float_t efficiency_cos[NBinsCos];

  registerPurity(h_dedx_p_Pi, h_dedx_p_K, h_dedx_p_p, purity_p, efficiency_p, p);
  registerPurity(h_dedx_cos_Pi, h_dedx_cos_K, h_dedx_cos_p, purity_cos, efficiency_cos, cos);

  float ave_purity = 0;
  for(auto i : purity_cos){
    ave_purity += i;
  }
  ave_purity /= NBinsCos;
  cout << "ave_purity = " << ave_purity << endl;
  
  float ave_efficiency = 0;
  for(auto i : efficiency_cos){
    ave_efficiency += i;
  }
  ave_efficiency /= NBinsCos;
  cout << "ave_efficiency = " << ave_efficiency << endl;

  TGraph *g_purity_p = new TGraph(NBinsP, p, purity_p);
  StyleGraph(g_purity_p, kGreen+1, 20);
  leg_dedx_p_purity_efficiency->AddEntry(g_purity_p, "Purity", "p");
  TGraph *g_efficiency_p = new TGraph(NBinsP, p, efficiency_p);
  StyleGraph(g_efficiency_p, kBlue+1, 21);
  leg_dedx_p_purity_efficiency->AddEntry(g_efficiency_p, "Efficiency", "p");

  // plot purity and efficiency
  pad_pureff_dedx_p->cd();
  g_purity_p->SetTitle(";p [GeV];Ratio [%]");
  g_purity_p->GetXaxis()->SetRangeUser(5,100);
  g_purity_p->GetYaxis()->SetRangeUser(0,100);
  g_purity_p->Draw("AP");
  g_efficiency_p->Draw("Psame");
  leg_dedx_p_purity_efficiency->Draw("same");

  TGraph *g_purity_cos = new TGraph(NBinsCos, cos, purity_cos);
  StyleGraph(g_purity_cos, kGreen+1, 20);
  leg_dedx_cos_purity_efficiency->AddEntry(g_purity_cos, "Purity", "p");
  TGraph *g_efficiency_cos = new TGraph(NBinsCos, cos, efficiency_cos);
  StyleGraph(g_efficiency_cos, kBlue+1, 21);
  leg_dedx_cos_purity_efficiency->AddEntry(g_efficiency_cos, "Efficiency", "p");
  
  pad_pureff_dedx_cos->cd();
  g_purity_cos->GetXaxis()->SetRangeUser(-1,1);
  g_purity_cos->GetYaxis()->SetRangeUser(0,100);
  g_purity_cos->SetTitle(";cos#theta;Ratio [%]");
  g_purity_cos->Draw("AP");
  g_efficiency_cos->Draw("Psame");
  leg_dedx_cos_purity_efficiency->Draw("same");


}

void dedxSeparation()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  if (!file->IsOpen()) return;

  dedxPurity();

}