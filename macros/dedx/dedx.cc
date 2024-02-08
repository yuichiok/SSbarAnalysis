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

void get_dEdx()
{
  for (auto iprod_mode : prod_modes ){
    TString dir_dedx_name = iprod_mode + "/dEdx/";
    for ( auto imode : PFO_mode ){
      for ( auto itype : PFO_type ){
        for ( auto icut : hdEdx_name ){
          TString hname = "h2_" + iprod_mode + "_" + imode + "_" + itype + "_" + icut;
          h2_dEdx_str[iprod_mode][imode][itype][icut] = dir_dedx_name + hname;
          // TH2F *h = (TH2F*) file->Get(dir_dedx_name + hname);
          // h2_dEdx[iprod_mode][imode][itype][icut] = h;
        }
      }
    }
  }
}

void dedxP()
{
  get_dEdx_eff();

  TCanvas *c_dedx_p_type = new TCanvas("c_dedx_p_type", "c_dedx_p_type", 800,800);
  TPad *pad_dedx_p_type  = new TPad("pad_dedx_p_type", "pad_dedx_p_type",0,0,1,1);
  StylePad(pad_dedx_p_type,0,0.15,0,0.17);
  TLegend *leg = new TLegend(0.5,0.76,0.75,0.85);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.8);

  int count_draw = 0;
  for ( auto itype : PFO_type ){

    TH2F *h_dedx_p;
    Int_t countAdd = 0;

    for ( auto iprod_mode : prod_modes ){

      // TH2F *h_dedx_p_tmp = h2_dEdx_eff.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at("nocut").at("dEdx_p");
      TString key = h2_dEdx_eff_str.at(iprod_mode).at("reco").at(LPFO_mode).at(itype).at("nocut").at("dEdx_p");
      TH2F *h_dedx_p_tmp = (TH2F*) file->Get(key);
      h_dedx_p_tmp->SetFillColor(type_color_map[itype]);
      h_dedx_p_tmp->SetLineColor(type_color_map[itype]);
      h_dedx_p_tmp->GetXaxis()->SetRangeUser(15,100);

      if(countAdd){
        h_dedx_p->Add(h_dedx_p_tmp);
      }else{
        h_dedx_p = (TH2F*) h_dedx_p_tmp->Clone();
      }

      countAdd++;
    }

    if(count_draw) h_dedx_p->Draw("box same");
    else{
      h_dedx_p->SetTitle(";p [GeV];dE/dx #times 10^{-6} GeV/mm");
      h_dedx_p->Draw("box");
    }
    leg->AddEntry(h_dedx_p,itype,"f");
    count_draw++;

  }

  leg->Draw();


}

void dedxDistCosProjType()
{
  get_dEdx_eff();

  TCanvas *c_type_dedx_cos_type = new TCanvas("c_type_dedx_cos_type", "c_type_dedx_cos_type", 1500,1000);
  c_type_dedx_cos_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_p_type = new TCanvas("c_type_dedx_p_type", "c_type_dedx_p_type", 1500,1000);
  c_type_dedx_p_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_dist_cos_type = new TCanvas("c_type_dedx_dist_cos_type", "c_type_dedx_dist_cos_type", 1500,1000);
  c_type_dedx_dist_cos_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_dist_cos_proj_type = new TCanvas("c_type_dedx_dist_cos_proj_type", "c_type_dedx_dist_cos_proj_type", 1500,1000);
  c_type_dedx_dist_cos_proj_type->Divide(cut_name.size(),PFO_type.size());

  TCanvas *c_type_dedx_dist_cos_proj_efficiency_type = new TCanvas("c_type_dedx_dist_cos_proj_efficiency_type", "c_type_dedx_dist_cos_proj_efficiency_type", 1500,1000);
  c_type_dedx_dist_cos_proj_efficiency_type->Divide(cut_name.size(),PFO_type.size());

  TH1F *PiID_eff;

  TPad *p1;

  int count = 0;
  for ( auto itype : PFO_type ){
    TH1F *h_denom = nullptr;
    for ( auto iname : cut_name ){
      
      count++;

      TString hTitle = itype + " | " + iname;
      // TH2F *h_dedx_cos           = h2_dEdx_eff.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_cos");
      // TH2F *h_dedx_p             = h2_dEdx_eff.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_p");
      // TH2F *h_dedx_dist_cos      = h2_dEdx_eff.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_dist_cos");
      TString key_dedx_cos = h2_dEdx_eff_str.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_cos");
      TString key_dedx_p   = h2_dEdx_eff_str.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_p");
      TString key_dedx_dist_cos = h2_dEdx_eff_str.at(prod_mode).at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_dist_cos");

      TH2F *h_dedx_cos           = (TH2F*) file->Get(key_dedx_cos);
      TH2F *h_dedx_p             = (TH2F*) file->Get(key_dedx_p);
      TH2F *h_dedx_dist_cos      = (TH2F*) file->Get(key_dedx_dist_cos);
      TH1F *h_dedx_dist_cos_proj = (TH1F*) h_dedx_dist_cos->ProjectionX();

      h_dedx_dist_cos->SetTitle(hTitle);
      h_dedx_dist_cos_proj->SetTitle(hTitle);

      c_type_dedx_cos_type->cd(count);
      h_dedx_cos->Draw("colz");
      // h_dedx_cos->Draw("cont3 same");

      c_type_dedx_p_type->cd(count);
      p1 = (TPad*)(c_type_dedx_p_type->cd(count));
      p1->SetLogx();
      h_dedx_p->GetXaxis()->SetRangeUser(15,100);
      h_dedx_p->Draw("colz");

      c_type_dedx_dist_cos_type->cd(count);
      h_dedx_dist_cos->Draw("colz");
      // h_dedx_dist_cos->Draw("cont3 same");
      
      c_type_dedx_dist_cos_proj_type->cd(count);
      StyleHist(h_dedx_dist_cos_proj,kBlue);
      h_dedx_dist_cos_proj->Draw("h");

      TH1F *h_efficiency;
      if(h_denom){
        h_efficiency = plotEfficiency(h_dedx_dist_cos_proj, h_denom);
        h_denom = (TH1F*) h_dedx_dist_cos_proj->Clone();
        c_type_dedx_dist_cos_proj_efficiency_type->cd(count);
        StyleHist(h_efficiency,kBlue);
        h_efficiency->Draw("h");
        if(itype == "Pi" && iname == "PID"){
          PiID_eff = (TH1F*) h_efficiency->Clone();
        }
      }else{
        h_denom = (TH1F*) h_dedx_dist_cos_proj->Clone();        
      }

    }

  }

}

void dedxDistCosProj()
{
  get_dEdx();

  TCanvas *c_dedx_dist_cos = new TCanvas("c_dedx_dist_cos", "c_dedx_dist_cos", 500,500);
  TCanvas *c_type_dedx_dist_cos = new TCanvas("c_type_dedx_dist_cos", "c_type_dedx_dist_cos", 1500,400);
  c_type_dedx_dist_cos->Divide(5,1);

  int count = 0;
  TH1F * h_dedx_dist_cos_proj_sum = new TH1F("h_dedx_dist_cos_proj_sum","h_dedx_dist_cos_proj_sum",100,-1,1);
  StyleHist(h_dedx_dist_cos_proj_sum,kBlue);
  for ( auto itype : PFO_type ){
    // TH2F *h_dedx_dist_cos = h2_dEdx.at(prod_mode).at(LPFO_mode).at(itype).at("dEdx_dist_cos");
    TString key_dedx_dist_cos = h2_dEdx_str.at(prod_mode).at(LPFO_mode).at(itype).at("dEdx_dist_cos");
    TH2F *h_dedx_dist_cos = (TH2F*) file->Get(key_dedx_dist_cos);
    TH1F *h_dedx_dist_cos_proj = (TH1F*) h_dedx_dist_cos->ProjectionX();
    StyleHist(h_dedx_dist_cos_proj,kBlue);

    h_dedx_dist_cos_proj_sum->Add(h_dedx_dist_cos_proj);

    c_type_dedx_dist_cos->cd(++count);
    h_dedx_dist_cos_proj->Draw("h");
  }

  c_dedx_dist_cos->cd();
  h_dedx_dist_cos_proj_sum->Draw("h");

}

void dedxCos()
{
  get_dEdx();

  TCanvas *c_type_dedx_cos = new TCanvas("c_type_dedx_cos", "c_type_dedx_cos", 1500,400);
  c_type_dedx_cos->Divide(5,1);

  int count = 0;
  for ( auto itype : PFO_type ){

    count++;

    // TH2F *h_dedx_cos = h2_dEdx.at(prod_mode).at(LPFO_mode).at(itype).at("dEdx_cos");
    TString key_dedx_cos = h2_dEdx_str.at(prod_mode).at(LPFO_mode).at(itype).at("dEdx_cos");
    TH2F *h_dedx_cos = (TH2F*) file->Get(key_dedx_cos);

    c_type_dedx_cos->cd(count);
    h_dedx_cos->Draw("colz");

  }

}

void dedxOffsetProjection()
{
  get_dEdx_eff();

  // legend
  TLegend *leg_dedx_p = new TLegend(0.62,0.67,0.80,0.83);
  leg_dedx_p->SetLineColor(0);
  leg_dedx_p->SetFillStyle(0);
  leg_dedx_p->SetMargin(0.8);

  TLegend *leg_dedx_cos = new TLegend(0.62,0.67,0.80,0.83);
  leg_dedx_cos->SetLineColor(0);
  leg_dedx_cos->SetFillStyle(0);
  leg_dedx_cos->SetMargin(0.8);

  Int_t count = 0;

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


  TCanvas *c_type_dedx_p = new TCanvas("c_type_dedx_p", "c_type_dedx_p", 800,800);
  TPad *pad_type_dedx_p  = new TPad("pad_type_dedx_p", "pad_type_dedx_p",0,0,1,1);
  StylePad(pad_type_dedx_p,0,0.15,0,0.17);
  pad_type_dedx_p->SetLogx();

  TCanvas *c_type_dedx_cos = new TCanvas("c_type_dedx_cos", "c_type_dedx_cos", 800,800);
  TPad *pad_type_dedx_cos  = new TPad("pad_type_dedx_cos", "pad_type_dedx_cos",0,0,1,1);
  StylePad(pad_type_dedx_cos,0,0.15,0,0.17);

  int count_draw = 0;
  for (int i=0; i < h_dedx_p_vec.size(); i++){
    TH2F *h_dedx_p   = h_dedx_p_vec.at(i);
    TH2F *h_dedx_cos = h_dedx_cos_vec.at(i);

    h_dedx_p->GetXaxis()->SetRangeUser(15,100);

    if(count_draw){
      pad_type_dedx_p->cd();
      h_dedx_p->Draw("box same");

      pad_type_dedx_cos->cd();
      h_dedx_cos->Draw("box same");
    }else{
      pad_type_dedx_p->cd();
      h_dedx_p->Draw("box");
      
      pad_type_dedx_cos->cd();
      h_dedx_cos->Draw("box");
    }
    count_draw++;

  }

  pad_type_dedx_p->cd();
  leg_dedx_p->Draw();

  pad_type_dedx_cos->cd();
  leg_dedx_cos->Draw();


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
    h_dedx_p_proj->SetTitle(";dE/dx #times 10^{-6} [GeV/mm];Entries");

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

void dedxOffsetMeanSigma()
{
  get_dEdx_eff();

  unordered_map<TString, TGraphErrors*> type_tge; // [type]

  for ( const auto ipair : type_color ){

    TString itype  = ipair.first;
    Color_t icolor = ipair.second;

    TH2F *h_dedx_cos;
    Int_t countAdd = 0;
    for ( const auto iprod_mode : prod_modes ){
      TString key_dedx_cos = h2_dEdx_eff_str.at(iprod_mode).at("reco").at(LPFO_mode).at(itype).at("nocut").at("dEdx_cos");
      TH2F *h_dedx_cos_tmp = (TH2F*) file->Get(key_dedx_cos);
      if(countAdd){
        h_dedx_cos->Add(h_dedx_cos_tmp);
      }else{
        h_dedx_cos = (TH2F*) h_dedx_cos_tmp->Clone();
      }
      countAdd++;
    }


    Int_t nbins = h_dedx_cos->GetNbinsX();
    Float_t tge_cos[nbins], tge_err_cos[nbins];
    Float_t tge_dedx[nbins], tge_err_dedx[nbins];
    for( int ibin=1; ibin <= nbins; ibin++ ){
      TH1F *h_dedx_cos_proj = (TH1F*) h_dedx_cos->ProjectionY("h_dedx_cos_proj",ibin,ibin);
      TF1 *f_dedx_cos = new TF1("f_dedx_cos","gaus",0.1,0.3);
      h_dedx_cos_proj->Fit(f_dedx_cos,"MNRSQ");
      Float_t cos  = h_dedx_cos->GetXaxis()->GetBinCenter(ibin);
      Float_t binw = h_dedx_cos->GetXaxis()->GetBinWidth(ibin);
      Float_t mean_dedx  = f_dedx_cos->GetParameter(1);
      Float_t sigma_dedx = f_dedx_cos->GetParameter(2);

      tge_cos[ibin-1]  = cos;
      tge_dedx[ibin-1] = mean_dedx;
      tge_err_cos[ibin-1]  = binw/2.;
      tge_err_dedx[ibin-1] = sigma_dedx;
    }

    if(itype == "Pi" || itype == "K" || itype == "p"){
      type_tge[itype] = new TGraphErrors(nbins,tge_cos,tge_dedx,tge_err_cos,tge_err_dedx);
      type_tge[itype]->SetName(itype);
      type_tge[itype]->SetTitle(itype+" dE/dx Projection;cos#theta;dE/dx");
    }

  }

  TFile *fout = new TFile("dedxOffsetMeanSigma.root","RECREATE");
  TCanvas *c_dedx_cos_proj = new TCanvas("c_dedx_cos_proj", "c_dedx_cos_proj", 800,800);
  TPad *pad_dedx_cos_proj  = new TPad("pad_dedx_cos_proj", "pad_dedx_cos_proj",0,0,1,1);
  StylePad(pad_dedx_cos_proj,0,0.15,0,0.17);

  // legend
  TLegend *leg_dedx_cos_proj = new TLegend(0.59,0.20,0.77,0.36);
  leg_dedx_cos_proj->SetLineColor(0);
  leg_dedx_cos_proj->SetMargin(0.8);

  Int_t counter = 0;
  for ( const auto ipair : type_color ){

    TString itype  = ipair.first;
    Color_t icolor = ipair.second;
    if(!(itype == "Pi" || itype == "K" || itype == "p")) continue;

    TGraphErrors *tge = type_tge[itype];

    pad_dedx_cos_proj->cd();
    tge->SetMarkerColor(icolor);
    tge->SetLineColor(icolor);
    tge->SetMarkerStyle(20);
    tge->SetMarkerSize(0.5);
    leg_dedx_cos_proj->AddEntry(tge,itype,"lep");
    if(counter) tge->Draw("P same");
    else{
      tge->SetTitle(";cos#theta;dE/dx #times 10^{-6} GeV/mm");
      tge->GetYaxis()->SetRangeUser(0.1,0.2);
      tge->Draw("AP");
    }

    fout->cd();
    tge->Write();

    counter++;
  }


  leg_dedx_cos_proj->Draw();

  fout->Close();

}

void dedx()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  if (!file->IsOpen()) return;

  // getHistograms();

  // dedxP();
  // dedxDistCosProj();
  // dedxDistCosProjType();
  // dedxCos();
  dedxOffsetProjection();
  // dedxOffsetMeanSigma();


}