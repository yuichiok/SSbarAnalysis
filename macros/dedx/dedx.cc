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

TString prod_mode = "dd";
TString chiral    = "eR.pL";
TString LPFO_mode = "Pi";

TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h." + chiral + "." + prod_mode + ".KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root","READ");

const vector<TString> PFO_mode  = {"K","Pi"};
const vector<TString> PFO_type  = {"K","Pi", "p", "e", "mu"};
vector< pair< TString, Color_t > > type_color = {{"Pi",kBlue},{"K",kRed},{"p",kGreen},{"e",kMagenta},{"mu",kCyan}};

TH1F* plotEfficiency(TH1F *h_num, TH1F *h_denom)
{
  gStyle->SetOptStat(0);

  TH1F *h_eff = (TH1F*) h_num->Clone();
  h_eff->Divide(h_denom);
  h_eff->GetYaxis()->SetRangeUser(0,1);

  return h_eff;
}

unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > > > h2_dEdx_eff;  // [GenReco][LPFO][TruthID][cut][hist]
vector<TString> gen_reco  = {"gen","reco"};
vector<TString> cut_name = {"nocut","momentum", "tpc_hits", "offset", "PID", "SPFO", "charge"};
vector<TString> heff_dedx_name = {"dEdx_p","dEdx_cos","dEdx_error_cos","dEdx_dist_cos"};

unordered_map< TString, unordered_map< TString, unordered_map< TString, TH2F* > > > h2_dEdx; // [LPFO][TruthID][hist]
vector<TString> hdEdx_name = {"dEdx_p","dEdx_cos","dEdx_dist_cos"};

void getHistograms()
{
  TString dir_eff_name = "efficiency/dEdx/";

  for ( auto igen_reco : gen_reco ){
    for ( auto imode : PFO_mode ){
      for ( auto itype : PFO_type ){
        for ( auto icut : cut_name ){
          for ( auto ihist : heff_dedx_name ){
            TString hname = "h2_" + igen_reco + "_" + imode + "_" + itype + "_" + icut + "_" + ihist;
            TH2F *h = (TH2F*) file->Get(dir_eff_name + hname);
            h2_dEdx_eff[igen_reco][imode][itype][icut][ihist] = h;
          }
        }
      }
    }
  }

  TString dir_dedx_name = "dEdx/";

  for ( auto imode : PFO_mode ){
    for ( auto itype : PFO_type ){
      for ( auto icut : hdEdx_name ){
        TString hname = "h2_" + imode + "_" + itype + "_" + icut;
        TH2F *h = (TH2F*) file->Get(dir_dedx_name + hname);
        h2_dEdx[imode][itype][icut] = h;
      }
    }
  }
}

void dedxDistCosProjType()
{
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
      TH2F *h_dedx_cos           = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_cos");
      TH2F *h_dedx_p             = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_p");
      TH2F *h_dedx_dist_cos      = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at(iname).at("dEdx_dist_cos");
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
  TCanvas *c_dedx_dist_cos = new TCanvas("c_dedx_dist_cos", "c_dedx_dist_cos", 500,500);
  TCanvas *c_type_dedx_dist_cos = new TCanvas("c_type_dedx_dist_cos", "c_type_dedx_dist_cos", 1500,400);
  c_type_dedx_dist_cos->Divide(5,1);

  int count = 0;
  TH1F * h_dedx_dist_cos_proj_sum = new TH1F("h_dedx_dist_cos_proj_sum","h_dedx_dist_cos_proj_sum",100,-1,1);
  StyleHist(h_dedx_dist_cos_proj_sum,kBlue);
  for ( auto itype : PFO_type ){
    TH2F *h_dedx_dist_cos = h2_dEdx.at(LPFO_mode).at(itype).at("dEdx_dist_cos");
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
  TCanvas *c_type_dedx_cos = new TCanvas("c_type_dedx_cos", "c_type_dedx_cos", 1500,400);
  c_type_dedx_cos->Divide(5,1);

  int count = 0;
  for ( auto itype : PFO_type ){

    count++;

    TH2F *h_dedx_cos = h2_dEdx.at(LPFO_mode).at(itype).at("dEdx_cos");

    c_type_dedx_cos->cd(count);
    h_dedx_cos->Draw("colz");

  }

}

void dedxOffsetProjection()
{
  TCanvas *c_type_dedx_p = new TCanvas("c_type_dedx_p", "c_type_dedx_p", 800,800);
  TPad *pad_type_dedx_p  = new TPad("pad_type_dedx_p", "pad_type_dedx_p",0,0,1,1);
  StylePad(pad_type_dedx_p,0,0.15,0,0.17);

  // legend
  TLegend *leg = new TLegend(0.62,0.67,0.80,0.83);
  leg->SetLineColor(0);
  leg->SetMargin(0.8);

  Int_t count = 0;
  for ( const auto ipair : type_color ){

    TString itype  = ipair.first;
    Color_t icolor = ipair.second;

    count++;
    TH2F *h_dedx_p   = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at("offset").at("dEdx_p");
    TH2F *h_dedx_cos = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at("offset").at("dEdx_cos");
    TH1F *h_dedx_p_proj   = (TH1F*) h_dedx_p->ProjectionY();
    StyleHist(h_dedx_p_proj,icolor);

    pad_type_dedx_p->cd();
    if(count){
      h_dedx_p_proj->Draw("hsame");
    }else{
      h_dedx_p_proj->SetTitle(";dE/dx;Entries");
      h_dedx_p_proj->Draw("h");
    }
    leg->AddEntry(h_dedx_p_proj,itype,"l");

  }

  pad_type_dedx_p->cd();
  leg->Draw();

}

void dedxOffsetMeanSigma()
{
  unordered_map<TString, TGraphErrors*> type_tge; // [type]

  for ( const auto ipair : type_color ){

    TString itype  = ipair.first;
    Color_t icolor = ipair.second;

    TH2F *h_dedx_cos = h2_dEdx_eff.at("reco").at(LPFO_mode).at(itype).at("offset").at("dEdx_cos");

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
  for ( const auto &[itype, tge] : type_tge ){
    tge->Write();
  }
  fout->Close();

}

void dedx()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  if (!file->IsOpen()) return;

  getHistograms();

  // dedxDistCosProj();
  dedxDistCosProjType();
  // dedxCos();
  // dedxOffsetProjection();
  // dedxOffsetMeanSigma();


}