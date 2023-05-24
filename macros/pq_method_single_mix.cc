#include "Styles.cc"

vector<TH1F*> extract_ps_hists(TFile *file)
{
  TH1F *h_gen_N_K_cos  = (TH1F*) file->Get("h_gen_N_Pi_cos");
  TH1F *h_reco_N_K_cos = (TH1F*) file->Get("h_reco_N_Pi_cos");
  TH1F *h_N_K_corr_cos = (TH1F*) file->Get("h_N_Pi_corr_cos");

  TH1F *h_stable_cos = (TH1F*) h_N_K_corr_cos->Clone();
  TH1F *h_purity_cos = (TH1F*) h_N_K_corr_cos->Clone();
  h_stable_cos->Divide(h_gen_N_K_cos);
  h_purity_cos->Divide(h_reco_N_K_cos);

  TH1F *h_weight = (TH1F*) h_stable_cos->Clone();
  h_weight->Divide(h_purity_cos);

  vector<TH1F*> hists;
  hists.push_back(h_stable_cos);
  hists.push_back(h_purity_cos);
  hists.push_back(h_weight);

  return hists;

}

void single_mix(TFile *file)
{
  enum MixProcess {kUU,kDD,kUD};
  gStyle->SetOptStat(0);

  // gen uu/dd polar
  TH1F *h_gen_uu_qcos = (TH1F*) files[kUU]->Get("h_gen_q_qcos");
  TH1F *h_gen_dd_qcos = (TH1F*) files[kDD]->Get("h_gen_q_qcos");
  TH1F *h_gen_ud_qcos = (TH1F*) files[kUD]->Get("h_gen_q_qcos");

  StyleHist(h_gen_uu_qcos,kBlue+1);
  h_gen_uu_qcos->SetFillStyle(0);
  h_gen_uu_qcos->SetLineStyle(2);
  StyleHist(h_gen_dd_qcos,kGreen+2);
  h_gen_dd_qcos->SetFillStyle(0);
  h_gen_dd_qcos->SetLineStyle(2);
  StyleHist(h_gen_ud_qcos,kBlack);
  h_gen_ud_qcos->SetFillStyle(0);
  h_gen_ud_qcos->SetLineStyle(2);

  // reco us polar
  TH1F *h_reco_ud_Pi_scos  = (TH1F*) files[kUD]->Get("h_reco_Pi_scos");
  TH1F *h_reco_ud_Pi_qcos  = (TH1F*) files[kUD]->Get("h_reco_Pi_qcos");
  TH1F *h_cheat_ud_Pi_qcos = (TH1F*) files[kUD]->Get("h_cheat_Pi_qcos");

  // efficiency correction
  TH1F *h_reco_ud_Pi_scos_eff_corr  = Efficiency_Correction2(h_reco_ud_Pi_scos,"scos_corr",files[kUD]);
  TH1F *h_reco_ud_Pi_qcos_eff_corr  = Efficiency_Correction2(h_reco_ud_Pi_qcos,"qcos_corr",files[kUD]);
  TH1F *h_cheat_ud_Pi_qcos_eff_corr = Efficiency_Correction2(h_cheat_ud_Pi_qcos,"cheat_qcos_corr",files[kUD]);
  // TH1F *h_reco_ud_Pi_scos_eff_corr = (TH1F*) h_reco_ud_Pi_scos->Clone();
  // TH1F *h_reco_ud_Pi_qcos_eff_corr = (TH1F*) h_reco_ud_Pi_qcos->Clone();

  // used for pq correction
  TH1F *h_acc_PiPi_cos  = (TH1F*) files[kUD]->Get("pq/h_acc_PiPi_cos");
  TH1F *h_rej_PiPi_cos  = (TH1F*) files[kUD]->Get("pq/h_rej_PiPi_cos");

  TH1F *h_acc_PiPi_cos_eff_corr = Efficiency_Correction2(h_acc_PiPi_cos,"acc_corr",files[kUD]);
  TH1F *h_rej_PiPi_cos_eff_corr = Efficiency_Correction2(h_rej_PiPi_cos,"rej_corr",files[kUD]);
  // TH1F *h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
  // TH1F *h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();

  StyleHist(h_reco_ud_Pi_scos_eff_corr,kBlue);
  h_reco_ud_Pi_scos_eff_corr->SetFillStyle(0);
  StyleHist(h_reco_ud_Pi_qcos_eff_corr,kRed+2);
  StyleHist(h_cheat_ud_Pi_qcos_eff_corr,kBlue);
  StyleHist(h_cheat_ud_Pi_qcos,kBlue);

  StyleHist(h_acc_PiPi_cos_eff_corr,kRed+2);
  StyleHist(h_rej_PiPi_cos_eff_corr,kBlue+2);

  const Int_t nbins = h_reco_ud_Pi_scos_eff_corr->GetNbinsX();

  TH1F *p_KK = new TH1F("p_KK", "p_KK", 50,0,1);
  p_KK->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_PiPi_cos_eff_corr, h_rej_PiPi_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_KK->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_KK->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_Pi_pq_cos = CorrectHist(h_reco_ud_Pi_qcos_eff_corr, p_vec);
  // TH1F *h_reco_Pi_pq_cos = (TH1F*) h_reco_ud_Pi_qcos_eff_corr->Clone();
  StyleHist(h_reco_Pi_pq_cos,kBlack);


}



void pq_method_single_mix()
{
  TGaxis::SetMaxDigits(3);

  try
  {
    TString process[3] = {"uu","dd","ud"};
    TFile* files[3];
    for ( int i=0; i<3; i++ ){
      files[i] = new TFile(TString::Format("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.%s.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root",process[i].Data()),"READ");
      if ( !files[i]->IsOpen() ) throw 0;
    }

    main_pq_BGFit( files );
  }
  catch ( int error_code ) {
    switch ( error_code ){
      default:
        cerr << "<< Error >>" << endl;
        break;
    }
  }

}