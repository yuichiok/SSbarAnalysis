#include <iostream>
#include <vector>

using std::cout; using std::endl;
using std::vector;

vector<TString> modeLPFO = {"K","Pi"};
vector<TString> modePID  = {"K","Pi","p","e","mu"};

Float_t bins_cos_fine[] = {-1.0,-0.98,-0.96,-0.94,-0.92,-0.90,-0.88,-0.86,-0.84,-0.82,-0.80,-0.75,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,
                            0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.0};
// Int_t   nbins_cos = sizeof(bins_cos_fine) / sizeof(Float_t) - 1;
Int_t   nbins_cos = 100;

Float_t bins_cos_fine_half[] = {0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.0};
Int_t   nbins_cos_half = sizeof(bins_cos_fine_half) / sizeof(Float_t) - 1;

TString getProductionMode( TFile *file )
{
  TString file_name = file->GetName();
  TObjArray *tx = file_name.Tokenize(".");
  TString prod_mode;
  for ( auto i : *tx ) {
    TString channel = ((TObjString *)(i))->String();
    if(channel=="uu" || channel=="dd" || channel=="ud"){
      prod_mode = channel;
      break;
    }
  }
  return prod_mode;
}

vector<Float_t> GetP( TH1F * h_accepted, TH1F * h_rejected )
{
  const Int_t nbins = h_accepted->GetNbinsX();
  vector<Float_t> result_error;
  vector<Float_t> result;

  for (int i = 1; i < nbins / 2 + 1; i++)
  {
    std::vector<float> result_j;
    for (int i1 = -1; i1 < 2; i1 += 2)
    {
      for (int i2 = -1; i2 < 2; i2 += 2)
      {
        for (int i3 = -1; i3 < 2; i3 += 2)
        {
          float accepted = h_accepted->GetBinContent(nbins + 1 - i) + i1 * sqrt(h_accepted->GetBinContent(nbins + 1 - i));
          float rejected = h_rejected->GetBinContent(nbins + 1 - i) + i2 * sqrt(h_rejected->GetBinContent(nbins + 1 - i));
          accepted += h_accepted->GetBinContent(i) + i3 * sqrt(h_accepted->GetBinContent(i));

          float a = 1;
          float b = -1;
          float c = rejected / (2 * (accepted + rejected));
          float p = (0.5 / a) * (-b + sqrt(b * b - 4 * a * c));
          float p2 = (0.5 / a) * (-b - sqrt(b * b - 4 * a * c));
          if (p > 0.99)
            p = 0;
          if (p2 > 0.99)
            p2 = 0;
          if (p > 0 || p2 > 0)
            result_j.push_back(max(p, p2));
        }
      }
    }

    float average = 0;
    float n = 0;
    for (unsigned j = 0; j < result_j.size(); j++)
    {
      if (result_j.at(j) > 0)
      {
        average += result_j.at(j);
        n++;
      }
    }
    average /= n;

    if ( (average != 0) && (average == average) )
    {
      result.push_back(average);
      float std_dev = 0;
      for (unsigned j = 0; j < result_j.size(); j++)
      {
        if (result_j.at(j) > 0)
        {
          std_dev += pow(result_j.at(j) - average, 2);
        }
      }
      std_dev = sqrt(std_dev / (n - 1));
      result_error.push_back(std_dev);
    }
    else
    {
      result_error.push_back(0);
      result.push_back(0);
    }
  }

  for (unsigned i = 0; i < result_error.size(); i++)
  {
    if (result_error.at(i) > 0)
      result.push_back(result_error.at(i));
    else
      result.push_back(0);
  }

  return result;

}

TH1F * CorrectHist( TString prodMode, TH1F * h_reco, vector<Float_t> p_vec)
{
  const Int_t nbins = h_reco->GetNbinsX();

  TString corrected_name = "corrected_" + (TString) h_reco->GetName() + "_" + prodMode;
  TH1F *corrected = new TH1F(corrected_name, corrected_name, 100,-1,1);
  corrected->Sumw2();
  for (int i = 1; i < nbins / 2 + 1; i++)
  {
    float p = p_vec.at(i - 1);
    // float p = 0.82; // good AFB for ud mix
    // float p = 0.73;    // average p for u & d
    float q = 1 - p;
    float weight = (p * p + q * q) / (q * q * q * q - p * p * p * p);

    // calcualte average
    float av_i = 0;
    float av_41i = 0;
    float n = 0;
    for (int i1 = -1; i1 < 2; i1 += 1)
    {
      for (int i2 = -1; i2 < 2; i2 += 1)
      {
        float nm_reco_error = h_reco->GetBinContent(i) + i1 * h_reco->GetBinError(i);
        float np_reco_error = h_reco->GetBinContent(nbins + 1 - i) + i2 * h_reco->GetBinError(nbins + 1 - i);
        av_i += (np_reco_error * q * q - nm_reco_error * p * p) * weight;
        av_41i += -(np_reco_error * p * p - nm_reco_error * q * q) * weight;
        n++;
      }
    }
    av_i /= n;
    av_41i /= n;

    bool abn = 0;
    if(av_i<0 || av_41i<0){
      abn = true;
      av_i   = h_reco->GetBinContent(i);
      av_41i = h_reco->GetBinContent(nbins + 1 - i);
    }

    corrected->SetBinContent(i, av_i);
    corrected->SetBinContent(nbins + 1 - i, av_41i);

    // calculate error
    float error_i = 0;
    float error_41i = 0;
    n = 0;
    for (int i1 = -1; i1 < 2; i1 += 1)
    {
      for (int i2 = -1; i2 < 2; i2 += 1)
      {
        float nm_reco_error = h_reco->GetBinContent(i) + i1 * h_reco->GetBinError(i);
        float np_reco_error = h_reco->GetBinContent(nbins + 1 - i) + i2 * h_reco->GetBinError(nbins + 1 - i);
        error_i += pow((np_reco_error * q * q - nm_reco_error * p * p) * weight - av_i, 2);
        error_41i += pow(-(np_reco_error * p * p - nm_reco_error * q * q) * weight - av_41i, 2);
        n++;
      }
    }
    error_i = sqrt(error_i / (n - 1.));
    error_41i = sqrt(error_41i / (n - 1.));

    if(abn){
      error_i   = h_reco->GetBinError(i);
      error_41i = h_reco->GetBinError(nbins + 1 - i);
    }


    corrected->SetBinError(i, error_i);
    corrected->SetBinError(nbins + 1 - i, error_41i);
  }

  return corrected;

}

TH1F * efficiencyCorrection( TH1F * h, TString LPFO_mode, TFile * file, TString category )
{
  // TString dirName = category + "/efficiency/";
  // TH1F *h_reco_offset = (TH1F*) file->Get(dirName + "h_" + category + "_reco_" + LPFO_mode + "_momentum");
  // TH1F *h_reco_PID    = (TH1F*) file->Get(dirName + "h_" + category + "_reco_" + LPFO_mode + "_charge");
  TFile *file_eff_weight = new TFile("../efficiency/eff_weight_" + LPFO_mode + ".root","READ");
  TH1F *h_reco_offset = (TH1F*) file_eff_weight->Get("h_reco_" + LPFO_mode + "_btag");
  TH1F *h_reco_PID    = (TH1F*) file_eff_weight->Get("h_reco_" + LPFO_mode + "_charge");


  TH1F *h_eff = (TH1F*) h_reco_PID->Clone();
  h_eff->Sumw2();
  h_eff->Divide(h_reco_offset);

  TH1F *h_eff_corr = (TH1F*) h->Clone();
  h_eff_corr->Divide(h_eff);

  return h_eff_corr;
}

TH1F * resolutionCorrection( TH1F * h, TString LPFO_mode, TFile * file, TString category )
{
  TH1F *h_gen_N_Pi_cos  = (TH1F*) file->Get( category + "/resolution/h_" + category + "_" + LPFO_mode + "_gen_N_cos");
  TH1F *h_reco_N_Pi_cos = (TH1F*) file->Get( category + "/resolution/h_" + category + "_" + LPFO_mode + "_reco_N_cos");
  TH1F *h_N_Pi_corr_cos = (TH1F*) file->Get( category + "/resolution/h_" + category + "_" + LPFO_mode + "_N_corr_cos");
  Int_t genN = h_gen_N_Pi_cos->Integral();
  TH1F *htmp = (TH1F*) h->Clone();
  if(!genN) return htmp;

  TH1F *h_stable_cos = (TH1F*) h_N_Pi_corr_cos->Clone();
  TH1F *h_purity_cos = (TH1F*) h_N_Pi_corr_cos->Clone();
  h_stable_cos->Divide(h_gen_N_Pi_cos);
  h_purity_cos->Divide(h_reco_N_Pi_cos);

  if( h->GetNbinsX() != h_stable_cos->GetNbinsX() ) throw std::logic_error("Error");
  if( h->GetNbinsX() != h_purity_cos->GetNbinsX() ) throw std::logic_error("Error");

  TH1F *h_weight = (TH1F*) h_stable_cos->Clone();
  h_weight->Divide(h_purity_cos);
  TH1F *corrected = (TH1F*) h->Clone();
  corrected->Divide(h_weight);

  return corrected;

}

void Func2Hist( TH1F * h, double * par )
{
  for ( int ibin=1; ibin < nbins_cos+1; ibin++ ){

    Float_t x = h->GetXaxis()->GetBinCenter(ibin);
    Float_t binc_h = par[0]*(1+x*x)+par[1]*x;
    h->SetBinContent(ibin,binc_h);

  }
}

Float_t AFB_calculation( TF1 * f )
{
  float N_forward  = f->Integral(0,1);
  float N_backward = f->Integral(-1,0);

  float AFB = (N_forward - N_backward) / (N_forward + N_backward);

  return AFB;

}


Float_t AFB_N_calculation( Float_t N_forward, Float_t N_backward )
{
  Float_t AFB = (N_forward - N_backward) / (N_forward + N_backward);

  return AFB;

}

Float_t AFB_calculation_fit( TF1 * f )
{
  float S = f->GetParameter(0);
  float A = f->GetParameter(1);
  float AFB = (3*A) / (8*S);

  return AFB;

}

Float_t AFB_error( Float_t AFB, Int_t N )
{
  Float_t AFB_error = sqrt((1-AFB*AFB)/(Float_t)N);

  return AFB_error;

}