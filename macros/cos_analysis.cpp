#include <iostream>
#include <vector>

using std::cout; using std::endl;
using std::vector;

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

TH1F * CorrectHist( TH1F * h_reco, vector<Float_t> p_vec)
{
  const Int_t nbins = h_reco->GetNbinsX();

  TH1F *corrected = new TH1F("corrected", "corrected", 100,-1,1);
  corrected->Sumw2();
  for (int i = 1; i < nbins / 2 + 1; i++)
  {
    float p = p_vec.at(i - 1);
    float q = 1 - p;
    if(p==0.5){
        corrected->SetBinContent(i, h_reco->GetBinContent(i));
        corrected->SetBinContent(nbins + 1 - i, h_reco->GetBinContent(nbins + 1 - i));

        corrected->SetBinError(i, h_reco->GetBinError(i));
        corrected->SetBinError(nbins + 1 - i, h_reco->GetBinError(nbins + 1 - i));
    }
    else
    {
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
        corrected->SetBinError(i, error_i);
        corrected->SetBinError(nbins + 1 - i, error_41i);
    }
  }

  return corrected;

}

void cos_analysis() {
    TFile* file_ss = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.PFOp5.yevhenii.all.root");
    TH1F* h_full_ss = (TH1F*)file_ss->Get("/cos_theta/cos_theta");
    TH1F* h_acc_ss = (TH1F*)file_ss->Get("/cos_theta/acc_cos_theta");
    TH1F* h_rej_ss = (TH1F*)file_ss->Get("/cos_theta/rej_cos_theta");
                
    vector<Float_t> p_ss = GetP(h_acc_ss, h_rej_ss);
    TH1F* h_ss = CorrectHist(h_acc_ss, p_ss);

    TFile* file_cc = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.cc.PFOp5.yevhenii.all.root");
    TH1F* h_full_cc = (TH1F*)file_cc->Get("/cos_theta/cos_theta");
    TH1F* h_acc_cc = (TH1F*)file_cc->Get("/cos_theta/acc_cos_theta");
    TH1F* h_rej_cc = (TH1F*)file_cc->Get("/cos_theta/rej_cos_theta");
                
    vector<Float_t> p_cc = GetP(h_acc_cc, h_rej_cc);
    TH1F* h_cc = CorrectHist(h_acc_cc, p_cc);

    // h_acc_ss->Scale(1./h_acc_ss->Integral(50,150));
    // h_ss->Scale(1./h_ss->Integral(50,150));

    // h_acc_cc->Scale(1./h_acc_cc->Integral(50,150));
    // h_cc->Scale(1./h_cc->Integral(50,150));


    TCanvas* canvas = new TCanvas("c","c",1600,800);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);
    canvas->SetFrameLineWidth(2);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderSize(0);
    canvas->SetLeftMargin(0.05);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);

    canvas->Divide(2,1);
    canvas->cd(1);
    TVirtualPad *pad1 = canvas->cd(1);
    pad1->SetGrid(1,1);
    h_ss->GetXaxis()->SetTitle("cos(#theta)");
    h_acc_ss->SetTitle("s#bar{s}");
    h_ss->SetFillColor(kBlue);
    h_ss->SetFillStyle(3001);
    h_ss->SetLineColor(kBlue);
    h_ss->SetLineWidth(1);
    h_acc_ss->GetYaxis()->SetRangeUser(0,9000);
    h_acc_ss->SetFillColor(kRed);
    h_acc_ss->SetFillStyle(3001);
    h_acc_ss->SetLineColor(kRed);
    h_acc_ss->SetLineWidth(1);
    // h_acc_ss->Rebin(4);
    // h_ss->Rebin(4);
    h_acc_ss->Draw("HIST");
    h_ss->Draw("HIST SAME");

    auto legend = new TLegend(0.1,0.8,0.3,0.9);
    legend->SetMargin(0.2);
    legend->AddEntry(h_ss, "KLPFO after pq-method", "f");
    legend->AddEntry(h_acc_ss, "KLPFO", "f");
    legend->Draw();

    canvas->cd(2);
    TVirtualPad *pad2 = canvas->cd(2);
    pad2->SetGrid(1,1);
    h_cc->GetXaxis()->SetTitle("cos(#theta)");
    h_acc_cc->SetTitle("c#bar{c}");
    h_cc->SetFillColor(kBlue);
    h_cc->SetFillStyle(3001);
    h_cc->SetLineColor(kBlue);
    h_cc->SetLineWidth(1);
    h_acc_cc->GetYaxis()->SetRangeUser(0,9000);
    h_acc_cc->SetFillColor(kRed);
    h_acc_cc->SetFillStyle(3001);
    h_acc_cc->SetLineColor(kRed);
    h_acc_cc->SetLineWidth(1);
    // h_acc_cc->Rebin(4);
    // h_cc->Rebin(4);
    h_acc_cc->Draw("HIST");
    h_cc->Draw("HIST SAME");

    auto legend1 = new TLegend(0.1,0.8,0.3,0.9);
    legend1->SetMargin(0.2);
    legend1->AddEntry(h_cc, "KLPFO after pq-method", "f");
    legend1->AddEntry(h_acc_cc, "KLPFO", "f");
    legend1->Draw();
}