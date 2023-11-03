#ifndef GUARD_HistManager_h
#define GUARD_HistManager_h

#include <iostream>
#include <map>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h> 

const int NSLABS = 15;

class HistManager
{
  public:
    HistManager();
    ~HistManager(){}

  // Methods
  virtual void InitializeHists();
  virtual void Hist2List();
  virtual void WriteLists( TFile * output );

  // Declear histograms
  // h1 hist

  enum class h_primary{
    d0,
    z0,
    d0_sigma,
    z0_sigma,
    d0_sigma_d0,
    z0_sigma_z0,
    dummy,
    last = dummy 
  };
  TH2F * h_primary[static_cast<int>(h_primary::last)];

  enum class h_secondary{
    d0,
    z0,
    d0_sigma,
    z0_sigma,
    d0_sigma_d0,
    z0_sigma_z0,
    dummy,
    last = dummy 
  };
  TH2F * h_secondary[static_cast<int>(h_secondary::last)];

  enum class h_garbage{
    d0,
    z0,
    d0_sigma,
    z0_sigma,
    d0_sigma_d0,
    z0_sigma_z0,
    dummy,
    last = dummy 
  };
  TH2F * h_garbage[static_cast<int>(h_garbage::last)];

  enum class h_total{
    d0,
    z0,
    d0_sigma,
    z0_sigma,
    d0_sigma_d0,
    z0_sigma_z0,
    dummy,
    last = dummy 
  };
  TH2F * h_total[static_cast<int>(h_total::last)];

  enum class h_general{
    cuts,
    ctag_primary,
    ctag_secondary,
    ctag_total,
    dummy,
    last = dummy 
  };
  TH1F * h_general[static_cast<int>(h_general::last)];

  enum class h_cos{
    cos_theta,
    cos_theta_acc,
    cos_theta_rej,
    cos_theta_primary,
    cos_theta_acc_primary,
    cos_theta_rej_primary,
    dummy,
    last=dummy
  };
  TH2F * h_cos[static_cast<int>(h_cos::last)];

  private:
    TList* hList_primary   = new TList();
    TList* hList_secondary = new TList();
    TList* hList_garbage   = new TList();
    TList* hList_total     = new TList();
    TList* hList_general   = new TList();
    TList* hList_cos       = new TList();
};

#endif
