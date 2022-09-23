
/*------------------------------------------------------------------------------
PFOTools.cpp
 Created : 2022-09-11  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <algorithm>
#include <TString.h>
#include <TFile.h> 

#include "PFOTools.hh"
#include "VectorTools.hh"

using std::cout;   using std::endl;

template<typename T>
void pop_front(std::vector<T>& vec)
{
    assert(!vec.empty());
    vec.front() = std::move(vec.back());
    vec.pop_back();
}

PFOTools::PFOTools() {}

PFOTools::PFOTools( PFO_QQbar *ptr )
: data(ptr)
{
    InitializePFOTools( data );
}

void PFOTools::InitializePFOTools( PFO_QQbar *data )
{

  for (int ipfo=0; ipfo < data->pfo_n; ipfo++)
  {
    // Make sure PFO belongs to either jet0 or jet 1
    if (data->pfo_match[ipfo] < 0 || 1 < data->pfo_match[ipfo]) continue;

    // Make suer PFO has only one reconstructed track to avoid (lambda/sigma)
    if (data->pfo_ntracks[ipfo] != 1) continue;

    VectorTools vt(data->pfo_px[ipfo], data->pfo_py[ipfo], data->pfo_pz[ipfo], data->pfo_E[ipfo]);

    // Initialize & categorize PFO variables with different jets
    const int ijet = data->pfo_match[ipfo];
    PFO = {
      vt,
      (Float_t) vt.v3().R(),
      data->pfo_E[ipfo],
      data->pfo_charge[ipfo],
      data->pfo_pdgcheat[ipfo],
      data->pfo_tpc_hits[ipfo],
      sqrt(data->pfo_d0[ipfo] * data->pfo_d0[ipfo] + data->pfo_z0[ipfo] * data->pfo_z0[ipfo]),
      data->pfo_dedx[ipfo],
      data->pfo_piddedx_k_dedxdist[ipfo],
      data->pfo_piddedx_pi_dedxdist[ipfo],
      data->pfo_piddedx_p_dedxdist[ipfo],
    };
    PFO.dEdx_dist_pdg = Get_dEdx_dist_PID( PFO.kdEdx_dist, PFO.pidEdx_dist, PFO.pdEdx_dist );
    PFO.cos           = std::cos(PFO.vt.v3().Theta());
    PFO.qcos          = (PFO.charge < 0) ? PFO.cos : -PFO.cos;

    PFO_jet[ijet].push_back(PFO);
    
  }

  if( ValidPFO() ){
    for (int ijet=0; ijet < 2; ijet++){
      LPFO[ijet]  = GetSortedJet(ijet).at(0);
      if( PFO_jet[ijet].size() > 1 ){
        SPFOs[ijet] = GetSortedJet(ijet);
        // SPFOs[ijet].erase(SPFOs[ijet].begin());
        pop_front(SPFOs[ijet]); // faster algorithm wise?

        std::copy_if(SPFOs[ijet].begin(), SPFOs[ijet].end(), std::back_inserter(SPFOs_K[ijet]), [](PFO_Info iPFO) {
            return isKaon(iPFO);
        });

        std::copy_if(SPFOs[ijet].begin(), SPFOs[ijet].end(), std::back_inserter(SPFOs_Pi[ijet]), [](PFO_Info iPFO) {
            return isPion(iPFO);
        });

      }
    }
  }

}

vector<PFO_Info> PFOTools::SortJet( vector<PFO_Info> jet )
{
    std::sort(jet.begin(), jet.end(),std::greater<PFO_Info>());
    return jet;
}

Bool_t PFOTools::ValidPFO()
{
    for ( auto iPFO_jet : PFO_jet ) {
        if( !iPFO_jet.size() ) return false;
    }
    return true;
}

vector<PFO_Info> PFOTools::GetJet( int ijet )
{
    return PFO_jet[ijet];
}

vector<PFO_Info> PFOTools::GetSortedJet( int ijet )
{
    vector<PFO_Info> sorted_jet = PFO_jet[ijet];
    sorted_jet = SortJet( sorted_jet );
    return sorted_jet;
}

Int_t PFOTools::Get_dEdx_dist_PID( Float_t kdEdx_dist, Float_t pidEdx_dist, Float_t pdEdx_dist )
{
    const static int nparticles = 3;
    Int_t        k_pi_p[nparticles] = { 321, 211, 2212 };
    Float_t k_pi_p_dEdx[nparticles] = { abs(kdEdx_dist), abs(pidEdx_dist), abs(pdEdx_dist) };

    Float_t  min_dEdx_dist = *std::min_element(k_pi_p_dEdx, k_pi_p_dEdx + nparticles );
    Int_t   imin_dEdx_dist =  std::find( k_pi_p_dEdx, k_pi_p_dEdx + nparticles, min_dEdx_dist ) - k_pi_p_dEdx;

    return k_pi_p[imin_dEdx_dist];

}

Bool_t PFOTools::isKaon( PFO_Info iPFO )
{
    return iPFO.dEdx_dist_pdg == 321;
}

Bool_t PFOTools::isPion( PFO_Info iPFO )
{
    return iPFO.dEdx_dist_pdg == 211;
}

Bool_t PFOTools::isProton( PFO_Info iPFO )
{
    return iPFO.dEdx_dist_pdg == 2212;
}

Bool_t PFOTools::is_charge_config( ChargeConfig cc )
{

  Int_t charge_config = LPFO[0].charge * LPFO[1].charge;

  switch (cc)
  {
    case kSame:
      if ( charge_config > 0 ) return true;
      break;
    
    case kOpposite:
      if ( charge_config < 0 ) return true;
      break;

    default:
      return false;
      break;
  }

  return false;

}

Bool_t PFOTools::PFO_Quality_checks( PFO_Info iPFO )
{
  vector<Bool_t> CutTrigger;

  CutTrigger.push_back( is_momentum( iPFO, 20.0, 60.0 ) );     // MIN/MAX momentum check
  CutTrigger.push_back( is_tpc_hits( iPFO, 210 ) );            // Number of TPC hit check
  CutTrigger.push_back( is_offset_small( iPFO, 1.0 ) );        // Offset distance check
  
  for (auto trigger : CutTrigger ){
    if (!trigger) { return false; }
  }

  return true;

}

Bool_t PFOTools::is_momentum( PFO_Info iPFO, Float_t MINP_CUT, Float_t MAXP_CUT )
{
  return ( MINP_CUT < iPFO.p_mag && iPFO.p_mag < MAXP_CUT);
}

Bool_t PFOTools::is_tpc_hits( PFO_Info iPFO, Int_t MIN_TPC_HITS )
{
  return ( MIN_TPC_HITS < iPFO.tpc_hits );
}

Bool_t PFOTools::is_offset_small( PFO_Info iPFO, Int_t MAX_OFFSET )
{
  return ( iPFO.pv < MAX_OFFSET );
}