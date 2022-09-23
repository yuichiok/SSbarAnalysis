
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

PFOTools::PFOTools( PFO_QQbar *ptr )
: data(ptr)
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
    Float_t k_pi_p_dEdx[nparticles] = { kdEdx_dist, pidEdx_dist, pdEdx_dist };

    Float_t  min_dEdx_dist = *std::min_element(k_pi_p_dEdx, k_pi_p_dEdx + nparticles );
    Int_t   imin_dEdx_dist =  std::find( k_pi_p_dEdx, k_pi_p_dEdx + nparticles, min_dEdx_dist ) - k_pi_p_dEdx;

    return k_pi_p[imin_dEdx_dist];

}