
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
            data->pfo_piddedx_p_dedxdist[ipfo],
            data->pfo_piddedx_pi_dedxdist[ipfo],
        };
        PFO.cos  = std::cos(PFO.vt.v3().Theta());
        PFO.qcos = (PFO.charge < 0) ? PFO.cos : -PFO.cos;

        PFO_jet[ijet].push_back(PFO);
    }

}

vector<PFO_Info> PFOTools::SortJet( vector<PFO_Info> jet )
{
    std::sort(jet.begin(), jet.end(),std::greater<PFO_Info>());
    return jet;
}

bool PFOTools::ValidPFO()
{
    bool isValid = false;
    for (int ijet=0; ijet < 2; ijet++) isValid = PFO_jet[ijet].size();
    return isValid;
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
