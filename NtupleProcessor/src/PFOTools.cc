
/*------------------------------------------------------------------------------
PFOTools.cpp
 Created : 2022-09-11  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <algorithm>
#include <TString.h>
#include <TFile.h> 
#include "../include/PFOTools.hh"

using std::cout;   using std::endl;

PFOTools::PFOTools( PFO_QQbar *ptr )
: data(ptr)
{
    int test[2] = {0};
    for(int ipfo=0; ipfo < data->pfo_n; ipfo++)
    {
        // Make sure PFO belongs to either jet0 or jet 1
        if (data->pfo_match[ipfo] < 0 || 1 < data->pfo_match[ipfo]) continue;

        // Make suer PFO has only one reconstructed track to avoid (lambda/sigma)
        if (data->pfo_ntracks[ipfo] != 1) continue;

        TVector3 p3(data->pfo_px[ipfo], data->pfo_py[ipfo], data->pfo_pz[ipfo]);

        // Initialize & categorize PFO variables with different jets
        const int ijet = data->pfo_match[ipfo];
        PFO = {
            p3,
            p3.Mag(),
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
        PFO.cos  = PFO.p3.CosTheta();
        PFO.qcos = (PFO.charge < 0) ? PFO.cos : -PFO.cos;

        PFO_jet[ijet].push_back(PFO);

    }

    cout << "Jet PFO (" << data->pfo_n << ") size" << endl;
    cout << "ext0="<< test[0] << ", ext1=" << test[1] << endl;
    cout << PFO_jet[0].size() << endl;
    cout << PFO_jet[1].size() << endl;

}

void PFOTools::PFOSort( vector<PFO_Info> jet )
{
    std::sort(jet.begin(), jet.end(),std::greater<PFO_Info>());
}

PFO_Info PFOTools::GetLPFO( int ijet )
{
    PFOSort( PFO_jet[ijet] );
    return PFO_jet[ijet].at(0);
}