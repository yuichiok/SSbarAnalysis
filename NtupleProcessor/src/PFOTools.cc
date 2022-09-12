
/*------------------------------------------------------------------------------
PFOTools.cpp
 Created : 2022-09-11  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <TString.h>
#include <TFile.h> 
#include "../include/PFOTools.hh"

using std::cout;   using std::endl;

PFOTools::PFOTools( PFO_QQbar *ptr )
: data(ptr)
{
    for(int ipfo=0; ipfo < data->pfo_n; ipfo++)
    {
        // Make sure PFO belongs to either jet0 or jet 1
        if (data->pfo_match[ipfo] < 0 || data->pfo_match[ipfo] == 2)
            continue;
        // Make suer PFO has only one reconstructed track to avoid (lambda/sigma)
        if (data->pfo_ntracks[ipfo] != 1)
            continue;

        // Initialize PFO variables
        PFO_Info PFO = {
            TVector3(data->pfo_px[ipfo], data->pfo_py[ipfo], data->pfo_pz[ipfo]),
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
        PFO.qcos = (PFO.chg < 0) ? PFO.cos : -PFO.cos;



    }

}