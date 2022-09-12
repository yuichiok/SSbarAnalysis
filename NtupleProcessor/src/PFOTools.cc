
/*------------------------------------------------------------------------------
PFOTools.cpp
 Created : 2022-09-11  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <TString.h>
#include <TFile.h> 
#include "../include/PFOTools.hh"

using std::cout;   using std::endl;

PFOTools::PFOTools( PFO_QQbar *ptr ) : data(ptr){}