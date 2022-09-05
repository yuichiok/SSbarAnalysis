
/*------------------------------------------------------------------------------
EventAnalyzer.cpp
 Created : 2022-09-05  okugawa
------------------------------------------------------------------------------*/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>
#include "../include/EventAnalyzer.hh"

using std::cout;   using std::endl;

EventAnalyzer::EventAnalyzer()
{
  // TEST output
    cout << "    NtupleProcessor: Created." << endl;

}