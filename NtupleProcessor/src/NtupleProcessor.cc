
/*------------------------------------------------------------------------------
NtupleProcessor.cpp
 Created : 2022-09-02  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <TString.h>
#include "../include/NtupleProcessor.hh"

using std::cout;   using std::endl;

NtupleProcessor::NtupleProcessor()
{
  // TEST output
    cout << "    NtupleProcessor: Created." << endl;


  // Takes input options and processes appropriately
  //   Options that can be specified:
  //     - Operating points
  //     - Differential variables (default: JetPT)
  //     - Dataset ntuple file to run on (default: all)
  //     - Desired number of entry to run on per file and entries to skip.
  //     - Modifications to the analysis (systematic analyses, etc)
  // For each dataset that needs to be run over...
  //   Opens the appropriate file and tree
  //   Creates a TreeIterator, EventHandler, and the desired HistogramMakers
  //   Runs the TreeIterator.

}