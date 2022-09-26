#include <iostream>
#include <TString.h>
#include <TFile.h> 

#include "FileSelector.hh"

using std::cout;   using std::endl;

FileSelector::FileSelector(){}

FileSelector::FileSelector(TString o)
:_full(o)
{
  cout << _full << endl;


}

