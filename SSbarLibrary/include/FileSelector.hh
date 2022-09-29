#ifndef GUARD_FileSelector_h
#define GUARD_FileSelector_h

#include <iostream>
#include <string>
#include <TString.h>
#include <TFile.h> 
#include <TObjArray.h>
#include <TObjString.h>
#include <map>

class FileSelector
{
  public:
    FileSelector();
    FileSelector(TString input);
    ~FileSelector(){}

    void   SetNames(TString input);
    TString GetOutName();
    TString GetOutName_withPath();

  private:
    TString _full;
    TString _name;
    TString _suffix = ".ss.tmp.root";

    TString _out_path = "rootfiles/tmp_root/";


};

#endif