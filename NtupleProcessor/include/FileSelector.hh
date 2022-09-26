#ifndef GUARD_FileSelector_h
#define GUARD_FileSelector_h

#include <iostream>
#include <string>
#include <TString.h>
#include <TFile.h> 
#include <TObjArray.h>
#include <TObjString.h>
#include <map>

using std::string;

class FileSelector
{
  public:
    FileSelector();
    FileSelector(string o);
    ~FileSelector(){}

    virtual void   SetNames(string o);
    virtual string GetOutName();
    virtual string GetOutName_withPath();

  private:
    string _full;
    string _name_ext;
    string _name;
    string _suffix = ".ss.tmp.root";

    string _out_path = "rootfiles/tmp_root/";


};

#endif