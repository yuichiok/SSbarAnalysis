#ifndef GUARD_FileSelector_h
#define GUARD_FileSelector_h

#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <map>

class FileSelector
{
  public:
    FileSelector();
    FileSelector(TString o="");
    ~FileSelector(){}

  private:
    TString _full;
    TString _path;
    TString _name;
    TString _ext;



};

#endif