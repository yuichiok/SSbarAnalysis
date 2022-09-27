#include "FileSelector.hh"

using std::cout;   using std::endl;

FileSelector::FileSelector() {}

FileSelector::FileSelector(TString input)
:_full(input)
{
  SetNames(_full);
}

void FileSelector::SetNames(TString input)
{
  _full     = input;
  std::string str_full = _full.Data();
  std::string name_ext = str_full.substr(str_full.find_last_of("/") + 1);

  std::string::size_type const p(name_ext.find_last_of('.'));
  _name = TString(name_ext.substr(0, p));
}

TString FileSelector::GetOutName()
{
  return _name + _suffix;
}

TString FileSelector::GetOutName_withPath()
{
  return _out_path + _name + _suffix;
}