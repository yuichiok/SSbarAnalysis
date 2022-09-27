#include "FileSelector.hh"

using std::cout;   using std::endl;
using std::string;

FileSelector::FileSelector() {}

FileSelector::FileSelector(string o)
:_full(o)
{
  SetNames(o);
}

void FileSelector::SetNames(string o)
{
  _full     = o;
  _name_ext = _full.substr(_full.find_last_of("/") + 1);

  string::size_type const p(_name_ext.find_last_of('.'));
  _name = _name_ext.substr(0, p);
}

string FileSelector::GetOutName()
{
  return _name + _suffix;
}

string FileSelector::GetOutName_withPath()
{
  return _out_path + _name + _suffix;
}