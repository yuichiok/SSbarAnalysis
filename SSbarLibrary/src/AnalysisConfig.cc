#include "AnalysisConfig.hh"

using std::cout;     using std::endl;   using std::stringstream;
using std::string;   using std::pair;   using std::vector;

AnalysisConfig::AnalysisConfig() {}

void AnalysisConfig::SetConfig(TString fnc)
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(fnc.Data(), pt);

  // Gen cuts
    gen_quark  = pt.get<int>("GENERATED.gen_quark");

  // PFO cuts
    PFO_TPCHits_max = pt.get<int>("PFO.PFO_TPCHits_max");
    PFO_p_min       = pt.get<float>("PFO.PFO_p_min");
    PFO_p_max       = pt.get<float>("PFO.PFO_p_max");
    PFO_offset_max  = pt.get<float>("PFO.PFO_offset_max");

  // LPFO cuts
    LPFO_p_min = pt.get<float>("LPFO.LPFO_p_min");
    LPFO_p_max = pt.get<float>("LPFO.LPFO_p_max");

    cout << LPFO_p_max << endl;

}