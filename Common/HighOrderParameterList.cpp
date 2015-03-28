#include "HighOrderParameterList.h"
#include "ParameterManip.h"

std::set<std::string> HighOrderParameterList::GetValidParameters() const{
  return {"order","targetMinQuality","useOptimizedPoints"};
}

void HighOrderParameterList::Populate(const parameter_map_type& pmap){
  ParameterManip::OverrideInteger(order,"order",pmap);
  ParameterManip::OverrideDouble(targetMinQuality,"targetMinQuality",pmap);
  ParameterManip::OverrideBool(useOptimizedPoints,"useOptimizedPoints",pmap);
}
