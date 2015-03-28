#include "BLParameterList.h"
#include "ParameterManip.h"

std::set<std::string> BLParameterList::GetValidParameters() const{
  return {"thickness","Re","distortionRatio","NLayers","flowType",
      "minSpacing","growthRatio"};
}

void BLParameterList::Populate(const parameter_map_type& pmap){
  ParameterManip::OverrideDouble(thickness,"thickness",pmap);
  ParameterManip::OverrideDouble(Re,"Re",pmap);
  ParameterManip::OverrideDouble(minSpacing,"minSpacing",pmap);
  ParameterManip::OverrideDouble(growthRatio,"growthRatio",pmap);
  ParameterManip::OverrideInteger(distortionRatio,"distortionRatio",pmap);
  ParameterManip::OverrideInteger(NLayers,"NLayers",pmap);
  ParameterManip::OverrideString(flowType,"flowType",pmap);
}
