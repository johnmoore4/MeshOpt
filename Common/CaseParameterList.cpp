#include "CaseParameterList.h"
#include "ParameterManip.h"

std::set<std::string> CaseParameterList::GetValidParameters() const{
  return {"gmshFileName","stepFileName","newFileName"};
}

void CaseParameterList::Populate(const parameter_map_type& pmap){
  ParameterManip::OverrideString(gmshFileName,"gmshFileName",pmap);
  ParameterManip::OverrideString(stepFileName,"stepFileName",pmap);
  ParameterManip::OverrideString(newFileName,"newFileName",pmap);
}
