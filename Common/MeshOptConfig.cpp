#include "MeshOptConfig.h"
#include "ConfigFileReader.h"
#include "CaseParameterList.h"
#include "BLParameterList.h"
#include "HighOrderParameterList.h"

int MeshOptConfig::Read(){
  ConfigFileReader reader;
  reader.ReadConfigFile("MeshOpt.config");
  all_parameter_type& params = reader.GetParameterList();
  auto it = params.find("CaseParameters");
  if(it != params.end()){
    case_parameters = std::make_shared<CaseParameterList>(it->second);
  }
  else{
    case_parameters = std::make_shared<CaseParameterList>();
  }

  it = params.find("BLParameters");
  if(it != params.end()){
    bl_parameters = std::make_shared<BLParameterList>(it->second);
  }
  else{
    bl_parameters = std::make_shared<BLParameterList>();
  }

  it = params.find("HighOrderParameters");
  if(it != params.end()){
    ho_parameters = std::make_shared<HighOrderParameterList>(it->second);
  }
  else{
    ho_parameters = std::make_shared<HighOrderParameterList>();
  }
}
