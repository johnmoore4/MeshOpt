#pragma once
#include "MyParameterList.h"

class BLParameterList: public MyParameterList{
 public:
  BLParameterList(){}
  BLParameterList(const parameter_map_type& pmap){ Initialize(pmap); }

  double thickness = 0.01;
  double Re = 1000;
  double minSpacing = 0.001;
  double growthRatio = 1.3;


  int distortionRatio = 2;
  int NLayers = 1;
  
  std::string flowType = "laminar";

 private:
  std::set<std::string> GetValidParameters() const;
  void Populate(const parameter_map_type& pmap);

  void Broadcast(){}
};
