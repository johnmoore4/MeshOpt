#pragma once
#include <memory>

class CaseParameterList;
class BLParameterList;
class HighOrderParameterList;

class MeshOptConfig{
 public:
  int Read();
  std::shared_ptr<CaseParameterList>& getCaseParameters(){ 
    return case_parameters; 
  }
  std::shared_ptr<BLParameterList>& getBLParameters(){
    return bl_parameters;
  }
  std::shared_ptr<HighOrderParameterList>& getHighOrderParameters(){
    return ho_parameters;
  }

 private:
  std::shared_ptr<CaseParameterList> case_parameters;
  std::shared_ptr<BLParameterList> bl_parameters;
  std::shared_ptr<HighOrderParameterList> ho_parameters;
  
};
