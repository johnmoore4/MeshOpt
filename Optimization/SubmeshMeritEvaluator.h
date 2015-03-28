#pragma once
#include "MeritEvaluator.h"
#include "GlobalDefines.h"
#include "OptElManager.h"
#include <vector>

class MNode;
class MEl;
//class OptElManager;


class SubmeshMeritEvaluator: public MeritEvaluator{
 private:
  MNode* nd;
  const std::vector<ElIdealPair>& elements;
  OptElManager& optel_manager;
  std::vector<double> merits;
  const std::map<const MEl*,double>& all_merits;
  double factor;
  std::vector<bool> to_optimize;
 public:
 SubmeshMeritEvaluator(MNode* nd_t, 
		       const std::vector<ElIdealPair>& elements_t,
		       OptElManager& optel_manager_t,  
		       const std::map<const MEl*,double>& all_merits_t, 
		       double factor_t=1.0,
		       std::vector<bool> to_opt = {1,1,1}): 
  nd(nd_t), elements(elements_t), optel_manager(optel_manager_t), 
    all_merits(all_merits_t),
    factor(factor_t),
    to_optimize(to_opt){}
  const std::vector<double>& getMerits() const{ return merits; }

  //int NumDOFs(){ return 2; }

  int NumDOFs(){ return nd->getType(); }
  void GetCurrentState(double* state);
  void SetCurrentState(const double* state);
  double EvaluateGradient(double* gradient);
  /*
  double EvaluateHessian(double* gradient, double* Hessian,
			 gind* indi, gind* ingj){};
  */

};
