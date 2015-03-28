#pragma once
#include "MeritEvaluator.h"
#include "GlobalDefines.h"
#include "OptElManager.h"
//#include "TrilinosTypedefs.h"

#include <vector>
#include <unordered_set>
#include <map>




class MNode;
class MEl;
//class OptElManager;


class LocalMeshMeritEvaluator: public MeritEvaluator{
 private:
 
  const std::map<gind,gind>& activenodes;
  //const std::unordered_set<MEl*>& elements;
  const std::vector<const MEl*>& elements;
  const std::map<const MEl*,arma::mat>& idealElements;
  OptElManager& optel_manager;
  std::vector<double> merits;
  const double factor;
  std::vector<bool> to_optimize;

 public:
 LocalMeshMeritEvaluator(const std::map<gind,gind>& activenodes_t,
			 const std::vector<const MEl*>& elements_t,
			 const std::map<const MEl*,arma::mat>& ideal_t,
			 OptElManager& optel_manager_t,
			 const double factor_t=1.0,
			 std::vector<bool> to_opt = {1,1,1}): 
  activenodes(activenodes_t), elements(elements_t), idealElements(ideal_t),
    optel_manager(optel_manager_t), factor(factor_t), to_optimize(to_opt){}

  const std::vector<double>& getMerits() const{ return merits; }

  int NumDOFs();
  void GetCurrentState(double* state);
  void SetCurrentState(const double* state);
  double EvaluateGradient(double* gradient);

  //typedef Tpetra::MultiVector<ST,LO,GO, node_type> multivector_type;
  //typedef Tpetra::CrsMatrix<ST,LO,GO,node_type> sparse_mat_type;

  //double EvaluateHessian(Teuchos::RCP<sparse_mat_type>& HESS, 
  //			 Teuchos::RCP<multivector_type>& GRAD);

  double EvaluateHessianDense(double* gradient, double* Hessian);

};
