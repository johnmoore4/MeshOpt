#include "NodeInterpolator.h"
#include "../InterpolationPoints/InterpolationPoints.h"
#include "../Polynomials/PolynomialBasis.h"

NodeInterpolator::NodeInterpolator(int interp_type): interp_type(interp_type){
  interpolation_factory = std::make_shared<InterpolationPointFactory>();
  basis_factory = std::make_shared<PolynomialBasisFactory>();
}

const arma::mat& NodeInterpolator::
getNodeInterpolationMatrix(int eltype, int pcurr, int pfinal){
  
  arma::mat& interp = interp_mats[eltype][pcurr][pfinal];
  
  if(interp.n_elem == 0){
    const PolynomialBasis& basis = basis_factory->getPolynomialBasis(eltype);

    const InterpolationPoints& points_computer = 
      interpolation_factory->getInterpolationPoints(eltype);

    arma::mat interp_curr = points_computer.ComputePoints(pcurr,interp_type);

    arma::mat interp_new = points_computer.ComputePoints(pfinal,interp_type);
    

    
  arma::mat A, f;
  arma::cube dA, df;
  arma::cube dA2, df2;
  basis.EvalBasis(pcurr,interp_curr,A,dA,dA2);
  basis.EvalBasis(pcurr,interp_new,f,df,df2);
 

  interp = arma::solve(A,f);
    
  }
  return interp;

}
