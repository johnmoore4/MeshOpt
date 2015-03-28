#pragma once

#include <armadillo>
#include <memory>

class InterpolationPointFactory;
class PolynomialBasisFactory;

class NodeInterpolator{
 public:
  NodeInterpolator(int interp_type);
  const arma::mat& getNodeInterpolationMatrix(int eltype, int pcurr, 
					      int pfinal);

 private:
  arma::mat interp_mats[8][10][11];

  std::shared_ptr<InterpolationPointFactory> interpolation_factory;
  std::shared_ptr<PolynomialBasisFactory> basis_factory;

  const int interp_type;

};
