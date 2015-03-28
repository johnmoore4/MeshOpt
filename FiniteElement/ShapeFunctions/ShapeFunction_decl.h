#pragma once

#include <armadillo>

class PolynomialBasis;
class InterpolationPoints;
class Quadrature;

template <class ScalarType>
class ShapeFunction{
 public:
  typedef arma::Mat<ScalarType> matType;
  typedef arma::Col<ScalarType> vecType;
  typedef arma::Cube<ScalarType> cubeType;

  ShapeFunction(int order, 
		int interp_type,
		const PolynomialBasis& basis, 
		const InterpolationPoints& points,
		const Quadrature& quadrature,
		int el_type,
		int quad_degree);

				  
  const matType& getQuadratureSF() const{ return sfquad; }
  //const cubeType& getQuadratureSFDeriv() const{ return dsfquad; }
  const matType& getSFTrans() const { return sfquadT; }

  const vecType& getQuadratureWeights() const{ return gweights; }
  
  const int getOrder() const{ return order; } 
  const int getInterpType() const{ return interp_type; }

 private:
  const int order;
  const int interp_type;
  const int quad_degree;
  const int el_type;

  matType sfquad;
  matType sfquadT;
  cubeType dsfquad;
  vecType gweights;
};
