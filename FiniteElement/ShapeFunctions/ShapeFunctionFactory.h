#pragma once

#include "ShapeFunction.h"

#include "../InterpolationPoints/InterpolationPoints.h"
#include "../Polynomials/PolynomialBasis.h"
#include "../Quadrature/Quadrature.h"

#include <memory>



template<class ScalarType>
class ShapeFunctionFactory{
 public:
ShapeFunctionFactory(){
  quadrature_factory = std::make_shared<QuadratureFactory>();
  
  interpolation_factory = std::make_shared<InterpolationPointFactory>();
  
  basis_factory = std::make_shared<PolynomialBasisFactory>();

}
  const ShapeFunction<ScalarType>* 
    getShapeFunction(int el_type, 
		     int p, 
		     int interp_type,
		     int quad_degree=-1){
    
    //auto& sf = ShapeFunctions[el_type][p];
  

    if(!ShapeFunctions[el_type][p]){
      //if(!sf){
      std::cout << "making shape function" << std::endl;
      const PolynomialBasis& basis = basis_factory->getPolynomialBasis(el_type);

      const InterpolationPoints& points = 
	interpolation_factory->getInterpolationPoints(el_type);

      const Quadrature& quadrature = quadrature_factory->getQuadrature(el_type);

      std::cout << "before make_shared SF" << std::endl;
      
      ShapeFunctions[el_type][p] = 
	std::make_shared<ShapeFunction<ScalarType> >(p,interp_type,basis,points,
							quadrature,el_type,
							quad_degree);
      
      std::cout << "after make shared sf" << std::endl;
    }
    return ShapeFunctions[el_type][p].get();
    //return sf.get();
  }

 private:

  //std::vector<std::vector<std::shared_ptr<ShapeFunction<ScalarType> > > > 
  //  ShapeFunctions;
  std::shared_ptr<ShapeFunction<ScalarType> > ShapeFunctions[8][11];
  

  std::shared_ptr<QuadratureFactory> quadrature_factory;
  std::shared_ptr<InterpolationPointFactory> interpolation_factory;
  std::shared_ptr<PolynomialBasisFactory> basis_factory;

};
