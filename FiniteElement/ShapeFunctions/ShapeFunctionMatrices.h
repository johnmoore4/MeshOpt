#pragma once
#include <armadillo>
#include <memory>
#include <vector>

class QuadratureFactory;
class InterpolationPointFactory;
class PolynomialBasisFactory;

class PolynomialBasis;
class InterpolationPoints;
class Quadrature;

class ShapeFunctionMatrices{


 private:

  //const PolynomialBasis& basis;
  //const InterpolationPoints& points;
  //const Quadrature& quadrature;

  const int order;
  const int interp_type;
  const int quad_degree;
  const int el_type;

  arma::mat sfquad;
  arma::mat sfquadT;
  arma::cube dsfquad;
  
  //arma::cube dsf2;
  arma::vec gweights;

    
 public:
 
  ShapeFunctionMatrices(int order, 
			int interp_type,
			const PolynomialBasis& basis, 
			const InterpolationPoints& points,
			const Quadrature& quadrature,
			int el_type,
			int quad_degree);

  /*
 ShapeFunctionMatrices(int orderT, int type, const PolynomialBasis& basisT,
		       const InterpolationPoints& pointsT, const Quadrature&
		       quadratureT, const int el_type): 
  order(orderT), type(type), basis(basisT), points(pointsT), 
    quadrature(quadratureT){
      //GenerateMatrices(order,type,el_type);
    }
  */
  const arma::mat& getQuadratureSF() const{ return sfquad; }
  const arma::cube& getQuadratureSFDeriv() const{ return dsfquad; }
  const arma::mat& getSFTrans() const { return sfquadT; }

  const arma::vec& getQuadratureWeights() const{ return gweights; }
  
  //void GenerateMatrices(const int p, const int type, const int el_type);
  const int getOrder() const{ return order; } 
  const int getInterpType() const{ return interp_type; }
};


class ShapeFunctionMatricesFactory{
 public:
  ShapeFunctionMatricesFactory();
  const ShapeFunctionMatrices* getShapeFunction(int el_type, int p, 
						int interp_type,
						int quad_degree=-1);

 private:
  std::shared_ptr<ShapeFunctionMatrices> ShapeFunctions[8][11];
  

  std::shared_ptr<QuadratureFactory> quadrature_factory;
  std::shared_ptr<InterpolationPointFactory> interpolation_factory;
  std::shared_ptr<PolynomialBasisFactory> basis_factory;
  
};


/*
class ShapeFunctionMatricesFactory{
 private:
  std::vector<std::unique_ptr<ShapeFunctionMatrices> > ShapeFunctions;
  QuadratureFactory quadrature_factory;
  InterpolationPointFactory interpolation_factory;
  PolynomialBasisFactory basis_factory;

 public:
  ShapeFunctionMatricesFactory(){ ShapeFunctions.resize(8); }
  const ShapeFunctionMatrices* getShapeFunction(const int el_type, const int p,
						const int interp_type){
    using namespace std;
    if(!ShapeFunctions[el_type]){
      const PolynomialBasis& basis = basis_factory.getPolynomialBasis(el_type);
      const InterpolationPoints& points = 
	interpolation_factory.getInterpolationPoints(el_type);
      const Quadrature& quadrature = quadrature_factory.getQuadrature(el_type);
      ShapeFunctions[el_type] = std::unique_ptr<ShapeFunctionMatrices>
      	(new ShapeFunctionMatrices(p,interp_type,basis,points,quadrature,el_type));

    }
    if(ShapeFunctions[el_type]->getOrder() != p ||
      ShapeFunctions[el_type]->getInterpType() != interp_type ){
      ShapeFunctions[el_type]->GenerateMatrices(p,interp_type,el_type);
    }
    return ShapeFunctions[el_type].get();
    
  }

};
*/
