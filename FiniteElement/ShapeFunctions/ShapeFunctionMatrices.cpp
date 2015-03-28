#include "ShapeFunctionMatrices.h"
#include "../InterpolationPoints/InterpolationPoints.h"
#include "../Polynomials/PolynomialBasis.h"
#include "../Quadrature/Quadrature.h"


ShapeFunctionMatrices::
ShapeFunctionMatrices(int order, 
			int interp_type,
			const PolynomialBasis& basis, 
			const InterpolationPoints& points,
			const Quadrature& quadrature,
			int el_type,
			int quad_degree):
  order(order), interp_type(interp_type), quad_degree(quad_degree), 
  el_type(el_type){

  arma::mat interp_pts = points.ComputePoints(order,interp_type);
  arma::mat gpoints;

  
  if(quad_degree == -1) quad_degree = 4*order;
  
  quadrature.ComputeQuadrature(quad_degree,gpoints,gweights);

  arma::mat A, f;
  arma::cube dA, df;
  arma::cube dA2, df2;
  basis.EvalBasis(order,interp_pts,A,dA,dA2);
  basis.EvalBasis(order,gpoints,f,df,df2);
 

  sfquad = arma::solve(A,f);
  
  sfquadT = sfquad.t();

  dsfquad.resize(f.n_rows,f.n_cols,df.n_slices);
  for(int i=0; i<df.n_slices; i++){
    dsfquad.slice(i) = arma::solve(A,df.slice(i));
  }
}

/*
void ShapeFunctionMatrices::GenerateMatrices(const int p, const int type,
					     const int el_type){
  using namespace std;

  arma::mat interp_pts = points.ComputePoints(p,type);
  arma::mat gpoints;

  int nq = 3*p;
  if(el_type != 4 && el_type != 2) nq = 4*p;

  quadrature.ComputeQuadrature(nq,gpoints,gweights);
  //quadrature.ComputeQuadrature(6*p-3,gpoints,gweights);
  //quadrature.ComputeQuadrature(3*p,gpoints,gweights);



  arma::mat A, f;
  arma::cube dA, df;
  arma::cube dA2, df2;
  basis.EvalBasis(p,interp_pts,A,dA,dA2);
  basis.EvalBasis(p,gpoints,f,df,df2);
 

  sfquad = arma::solve(A,f);
  
  sfquadT = sfquad.t();

  dsfquad.resize(f.n_rows,f.n_cols,df.n_slices);
  for(int i=0; i<df.n_slices; i++){
    dsfquad.slice(i) = arma::solve(A,df.slice(i));
  }


}
*/

ShapeFunctionMatricesFactory::ShapeFunctionMatricesFactory(){
  quadrature_factory = std::make_shared<QuadratureFactory>();
  
  interpolation_factory = std::make_shared<InterpolationPointFactory>();
  
  basis_factory = std::make_shared<PolynomialBasisFactory>();
}

const ShapeFunctionMatrices* ShapeFunctionMatricesFactory::
getShapeFunction(int el_type, int p, int interp_type, int quad_degree){
  auto& sf = ShapeFunctions[el_type][p];
  if(!sf){
    const PolynomialBasis& basis = basis_factory->getPolynomialBasis(el_type);

    const InterpolationPoints& points = 
      interpolation_factory->getInterpolationPoints(el_type);

    const Quadrature& quadrature = quadrature_factory->getQuadrature(el_type);

    sf = std::make_shared<ShapeFunctionMatrices>(p,interp_type,basis,points,
    						 quadrature,el_type,
    						 quad_degree);
  }

  return sf.get();

}
