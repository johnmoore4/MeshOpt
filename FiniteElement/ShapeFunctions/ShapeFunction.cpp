#include "ShapeFunction.h"
//#include "../InterpolationPoints/InterpolationPoints.h"
//#include "../Polynomials/PolynomialBasis.h"
//#include "../Quadrature/Quadrature.h"

template class ShapeFunction<double>;
template class ShapeFunction<float>;

/*
template<class T>
int ShapeFunction<T>::
computeDoubleShapeFunction(arma::mat& sf, 
			   arma::cube& dsf, 
			   arma::vec& gw,
			   const PolynomialBasis& basis, 
			   const InterpolationPoints& points,
			   const Quadrature& quadrature){

  arma::mat interp_pts = points.ComputePoints(order,interp_type);
  arma::mat gpoints;

  
  if(quad_degree == -1) quad_degree = 4*order;
  
  quadrature.ComputeQuadrature(quad_degree,gpoints,gw);

  arma::mat A, f;
  arma::cube dA, df;
  arma::cube dA2, df2;
  basis.EvalBasis(order,interp_pts,A,dA,dA2);
  basis.EvalBasis(order,gpoints,f,df,df2);
 

  sf = arma::solve(A,f);
  
  //sfquadT = sfquad.t();

  dsf.resize(f.n_rows,f.n_cols,df.n_slices);
  for(int i=0; i<df.n_slices; i++){
    dsf.slice(i) = arma::solve(A,df.slice(i));
  }

}
*/
