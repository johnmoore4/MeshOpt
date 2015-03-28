#pragma once
#include "../InterpolationPoints/InterpolationPoints.h"
#include "../Polynomials/PolynomialBasis.h"
#include "../Quadrature/Quadrature.h"

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
	       int quad_degree):
  order(order), interp_type(interp_type), quad_degree(quad_degree), 
    el_type(el_type){

      arma::mat interp_pts = points.ComputePoints(order,interp_type);
      arma::mat gpoints;

  
      if(quad_degree == -1) quad_degree = 4*order;
  
      arma::vec gw;
      quadrature.ComputeQuadrature(quad_degree,gpoints,gw);
      
      arma::mat A, f;
      arma::cube dA, df;
      arma::cube dA2, df2;
      basis.EvalBasis(order,interp_pts,A,dA,dA2);
      basis.EvalBasis(order,gpoints,f,df,df2);
 

      arma::mat sf;
      //arma::cube dsf;

      sf = arma::solve(A,f);
  
      //sfquadT = sfquad.t();

      //dsf.resize(f.n_rows,f.n_cols,df.n_slices);
      arma::mat dsf[3];
      for(int i=0; i<df.n_slices; i++){
	//dsf.slice(i) = arma::solve(A,df.slice(i)).t();
	dsf[i] = arma::solve(A,df.slice(i));
      }

      sfquad = arma::conv_to<matType>::from(sf);
      sfquadT = sfquad.t();
      gweights = arma::conv_to<vecType>::from(gw);
      dsfquad.resize(df.n_slices);
      dsfquadTrans.resize(dsfquad.size());
      for(int i = 0; i < df.n_slices; i++){
	dsfquad[i] = arma::conv_to<matType>::from(dsf[i]);
	dsfquadTrans[i] = dsfquad[i].t();
      }
    }


				  
  const matType& getQuadratureSF() const{ return sfquad; }
  //const cubeType& getQuadratureSFDeriv() const{ return dsfquad; }
  const std::vector<matType>& getQuadratureSFDeriv() const{ return dsfquad; }
  const std::vector<matType>& getQuadratureSFDerivTrans() const{ 
    return dsfquadTrans; 
  }
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
  std::vector<matType> dsfquad;
  std::vector<matType> dsfquadTrans;
  //cubeType dsfquad;
  vecType gweights;
};
