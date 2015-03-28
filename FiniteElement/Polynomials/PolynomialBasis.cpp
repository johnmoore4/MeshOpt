#include "PolynomialBasis.h"
#include "JacobiPolynomials.h"
#include "KoornwinderPolynomials.h"
#include "assert.h"

/*
// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FunctionSpaceToolsInPlace.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"
#include <Intrepid_PointTools.hpp>
#include <Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp>
#include <Intrepid_HGRAD_TRI_Cn_FEM_ORTH.hpp>
#include <Intrepid_HGRAD_PYR_C1_FEM.hpp>
#include <Intrepid_Basis.hpp>


// Shards includes
#include "Shards_CellTopology.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

using namespace Intrepid;
using Teuchos::RCP;


typedef FieldContainer<double> array_type;
*/

void LinePolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				    arma::mat& f, arma::cube& df,
				    arma::cube& df2) const{
  Jacobi1D(p,evalpts,f,df);  
}

void TriPolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				   arma::mat& f, arma::cube& df,
				   arma::cube& df2) const{
  Koornwinder2D koornwinder;
  koornwinder.EvaluatePolynomial(p,evalpts,f,df);

  /*
  const int ndof = (p+1)*(p+2)/2;
  const int neval = evalpts.n_rows;

  Basis_HGRAD_TRI_Cn_FEM_ORTH<double,array_type> tri_orth_basis(p);
  array_type eval_vals(ndof,neval);
  array_type eval_grads(ndof,neval,2);
  array_type eval_pts_array(neval,2);

  for(int i=0; i<neval; i++){
    for(int j=0; j<2; j++) eval_pts_array(i,j) = evalpts(i,j);
  }

 
  tri_orth_basis.getValues(eval_vals,eval_pts_array,OPERATOR_VALUE);
  tri_orth_basis.getValues(eval_grads,eval_pts_array,OPERATOR_GRAD);
  f.resize(ndof,neval);
  df.resize(ndof,neval,2);
  for(int i=0; i<ndof; i++){
    for(int j=0; j<neval; j++){
      f(i,j) = eval_vals(i,j);
      for(int k=0; k<2; k++){
	df(i,j,k) = eval_grads(i,j,k);
      }
    }
  }
  */
}

void QuadPolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				    arma::mat& f, arma::cube& df,
				    arma::cube& df2) const{

  
 
  arma::mat fi;
  arma::cube dfi;
  Jacobi1D(p,evalpts.unsafe_col(0),fi,dfi);
  
  arma::mat fj;
  arma::cube dfj;
  Jacobi1D(p,evalpts.unsafe_col(1),fj,dfj);

  const int neval = evalpts.n_rows;
  const int ndof = (p+1)*(p+1);
  f.resize(ndof,neval);
  df.resize(ndof,neval,2);
  
  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<=p; i++,cnt++){
      for(int k=0; k<neval; k++){
	f(cnt,k) = fi(i,k)*fj(j,k);
	df(cnt,k,0) = dfi(i,k,0)*fj(j,k);
	df(cnt,k,1) = dfj(j,k,0)*fi(i,k);
      }
    }
  }
}

void TetPolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				   arma::mat& f, arma::cube& df,
				   arma::cube& df2) const{
  
  Koornwinder3D koornwinder;
  koornwinder.EvaluatePolynomial(p,evalpts,f,df);
  
  /*
  using namespace std;
  const int ndof = (p+1)*(p+2)*(p+3)/6;
  const int neval = evalpts.n_rows;
  

  Basis_HGRAD_TET_Cn_FEM_ORTH<double,array_type> tet_orth_basis(p);
  array_type eval_vals(ndof,neval);
  array_type eval_grads(ndof,neval,3);;
  array_type eval_pts_array(neval,3);


  for(int i=0; i<neval; i++){
    for(int j=0; j<3; j++) eval_pts_array(i,j) = evalpts(i,j);
  }

 
  tet_orth_basis.getValues(eval_vals,eval_pts_array,OPERATOR_VALUE);
  tet_orth_basis.getValues(eval_grads,eval_pts_array,OPERATOR_GRAD);
  f.resize(ndof,neval);
  df.resize(ndof,neval,3);
  df2.resize(ndof,neval,9);
  for(int i=0; i<ndof; i++){
    for(int j=0; j<neval; j++){
      f(i,j) = eval_vals(i,j);
      for(int k=0; k<3; k++){
	df(i,j,k) = eval_grads(i,j,k);
      }
    }
  }
  */
}

void HexPolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				   arma::mat& f, arma::cube& df,
				   arma::cube& df2) const{

  const int neval = evalpts.n_rows;
  const int ndof = std::pow(p+1,3);

  arma::mat fi;
  arma::cube dfi;
  Jacobi1D(p,evalpts.unsafe_col(0),fi,dfi);
  
  arma::mat fj;
  arma::cube dfj;
  Jacobi1D(p,evalpts.unsafe_col(1),fj,dfj);

  arma::mat fk;
  arma::cube dfk;
  Jacobi1D(p,evalpts.unsafe_col(2),fk,dfk);

  f.resize(ndof,neval);
  df.resize(ndof,neval,3);
  
  for(int k = 0, cnt = 0; k<=p; k++){
    for(int j=0; j<=p; j++){
      for(int i=0; i<=p; i++,cnt++){
	for(int e=0; e<neval; e++){
	  f(cnt,e) = fi(i,e)*fj(j,e)*fk(k,e);
	  df(cnt,e,0) = dfi(i,e,0)*fj(j,e)*fk(k,e);
	  df(cnt,e,1) = fi(i,e)*dfj(j,e,0)*fk(k,e);
	  df(cnt,e,2) = fi(i,e)*fj(j,e)*dfk(k,e,0);
	}
      }
    }
  }

  
}

void PrismPolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				     arma::mat& f, arma::cube& df, 
				     arma::cube& df2) const{
  arma::mat trieval = evalpts.cols(0,1);

  //Koornwinder2D koornwinder;
  arma::mat ftri;
  arma::cube dftri;
  arma::cube temp;
  TriPolynomialBasis tri;
  tri.EvalBasis(p,trieval,ftri,dftri,temp);
  //koornwinder.EvaluatePolynomial(p,trieval,ftri,dftri);

  arma::mat lineeval = evalpts.cols(2,2);
  
  
  arma::mat fline;
  arma::cube dfline;
  Jacobi1D(p,lineeval,fline,dfline);



 
  const int doftri = (p+1)*(p+2)/2;
  const int dofprism = doftri*(p+1);
  const int neval = evalpts.n_rows;


  f.resize(dofprism,neval);
  df.resize(dofprism,neval,3);
 
  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<doftri; i++,cnt++){
      for(int k=0; k<neval; k++){
	f(cnt,k) = fline(j,k)*ftri(i,k);
	df(cnt,k,0) = dftri(i,k,0)*fline(j,k);
	df(cnt,k,1) = dftri(i,k,1)*fline(j,k);
	df(cnt,k,2) = ftri(i,k)*dfline(j,k,0);
      }
    }
  }
  

}

void PyramidPolynomialBasis::EvalBasis(const int p, const arma::mat& evalpts, 
				     arma::mat& f, arma::cube& df, 
				     arma::cube& dfold) const{


  arma::mat quadeval = evalpts.cols(0,1);
  arma::mat lineeval = evalpts.cols(2,2);



  LinePolynomialBasis line;
  QuadPolynomialBasis quad;
  
  arma::mat fline, fquad;
  arma::cube dfline, dfquad;
  arma::cube temp;

  const int neval = evalpts.n_rows;


  Jacobi1D(p,lineeval,fline,dfline);
  
  quad.EvalBasis(p,quadeval,fquad,dfquad,temp);
  

  int dofquad = std::pow(p+1,2);
  int dofline = p+1;
  int dofpyramid = (p+1)*(p+2)*(2*(p+1)+1)/6;
  
  f.resize(dofpyramid,neval);
  df.resize(dofpyramid,neval,3);
  

  for(int k = 0, cnt = 0; k <=p; k++){
    int pk = p-k;
    for(int j = 0; j <= pk; j++){
      for(int i = 0; i <= pk; i++, cnt++){
	for(int e = 0; e < neval; e++){
	  int qdof = i + j*(p+1);
	  f(cnt,e) = fquad(qdof,e)*fline(k,e);
	  df(cnt,e,0) = dfquad(qdof,e,0)*fline(k,e);
	  df(cnt,e,1) = dfquad(qdof,e,1)*fline(k,e);
	  df(cnt,e,2) = fquad(qdof,e)*dfline(k,e,0);
	}
      }
    }
  }
  

  /*
  int ndof = (p+1)*(p+2)*(2*(p+1)+1)/6;
  Basis_HGRAD_PYR_C1_FEM<double,array_type> Pyr_basis;
  array_type eval_vals(ndof,neval);
  array_type eval_grads(ndof,neval,3);;
  array_type eval_pts_array(neval,3);

  for(int i=0; i<neval; i++){
    for(int j=0; j<3; j++) eval_pts_array(i,j) = evalpts(i,j);
  }

  Pyr_basis.getValues(eval_vals,eval_pts_array,OPERATOR_VALUE);
  Pyr_basis.getValues(eval_grads,eval_pts_array,OPERATOR_GRAD);
  //arma::mat f2(ndof,neval);
  //arma::cube df2(ndof,neval,9);
  f.resize(ndof,neval);
  df.resize(ndof,neval,3);
  //df2.resize(ndof,neval,9);
  for(int i=0; i<ndof; i++){
    for(int j=0; j<neval; j++){
      f(i,j) = eval_vals(i,j);
      for(int k=0; k<3; k++){
	df(i,j,k) = eval_grads(i,j,k);
      }
    }
  }
  */

}
