#include "JacobiPolynomials.h"

using namespace arma;


void polyval(const arma::vec &c, const arma::vec &x, arma::vec& val){
  const int p = c.n_elem-1;
  const int ne = x.n_elem;
  val.zeros();
  for(int i=0; i<=p; i++) val+= c(i)*pow(x,i);

}

void polyder(const arma::vec &c, arma::vec& dc){
  const int nc = c.n_elem-1;
  if(dc.n_elem != nc) dc.resize(nc);

  for(int i=1; i<=nc; i++) dc(i-1) = i*c(i);
}

arma::vec jacobi(const int n, const double a, const double b){
  vec pp(n+1);
  vec p0(1);
  vec p1(2);
  vec2 tmp;
  vec p2, q;

  p0(0) = 1.0;
  p1(0) = (a+b+2.0)/2.0;
  p1(1) = (a-b)/2.0;
  if(n==0){
    pp(0) = 1.0;
  }
  else if(n==1){
    pp(0) = p1(0);
    pp(1) = p1(1);
  }
  else{
    for(int i=1; i<n; i++){
      
        const double a1 = 2.*(i+1.)*(i+a+b+1.)*(2.*i+a+b);
        const double a2 = (2.*i+a+b+1.)*(a*a-b*b);
        const double a3 = (2.*i+a+b)*(2.*i+a+b+1.)*(2.*i+a+b+2.);
        const double a4 = 2.*(i+a)*(i+b)*(2.*i+a+b+2.);
	tmp[0] = a3;
	tmp[1] = a2;
	p2  = conv(tmp,p1);
	q.zeros(i+2);
	for(int j=2; j<=i+1; j++) q(j) = a4*p0(j-2);
	//q.subvec(2,2+i-1) = a4*p0;
	pp = (p2-q)/a1;
	p0 = p1;
	p1 = pp;
    }
  }
  if(n>0){
    for(int i=0; i<n+1; i++) pp(i) = p1(n-i);
  }
  return pp;
}

void Jacobi1D(const int p, const arma::mat& xin, arma::mat &f, 
	      arma::cube& df){
  const arma::mat x = 2.0*xin-1.0;
  const int nx = x.n_rows;

  f.resize(p+1,nx);
  df.resize(p+1,nx,1);


  vec pval(nx);
  vec dpval(nx);
  vec dpp, pp;
  for(int i=0; i<=p; i++){
    pp = jacobi(i,0.0,0.0);
    pp = pp*(double)sqrt(2*i+1);
    polyder(pp,dpp);
    polyval(pp,x,pval);
    polyval(dpp,x,dpval);
    for(int j=0; j<nx; j++){
      f(i,j) = pval(j);
      df(i,j,0) = 2.0*dpval(j);
    }
  }
}

