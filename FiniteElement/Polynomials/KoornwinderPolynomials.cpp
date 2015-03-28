#include "KoornwinderPolynomials.h"
#include "JacobiPolynomials.h"

arma::umat pascalindex(const int p){
  const int ndof = (p+1)*(p+2)/2;

  arma::umat pind(ndof,2);
  for(int i=0, l=0; i<=p; i++){
    for(int j=0; j<=i; j++){
      pind(l,0) = i-j;
      pind(l,1) = j;
      l++;
    }
  }
  return pind;
}

arma::imat tetindex(const int p){
  const int ndof = (p+1)*(p+2)*(p+3)/6;
  arma::imat index(ndof,3);
  
  for(int i=0, cnt=0; i<=p; i++){
    for(int j=0; j<=i; j++){
      for(int k=0; k<=j; k++, cnt++){
	index(cnt,0) = i-j;
	index(cnt,1) = j-k;
	index(cnt,2) = k;
      }
    }
  }
  return index;
}

void Koornwinder2D::EvaluatePolynomial(const int p, const arma::mat& xin, 
				       arma::mat &f, arma::cube &df){

  using namespace arma;

 const int ndof = (p+1)*(p+2)/2;
 const int npt = xin.n_rows;

 arma::mat x = 2.*xin-1.;
 vec xi = x.col(0);
 vec eta = x.col(1);

 f.resize(ndof,npt);
 df.resize(ndof,npt,2);


 const umat pind = pascalindex(p);
 arma::vec r = 2.*(1.0+xi)/(1.-eta)-1.;
 const arma::vec s = eta;
  for(int i=0; i<npt; i++){
    if(x(i,1)==1){
      r(i) = -1.0;
    }
  }


  const arma::vec dr_dxi = 2./(1.0-x.col(1));
  const arma::vec dr_deta = 2.*(1.0+x.col(0))/pow((1.-x.col(1)),2.0);
  const arma::vec2 cv = {0.5,-0.5};

  arma::vec pval(npt), qval(npt), dpval(npt), dqval(npt);
  arma::vec dpp,dqp;
  arma::vec pp, qp;
  for(int i=0; i<ndof; i++){
    pp = jacobi(pind(i,0),0.0,0.0);
    qp = jacobi(pind(i,1),2.*pind(i,0)+1.,0.0);

    for(int j=0; j<pind(i,0); j++){
      qp = conv(cv,qp);
    }
    
    
    polyval(pp,r,pval);
    polyval(qp,s,qval);

    const double fc = sqrt((2.*pind(i,0)+1.)*2.*(pind(i,0)+pind(i,1)+1.));

  
    polyder(pp,dpp);
    polyder(qp,dqp);

    polyval(dpp,r,dpval);
    polyval(dqp,s,dqval);

    for(int j=0; j<npt; j++){
      f(i,j) = fc*pval(j)*qval(j);
      df(i,j,0) = 2.0*fc*dpval(j)*qval(j)*dr_dxi(j);
      df(i,j,1) = 2.0*fc*(dpval(j)*qval(j)*dr_deta(j) + pval(j)*dqval(j));
    }
  }
  


 
}

void Koornwinder3D::EvaluatePolynomial(const int p, const arma::mat& xin, 
				       arma::mat &f, arma::cube &df){

  using namespace arma;

  const int ndof = (p+1)*(p+2)*(p+3)/6;
  const int npts = xin.n_rows;
  
  f.resize(ndof,npts);
  df.resize(ndof,npts,3);



  const mat x = 2.*xin-1.;
  const vec xi = x.unsafe_col(0);
  const vec eta = x.unsafe_col(1);
  const vec zeta = x.unsafe_col(2);

  vec r = -2.*(1.+xi)/(eta+zeta) - 1.;
  vec s = 2.*(1.+eta)/(1.-zeta) - 1.;
  vec t = zeta;

  for(int i=0; i<r.n_rows; i++){
    if(std::abs(zeta(i)-1.0) < 1.0e-10) s(i) = -1.0;
    if(std::abs(eta(i)+zeta(i)) < 1.0e-10) r(i) = -1.0;
 
  }


  const vec dr_dxi = -2./(eta+zeta);
  const vec dr_deta = 2.*(1.+xi)/pow(eta+zeta,2.0);
  const vec dr_dzeta = 2.*(1.+xi)/pow(eta+zeta,2.0);
  
  const vec ds_deta = 2./(1. - zeta);
  const vec ds_dzeta = -2.*(1. + eta)/pow(1. - zeta,2.0);


  const imat mnl = tetindex(p);



  const vec2 pm_half = {0.5, -0.5}; 
  vec pm, pn, pl;
  vec dpm, dpn, dpl;
  vec mval(npts), nval(npts), lval(npts);
  vec dmval(npts), dnval(npts), dlval(npts);

  for(int i=0; i<ndof; i++){
    const int m = mnl(i,0);
    const int n = mnl(i,1);
    const int l = mnl(i,2);
    pm = jacobi(m,0.,0.);
    pn = jacobi(n,2.*m+1.,0.);
    pl = jacobi(l,2.*(m+n+1.),0.);
    for(int j=0; j<m; j++){
      pn = conv(pm_half,pn);
    }
    for(int j=0; j<(m+n); j++){
      pl = conv(pm_half,pl);
    }

    polyder(pm,dpm);
    polyder(pn,dpn);
    polyder(pl,dpl);

    polyval(pm,r,mval);
    polyval(pn,s,nval);
    polyval(pl,t,lval);

    polyval(dpm,r,dmval);
    polyval(dpn,s,dnval);
    polyval(dpl,t,dlval);
    

    double fc = sqrt((2.*m+1.)*(2.*n+1+2.*m+1.)*(2.*l+1+2.*(m+n+1.)));

    for(int j=0; j<npts; j++){
      f(i,j) = fc*mval(j)*nval(j)*lval(j);
      df(i,j,0) = 2.0*fc*dmval(j)*dr_dxi(j)*nval(j)*lval(j);
      df(i,j,1) = 2.0*fc*(dmval(j)*dr_deta(j)*nval(j)*lval(j) + mval(j)*
			  dnval(j)*ds_deta(j)*lval(j));
      df(i,j,2) = 2.0*fc*(dmval(j)*dr_dzeta(j)*nval(j)*lval(j) + mval(j)*
			  (dnval(j)*ds_dzeta(j)*lval(j) + nval(j)*dlval(j)));
									
    }
  }
}
