#include "InterpolationPoints.h"
#include "assert.h"

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void JacobiGQ(arma::vec &x, arma::vec &w, double alpha, double beta, int N){

  arma::mat J(N+1,N+1,arma::fill::zeros);

  for(int i=0; i<=N; i++){
    double h1 = 2.0*i + alpha+beta;
    double jd1 = -1.0/((h1+2)*h1)*(alpha*alpha-beta*beta);
    if(alpha + beta > 10.0*std::numeric_limits<double>::epsilon()) J(i,i) = jd1;

    double ip1 = i+1;
    double jd2 = 2.0/(h1+2.0)*sqrt(ip1*(ip1+alpha+beta)*(ip1+alpha)*
				   (ip1+beta)/((h1+1.)*(h1+3.)));
    if(i < N){
      J(i,ip1) = jd2;
      J(ip1,i) = jd2;
    }
    

  }

  arma::mat eigvec;
  arma::eig_sym(x,eigvec,J);


  w = pow(eigvec.row(0).t(),2.0)*pow(2.0,alpha+beta+1.)/
    (alpha+beta+1.0)*factorial(alpha)*factorial(beta)/factorial(alpha+beta);
  

}

arma::vec jacobiGL(int alpha, int beta, int N){

  
  arma::vec xtemp, wtemp;
  JacobiGQ(xtemp,wtemp,alpha+1,beta+1,N-2);

  arma::vec x2(N+1);
  x2(0) = -1.0;
  x2.subvec(1,N-1) = xtemp;
  x2(N) = 1.0;


  return x2;

}



arma::mat UniformTriPts(int p){
  const int ndof = (p+1)*(p+2)/2;
  const int nel = p*p;

  arma::mat plocal(ndof,2);
  int cnt = 0;
  for(int i=0; i<=p; i++){
    for(int j=0; j<=(p-i); j++){
      plocal(cnt,0) = (double)j/p;
      plocal(cnt,1) =  (double)i/p;
      cnt++;
    }
  }
  

  return plocal;
}

arma::vec jacobiP(arma::vec x, double alpha, double beta, int p){
  const int N = x.n_elem;
  arma::mat PL(N,p+1);

  
  double gamma0 = pow(2.0,alpha+beta+1.0)/(alpha+beta+1.0)*factorial(alpha)
    *factorial(beta)/factorial(alpha+beta);
 
  
  PL.col(0) = 1.0/sqrt(gamma0)*arma::ones<arma::vec>(N);

 
  if(p==0) return PL.col(0);
  double gamma1 = (alpha+1.0)*(beta+1.0)/(alpha+beta+3.0)*gamma0;
  PL.col(1) = ( (alpha+beta+2.0)/2.0*x + (alpha-beta)/2.0)/
    sqrt(gamma1)%arma::ones<arma::vec>(N);
  if(p == 1) return PL.col(1);
  
  double aold = 2.0/(2.0+alpha+beta)*
    sqrt((alpha+1.0)*(beta+1.0)/(alpha+beta+3.0));

  for(int i=1; i<p; i++){
    double h1 = 2.0*i + alpha + beta;
    double anew = 2.0/(h1+2.0)*sqrt((i+1.0)*(i+1.0+alpha+beta)*(i+1.0+alpha)*
				    (i+1.0+beta)/((h1+1.0)*(h1+3.0)));
    double bnew = -(pow(alpha,2.0) - pow(beta,2.0))/(h1*(h1+2.0));
    PL.col(i+1) = 1.0/anew*( -aold*PL.col(i-1) + (x-bnew)%PL.col(i) );
    aold = anew;
  }
  return PL.col(p);
  
}

arma::mat vandermonde1D(int p, arma::vec r){
  const int N = r.n_elem;

  arma::mat V1D(N,p+1);
  for(int i=0; i<=p; i++){
    V1D.col(i) = jacobiP(r,0.,0.,i);
  }

  return V1D;
}

arma::vec warpFactor(int p,arma::vec rout){
  const int nr = rout.n_elem;

  arma::vec LGLr = jacobiGL(0,0,p);


  arma::vec req = arma::linspace<arma::vec>(-1.0,1.0,p+1);

  arma::mat Veq = vandermonde1D(p,req);
 

  arma::mat Pmat(nr,p+1);
  for(int i=0; i<=p; i++){
    Pmat.col(i) = jacobiP(rout,0,0,i);
  }
  
 
  arma::mat Lmat = solve(Veq.t(),Pmat.t());
  arma::vec warp = Lmat.t()*(LGLr-req);
  
  

  arma::vec zerof(nr,arma::fill::zeros);
  arma::uvec ind = find(abs(rout) < 1.0-1.0e-12);
  zerof(ind) = arma::ones<arma::vec>(ind.n_elem);
  arma::vec sf = 1.0 - pow(zerof%rout,2.0);
  warp = warp/sf + warp%(zerof-1.0);
  return warp;

}

arma::mat ElectrostaticTriPts(int p){
  
  arma::vec alpopt;
  alpopt << 0.0 << 0.0 << 1.4152 << 0.1001 << 0.2751 << 0.9800 << 1.0999 <<
    1.2832 << 1.3648 << 1.4773 << 1.4959 << 1.5743 << 1.5770 << 1.6223 << 
    1.6258;

  double alpha;
  if(p < 16) alpha = alpopt(p-1);
  else alpha = 5.0/3.0;

  const int ndof = (p+1)*(p+2)/2;

 

  arma::vec L1(ndof), L2(ndof), L3(ndof);
  int cnt = 0;
  for(int n=0; n<=p; n++){
    for(int m=0; m< p+1-n; m++){
      L1(cnt) = (double)n/(double)p;
      L3(cnt) = (double)m/(double)p;
      L2(cnt) = 1.0-L1(cnt)-L3(cnt);
      cnt++;
    }
  }
  
 
  
  arma::vec x = L3-L2;
  arma::vec y = (2.0*L1-L2-L3)/sqrt(3.0);

  if(p > 1){
    arma::vec warpf1 = warpFactor(p,L3-L2);
    arma::vec warpf2 = warpFactor(p,L1-L3);
    arma::vec warpf3 = warpFactor(p,L2-L1);
  
  
    arma::vec blend1 = 4.0*L2%L3;
    arma::vec blend2 = 4.0*L1%L3;
    arma::vec blend3 = 4.0*L1%L2;

 

    arma::vec warp1 = blend1%warpf1%(1.0 + pow(alpha*L1,2.0));
    arma::vec warp2 = blend2%warpf2%(1.0 + pow(alpha*L2,2.0));
    arma::vec warp3 = blend3%warpf3%(1.0 + pow(alpha*L3,2.0));
 
 
    double pi23 = 2.0/3.0*arma::datum::pi;
    double pi43 = 2.0*pi23;

    x+= 1.0*warp1 + cos(pi23)*warp2 + cos(pi43)*warp3;
    y+= 0.0*warp1 + sin(pi23)*warp2 + sin(pi43)*warp3;
  }


  arma::mat xy(ndof,2);
  xy.col(0) = x;
  xy.col(1) = y;

  arma::mat map(2,2);
  map(0,0) = 2.0; map(0,1) = 1.0;
  map(1,0) = 0.0; map(1,1) = sqrt(3);

  arma::mat invMap = inv(map);

  
  arma::mat r = xy;
  for(int i=0; i<r.n_rows; i++){
    r(i,0)+= 1.0;
    r(i,1)+= sqrt(3)/3.0;
  }
  xy = r*invMap.t();
  
  return xy;

}


arma::mat UniformQuadPts(const int p){
  const int ndof = (p+1)*(p+1);
  arma::mat pts(ndof,2);
  LinePoints Line;
  arma::mat linepts = Line.ComputePoints(p,0);
  
  for(int i=0,cnt=0; i<=p; i++){
    for(int j=0; j<=p; j++,cnt++){
      pts(cnt,0) = linepts(j);
      pts(cnt,1) = linepts(i);
    }
  }

  return pts;
}

arma::mat JacobiQuadPts(const int p){
  const int ndof = (p+1)*(p+1);
  arma::mat pts(ndof,2);
  LinePoints Line;
  arma::mat linepts = Line.ComputePoints(p,1);

  for(int i=0,cnt=0; i<=p; i++){
    for(int j=0; j<=p; j++,cnt++){
      pts(cnt,0) = linepts(j);
      pts(cnt,1) = linepts(i);
    }
  }

  return pts;
}


arma::mat UniformTetPts(int p){
  const int ndof = (p+1)*(p+2)*(p+3)/6;
  arma::mat pts(ndof,3);
  for(int k=0, cnt=0; k<=p; k++){
    for(int j=0; j<=(p-k); j++){
      for(int i=0; i<=(p-j-k); i++,cnt++){
	pts(cnt,0) = (double)i/(double)p;
	pts(cnt,1) = (double)j/(double)p;
	pts(cnt,2) = (double)k/(double)p;
      }
    }
  }

  return pts;
}

arma::vec evalwarp(const double p, const arma::vec& xnodes, const arma::vec& xout){
  const int nn = xout.n_elem;

  arma::vec warp = 0.0*xout;
  arma::vec xeq = 0.0*xout;

  for(int i=0; i<=p; i++){
    xeq(i) = -1.0 + 2.0*(p-i)/p;
  }

  for(int i=0; i<=p; i++){
    arma::vec d = (xnodes(i)-xeq(i))*arma::ones<arma::vec>(nn);
    for(int j=1; j<p; j++){
      if(i != j) d = d%(xout-xeq(j))/(xeq(i)-xeq(j));
    }
    if(i != 0) d = -d/(xeq(i)-xeq(0));
    if(i != p) d = d/(xeq(i)-xeq(p));
    
    warp+= d;

  }

  return warp;

}
void EvalShift(int p, double pval, const arma::vec& L1, const arma::vec& L2,
		 const arma::vec& L3,arma::vec &warpx, arma::vec &warpy){

  arma::vec gaussX = -jacobiGL(0.0,0.0,p);

  //gaussX.save("gaussX",raw_ascii);

  arma::vec blend1 = L2%L3;
  arma::vec blend2 = L1%L3;
  arma::vec blend3 = L1%L2;

  arma::vec warpfactor1 = 4.0*evalwarp(p,gaussX,L3-L2);
  arma::vec warpfactor2 = 4.0*evalwarp(p,gaussX,L1-L3);
  arma::vec warpfactor3 = 4.0*evalwarp(p,gaussX,L2-L1);

  
  arma::vec warp1 = blend1%warpfactor1%(1.0+pow(pval*L1,2.0));
  arma::vec warp2 = blend2%warpfactor2%(1.0+pow(pval*L2,2.0));
  arma::vec warp3 = blend3%warpfactor3%(1.0+pow(pval*L3,2.0));

  double pi = arma::datum::pi;
  warpx = warp1 + cos(2.0/3.0*pi)*warp2 + cos(4.0/3.0*pi)*warp3;
  warpy = sin(2.0/3.0*pi)*warp2 + sin(4.0/3.0*pi)*warp3;

}

arma::mat ElectrostaticTetPts(int p){
  
  arma::vec alpha_param = {0.0, 0.0, 0.0, 0.1002, 1.1332, 1.5608, 1.3413, 
			   1.2577, 1.1603, 1.10153, 0.6080, 0.4523, 0.8856, 
			   0.8717, 0.9655};
  

  double alpha;
  if(p < 16) alpha = alpha_param(p-1);
  else alpha = 1.0;
  
  const int ndof = (p+1)*(p+2)*(p+3)/6;
  
  arma::mat rst = UniformTetPts(p);
  if(p < 2) return rst;
  
  rst = 2.0*rst-1.0;

  //rst.save("rst",raw_ascii);

  arma::vec L1 = (1.0+rst.col(2))/2.0;
  arma::vec L2 = (1.0+rst.col(1))/2.0;
  arma::vec L3 = -(1.0+sum(rst,1))/2.0;
  arma::vec L4 = (1.0+rst.col(0))/2.0;

  arma::rowvec3 v1 = {-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0)};
  arma::rowvec3 v2 = {1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0)};
  arma::rowvec3 v3 = {0.0, 2.0/sqrt(3.0), -1.0/sqrt(6.0)};
  arma::rowvec3 v4 = {0.0, 0.0, 3.0/sqrt(6.0)};

  arma::mat::fixed<4,3> t1;
  arma::mat::fixed<4,3> t2;

  t1.row(0) = v2-v1;
  t1.row(1) = v2-v1;
  t1.row(2) = v3-v2;
  t1.row(3) = v3-v1;

  t2.row(0) = v3 - 0.5*(v1+v2);
  t2.row(1) = v4 - 0.5*(v1+v2);
  t2.row(2) = v4 - 0.5*(v2+v3);
  t2.row(3) = v4 - 0.5*(v1+v3);

  for(int i=0; i<4; i++){
    t1.row(i)/= norm(t1.row(i),2.0);
    t2.row(i)/= norm(t2.row(i),2.0);    
  }

  arma::mat XYZ = L3*v1 + L4*v2 + L2*v3 + L1*v4;
  arma::mat shift = 0.0*XYZ;

  arma::vec La, Lb, Lc, Ld;
  arma::vec warp1, warp2;
  for(int i=0; i<4; i++){
    if(i==0){
      La = L1; Lb = L2; Lc = L3; Ld = L4;
    }
    else if(i==1){
      La = L2; Lb = L1; Lc = L3; Ld = L4;
    }
    else if(i==2){
      La = L3; Lb = L1; Lc = L4; Ld = L2;
    }
    else{
      La = L4; Lb = L1; Lc = L3; Ld = L2;
    }
    

    arma::mat Lbcd(Lb.n_elem,3);
    Lbcd.col(0) = Lb;
    Lbcd.col(1) = Lc;
    Lbcd.col(2) = Ld;

    EvalShift(p,alpha,Lb,Lc,Ld,warp1,warp2);


    arma::vec blend = Lb%Lc%Ld;

    arma::vec denom = (Lb + 0.5*La)%(Lc + 0.5*La)%(Ld + 0.5*La);

    double tol = 1.0e-10;
    arma::uvec ind = find(denom > tol);
    blend(ind) = (1.0 + pow(alpha*La(ind),2.0))%blend(ind)/denom(ind);

    shift+= (blend%warp1)*t1.row(i) + (blend%warp2)*t2.row(i);


    int cnt=0;
    ind.clear();
    for(int j=0; j<La.n_elem; j++){
      if(La(j) < tol && ((Lb(j)>tol) + (Lc(j)>tol) + (Ld(j)>tol) ) < 3){
	ind.resize(cnt+1);
	ind(cnt) = j;
	cnt++;
      }
    }
    ind.resize(cnt);
    
    shift.rows(ind) = warp1(ind)*t1.row(i) + warp2(ind)*t2.row(i);
  }


  XYZ+= shift;

  arma::mat33 R;
  R(0,0) = 2.0; R(1,0) = 0.0;           R(2,0) = 0.0;
  R(0,1) = 1.0; R(1,1) = sqrt(3.0);     R(2,1) = 0.0;
  R(0,2) = 1.0; R(1,2) = 1.0/sqrt(3.);  R(2,2) = 4.0/(sqrt(2.)*sqrt(3.));
  
  
  XYZ.col(0)-= -1.0;
  XYZ.col(1)-= -1.0/sqrt(3.);
  XYZ.col(2)-= -1.0/(sqrt(3.)*sqrt(2.));

  XYZ*= arma::inv(R.t());
  

  return XYZ;

}

arma::mat LinePoints::ComputePoints(const int p, const int type) const{
  if(type == 0){ // Uniform
    return arma::linspace(0.0,1.0,p+1);
  }
  else if(type == 1){ // Jacobi roots
    if(p == 1) return arma::linspace(0.0,1.0,2);

    arma::vec pts1d = jacobiGL(0,0,p);
  
    return (pts1d+1.0)/2.0;
  }
  else assert(type);
}

arma::mat TriPoints::ComputePoints(const int p, const int type) const{
  if(type == 0){ // Uniform
    return UniformTriPts(p);
  }
  else if(type == 1){ // Electrostatic
    return ElectrostaticTriPts(p);
  }
  else if(type == 2){ // Electristatic Tet
    TetPoints Tet;
    arma::mat pts = Tet.ComputePoints(p,1);
    pts.resize((p+1)*(p+2)/2,2);
    return pts;
  }
  else assert(type);
}

arma::mat QuadPoints::ComputePoints(const int p, const int type) const{
  if(type == 0){ // Uniform
    return UniformQuadPts(p);
  }
  else if(type == 1){ // Electrostatic
    return JacobiQuadPts(p);
  }
  else assert(type);
}


arma::mat TetPoints::ComputePoints(const int p, const int type) const{
  if(type == 0){ // Uniform
    return UniformTetPts(p);
  }
  else if(type == 1){ // Electrostatic
    return ElectrostaticTetPts(p);
  }
  else assert(type);
}

arma::mat HexPoints::ComputePoints(const int p, const int type) const{
  arma::mat pts(std::pow(p+1,3),3);
  if(type == 0){ // Uniform
    for(int k = 0, cnt = 0; k<=p; k++){
      for(int j = 0; j<=p; j++){
	for(int i = 0; i<=p; i++, cnt++){
	  pts(cnt,0) = double(i)/p;
	  pts(cnt,1) = double(j)/p;
	  pts(cnt,2) = double(k)/p;
	}
      }
    }
  }
  else if(type == 1){ // Jacobi
    throw std::runtime_error("Optimal node distribution not yet implemented for Hexahedron!");
    
  }
  return pts;

}

arma::mat PrismPoints::ComputePoints(const int p, const int type) const{
  const int doftri = (p+1)*(p+2)/2;
  const int ndof = doftri*(p+1); 
  arma::mat pts(ndof,3);

  LinePoints Line;
  TriPoints Tri;

  arma::mat tripts, linepts;
  if(type == 0){
    tripts = Tri.ComputePoints(p,0);
    linepts = Line.ComputePoints(p,0);
  }
  else if(type == 1){
    tripts = Tri.ComputePoints(p,2);
    linepts = Line.ComputePoints(p,1);
  }
  else assert(type);

  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<doftri; i++,cnt++){
      pts(cnt,0) = tripts(i,0);
      pts(cnt,1) = tripts(i,1);
      pts(cnt,2) = linepts(j);
    }
  }
 
  return pts;

 
}

arma::mat PyramidPoints::ComputePoints(const int p, const int type) const{
  std::cout << "Begeinning pyramid points" << std::endl;
  int dofpyramid = (p+1)*(p+2)*(2*(p+1)+1)/6;
  arma::mat pts(dofpyramid,3);
  if(type == 0){
    for(int k = 0, cnt = 0; k <=p; k++){
      int pk = p-k;
      
      for(int j = 0; j <=pk; j++){
	for(int i = 0; i <=pk; i++, cnt++){
	  pts(cnt,0) = double(i)/p;
	  pts(cnt,1) = double(j)/p;
	  pts(cnt,2) = double(k)/p;
	}
      }
    }
  }
  else{
    throw std::runtime_error("Optimal node distribution not yet implemented for Pyramid!");
  }
  std::cout << "End pyramid points" << std::endl;
  return pts;

}
