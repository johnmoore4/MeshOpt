
#include "Quadrature.h"
#include "JacobiPolynomials.h"
#include "quad.h"

#include <armadillo>

#define PHG_TET_RULE(x) QUAD_3D_P##x##_

QUAD getPHGQuad(int p){
  if(p == 1) return PHG_TET_RULE(1);
  if(p == 2) return PHG_TET_RULE(2);
  if(p == 3) return PHG_TET_RULE(3);
  if(p == 4) return PHG_TET_RULE(4);
  if(p == 5) return PHG_TET_RULE(5);
  if(p == 6) return PHG_TET_RULE(6);
  if(p == 7) return PHG_TET_RULE(7);
  if(p == 8) return PHG_TET_RULE(8);
  if(p == 9) return PHG_TET_RULE(9);
  if(p == 10) return PHG_TET_RULE(10);
  if(p == 11) return PHG_TET_RULE(11);
  if(p == 12) return PHG_TET_RULE(12);
  if(p == 13) return PHG_TET_RULE(13);
  if(p == 14) return PHG_TET_RULE(14);
}

arma::vec roots(arma::vec c){
  int n = c.n_elem-1;

  c = flipud(c);

  arma::mat A(n,n,arma::fill::zeros);
  for(int i=0; i<n-1; i++) A(i+1,i) = 1.0;
  A.row(0) = (-c.subvec(1,n)/c(0)).t();

  arma::cx_vec eval = arma::eig_gen(A);

  return real(eval);


}

void generateFromOrbit(const arma::mat &p0, const arma::vec &w0, const int np1,
		       const int np2, const int np3,arma::mat &x,arma::mat &w){

  using namespace arma;
  double cospi23 = cos(2.0*datum::pi/3.0);
  double sinpi23 = sin(2.0*datum::pi/3.0);
  
  double cospi43 = cos(4.0*datum::pi/3.0);
  double sinpi43 = sin(4.0*datum::pi/3.0);

  mat A = {
    cospi23, -sinpi23,
    sinpi23, cospi23};
  A.reshape(2,2);

  mat B = {
    cospi43, -sinpi43,
    sinpi43, cospi43};
  B.reshape(2,2);
  
  mat C(2,2,fill::zeros);
  C(0,0) = 1.0;
  C(1,1) = -1.0;

  mat D = {
    cospi23, -sinpi23, 
    -sinpi23, -cospi23};
  D.reshape(2,2);

  mat E = eye<mat>(2,2);

  mat F = {
    cospi43, -sinpi43, 
    -sinpi43, -cospi43};
  F.reshape(2,2);


  const int npts = np1 + 3*np2 + 6*np3;
  x.resize(2,npts);
  w.resize(npts);

  for(int i=0; i<np1; i++){
    x.col(0) = p0.col(0);
    w(0) = w0(0);
  }

  for(int i=0; i<np2; i++){
    const int start = np1;
    x.col(0+start+i*3) = E.t()*p0.col(start+i);
    x.col(1+start+i*3) = D.t()*p0.col(start+i);
    x.col(2+start+i*3) = F.t()*p0.col(start+i);
    for(int j=0; j<3; j++){
      w(j+start+i*3) = w0(start+i);
    }
  }

  for(int i=0; i<np3; i++){
    const int start = np1 + 3*np2;
    x.col(0+start+i*6) = A.t()*p0.col(np1+np2+i);
    x.col(1+start+i*6) = B.t()*p0.col(np1+np2+i);  
    x.col(2+start+i*6) = C.t()*p0.col(np1+np2+i);  
    x.col(3+start+i*6) = D.t()*p0.col(np1+np2+i);  
    x.col(4+start+i*6) = E.t()*p0.col(np1+np2+i);  
    x.col(5+start+i*6) = F.t()*p0.col(np1+np2+i);  
    for(int j=0; j<6; j++){
      w(j+start+i*6) = w0(np1+np2+i);
    }
  }
  
  mat::fixed<2,2> T;
  T(0,0) =  1.5;
  T(0,1) = 0.0;
  T(1,0) = sqrt(3.0)/2.0;
  T(1,1) = sqrt(3.0);

  x.row(0)+= 0.5;
  x.row(1)+= sqrt(3.0)/2.0;

  x = inv(T)*x;

  x = x.t();

  w/=2.0;
  


}

void LineQuadrature::ComputeQuadrature(const int p, arma::mat& x, arma::vec& w)
const{
  using namespace arma;
  using namespace std;

 
  const int n = ceil((p+1)/2.);
  vec P = jacobi(n,0.,0.);

 
  x = arma::sort(roots(P));
 

  mat A(n,n,fill::zeros);
  arma::vec temp(x.n_elem);
  for(int i=0; i<n; i++){
    vec P2 = jacobi(i,0.,0.);
    polyval(P2,x,temp);
    A.unsafe_col(i) = temp;
  }
 
  vec rhs(n,fill::zeros);
  rhs(0) = 2.0;
  w = arma::solve(A.t(),rhs);
  
  x = (x+1.0)/2.0;
  w = w/2.0;

}

void TriQuadrature::ComputeQuadrature(const int p, arma::mat& x, arma::vec& w)
const{

  using namespace arma;
  
  if(p<=1){
    x = {0.3333333333333333, 0.3333333333333333};
    w = {0.5};
    //x.reshape(2,1);
    
  }
  else if(p == 2){
    x = {0.6666666666666666, 0.166666666666666, 
	      0.166666666666666, 0.6666666666666666,
	      0.166666666666666, 0.166666666666666};
    w = {0.166666666666666, 0.166666666666666, 0.166666666666666};
    x.reshape(2,3);
    x = x.t();
  }
  
  else if(p  == 3){
    
    x = {0.3333333333333333, 0.3333333333333333,
	      0.6, 0.2,
	      0.2, 0.6,
	      0.2, 0.2};
    w = {-0.2812500000000000, 0.2604166666666667, 0.2604166666666667,  0.2604166666666667};

    x.reshape(2,4);
    x = x.t();
  }
  else if(p<=4){
    x = {0.1081030181680700,  0.4459484909159650,
	      0.4459484909159650, 0.1081030181680700,
	      0.4459484909159650, 0.4459484909159650,
	      0.8168475729804590, 0.0915762135097710,
	      0.0915762135097710, 0.8168475729804590,
	      0.0915762135097710, 0.0915762135097710};
    w = {0.1116907948390055, 0.1116907948390055, 0.1116907948390055, 0.0549758718276610, 0.0549758718276610, 0.0549758718276610};
    x.reshape(2,6);
    x = x.t();
  }
  else if(p <= 5){
    mat p0 = {
      0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4104261923153453E+00,   0.0000000000000000E+00,
      0.6961404780296310E+00,   0.0000000000000000E+00  };
 
    
    vec w0 = {
      0.2250000000000000E+00,
      0.1323941527885062E+00,
      0.1259391805448271E+00 };
   
    p0.reshape(2,w0.n_elem);
    generateFromOrbit(p0,w0,1,2,0,x,w);

    //cout << "x: " << endl;
    //cout << x << endl;


  }
  else if(p <=10){
    mat p0 = {
      0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4935962988634245E+00,   0.0000000000000000E+00,
      -0.2840373491871686E+00,   0.0000000000000000E+00,
      0.4457307617703263E+00,   0.0000000000000000E+00,
      0.9385563442849673E+00,   0.0000000000000000E+00,
      -0.4474955151540920E+00,  -0.5991595522781586E+00,
      -0.4436763946123360E+00,  -0.2571781329392130E+00  };
 


    vec w0 = {
      0.8352339980519638E-01,
      0.7229850592056743E-02,
      0.7449217792098051E-01,
      0.7864647340310853E-01,
      0.6928323087107504E-02,
      0.2951832033477940E-01,
      0.3957936719606124E-01 };
    
    p0.reshape(2,w0.n_elem);
    generateFromOrbit(p0,w0,1,4,2,x,w);


    
  }
  else if(p <=15){
    mat p0 = {
      -0.3748423891073751E+00,   0.0000000000000000E+00,
      -0.2108313937373917E+00,   0.0000000000000000E+00,
      0.1204084962609239E+00,   0.0000000000000000E+00,
      0.5605966391716812E+00,   0.0000000000000000E+00,
      0.8309113970031897E+00,   0.0000000000000000E+00,
      0.9502746194248890E+00,   0.0000000000000000E+00,
      -0.4851316950361628E+00,  -0.4425551659467111E+00,
      -0.4762943440546580E+00,  -0.1510682717598242E+00,
      -0.4922845867745440E+00,  -0.6970224211436132E+00,
      -0.4266165113705168E+00,  -0.5642774363966393E+00,
      -0.3968468770512212E+00,  -0.3095105740458471E+00,
      -0.2473933728129512E+00,  -0.2320292030461791E+00  };

 

    vec w0 = {
      0.3266181884880529E-01,
      0.2741281803136436E-01,
      0.2651003659870330E-01,
      0.2921596213648611E-01,
      0.1058460806624399E-01,
      0.3614643064092035E-02,
      0.8527748101709436E-02,
      0.1391617651669193E-01,
      0.4291932940734835E-02,
      0.1623532928177489E-01,
      0.2560734092126239E-01,
      0.3308819553164567E-01 };

    p0.reshape(2,w0.n_elem);
    generateFromOrbit(p0,w0,0,6,6,x,w);
  }
  else if(p <= 20){
    mat p0 = {
      0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4977490260133565E+00,   0.0000000000000000E+00,
      -0.3587903720915737E+00,   0.0000000000000000E+00,
      -0.1932918138657104E+00,   0.0000000000000000E+00,
      0.2064993924016380E+00,   0.0000000000000000E+00,
      0.3669431077237697E+00,   0.0000000000000000E+00,
      0.6767931784861860E+00,   0.0000000000000000E+00,
      0.8827927364865920E+00,   0.0000000000000000E+00,
      0.9664768608120111E+00,   0.0000000000000000E+00,
      -0.4919755727189941E+00,  -0.7513212483763635E+00,
      -0.4880677744007016E+00,  -0.5870191642967427E+00,
      -0.4843664025781043E+00,  -0.1717270984114328E+00,
      -0.4835533778058150E+00,  -0.3833898305784408E+00,
      -0.4421499318718065E+00,  -0.6563281974461070E+00,
      -0.4466292382741727E+00,  -0.6157647932662624E-01,
      -0.4254937754558538E+00,  -0.4783124082660027E+00,
      -0.4122204123735024E+00,  -0.2537089901614676E+00,
      -0.3177533194934086E+00,  -0.3996183176834929E+00,
      -0.2889337325840919E+00,  -0.1844183967233982E+00  };

  

    vec w0 = {
      0.2761042699769952E-01,
      0.1779029547326740E-02,
      0.2011239811396117E-01,
      0.2681784725933157E-01,
      0.2452313380150201E-01,
      0.1639457841069539E-01,
      0.1479590739864960E-01,
      0.4579282277704251E-02,
      0.1651826515576217E-02,
      0.2349170908575584E-02,
      0.4465925754181793E-02,
      0.6099566807907972E-02,
      0.6891081327188203E-02,
      0.7997475072478163E-02,
      0.7386134285336024E-02,
      0.1279933187864826E-01,
      0.1725807117569655E-01,
      0.1867294590293547E-01,
      0.2281822405839526E-01 };

    p0.reshape(2,w0.n_elem);
    generateFromOrbit(p0,w0,1,8,10,x,w);

  }
  else if(p <= 25){

    mat p0 = {
      -0.4580802753902387E+00,   0.0000000000000000E+00,
      -0.3032320980085228E+00,   0.0000000000000000E+00,
      -0.1696674057318916E+00,   0.0000000000000000E+00,
      0.1046702979405866E+00,   0.0000000000000000E+00,
      0.2978674829878846E+00,   0.0000000000000000E+00,
      0.5455949961729473E+00,   0.0000000000000000E+00,
      0.6617983193620190E+00,   0.0000000000000000E+00,
      0.7668529237254211E+00,   0.0000000000000000E+00,
      0.8953207191571090E+00,   0.0000000000000000E+00,
      0.9782254461372029E+00,   0.0000000000000000E+00,
      -0.4980614709433367E+00,  -0.4713592181681879E+00,
      -0.4919004480918257E+00,  -0.1078887424748246E+00,
      -0.4904239954490375E+00,  -0.3057041948876942E+00,
      -0.4924576827470104E+00,  -0.7027546250883238E+00,
      -0.4897598620673272E+00,  -0.7942765584469995E+00,
      -0.4849757005401057E+00,  -0.5846826436376921E+00,
      -0.4613632802399150E+00,  -0.4282174042835178E+00,
      -0.4546581528201263E+00,  -0.2129434060653430E+00,
      -0.4542425148392569E+00,  -0.6948910659636692E+00,
      -0.4310651789561460E+00,  -0.5691146659505208E+00,
      -0.3988357991895837E+00,  -0.3161666335733065E+00,
      -0.3949323628761341E+00,  -0.1005941839340892E+00,
      -0.3741327130398251E+00,  -0.4571406037889341E+00,
      -0.3194366964842710E+00,  -0.2003599744104858E+00,
      -0.2778996512639500E+00,  -0.3406754571040736E+00,
      -0.2123422011990124E+00,  -0.1359589640107579E+00  };

    vec w0 = {
      0.8005581880020417E-02,
      0.1594707683239050E-01,
      0.1310914123079553E-01,
      0.1958300096563562E-01,
      0.1647088544153727E-01,
      0.8547279074092100E-02,
      0.8161885857226492E-02,
      0.6121146539983779E-02,
      0.2908498264936665E-02,
      0.6922752456619963E-03,
      0.1248289199277397E-02,
      0.3404752908803022E-02,
      0.3359654326064051E-02,
      0.1716156539496754E-02,
      0.1480856316715606E-02,
      0.3511312610728685E-02,
      0.7393550149706484E-02,
      0.7983087477376558E-02,
      0.4355962613158041E-02,
      0.7365056701417832E-02,
      0.1096357284641955E-01,
      0.1174996174354112E-01,
      0.1001560071379857E-01,
      0.1330964078762868E-01,
      0.1415444650522614E-01,
      0.1488137956116801E-01 };

    p0.reshape(2,w0.n_elem);
    generateFromOrbit(p0,w0,0,10,16,x,w);

  }
  else if(p <= 30){
    mat p0 = {
      0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4890048253508517E+00,   0.0000000000000000E+00,
      -0.3755064862955532E+00,   0.0000000000000000E+00,
      -0.2735285658118844E+00,   0.0000000000000000E+00,
      -0.1461412101617502E+00,   0.0000000000000000E+00,
      0.1570364626117722E+00,   0.0000000000000000E+00,
      0.3179530724378968E+00,   0.0000000000000000E+00,
      0.4763226654738105E+00,   0.0000000000000000E+00,
      0.6302247183956902E+00,   0.0000000000000000E+00,
      0.7597473133234094E+00,   0.0000000000000000E+00,
      0.8566765977763036E+00,   0.0000000000000000E+00,
      0.9348384559595755E+00,   0.0000000000000000E+00,
      0.9857059671536891E+00,   0.0000000000000000E+00,
      -0.4986119432099803E+00,  -0.1459114994581331E+00,
      -0.4979211112166541E+00,  -0.7588411241269780E+00,
      -0.4944763768161339E+00,  -0.5772061085255766E+00,
      -0.4941451648637610E+00,  -0.8192831133859931E+00,
      -0.4951501277674842E+00,  -0.3331061247123685E+00,
      -0.4902988518316453E+00,  -0.6749680757240147E+00,
      -0.4951287867630010E+00,  -0.4649148484601980E+00,
      -0.4869873637898693E+00,  -0.2747479818680760E+00,
      -0.4766044602990292E+00,  -0.7550787344330482E+00,
      -0.4730349181194722E+00,  -0.1533908770581512E+00,
      -0.4743136319691660E+00,  -0.4291730489015232E+00,
      -0.4656748919801272E+00,  -0.5597446281020688E+00,
      -0.4508936040683500E+00,  -0.6656779209607333E+00,
      -0.4492684814864886E+00,  -0.3024354020045064E+00,
      -0.4466785783099771E+00,  -0.3733933337926417E-01,
      -0.4241903145397002E+00,  -0.4432574453491491E+00,
      -0.4144779276264017E+00,  -0.1598390022600824E+00,
      -0.4037707903681949E+00,  -0.5628520409756346E+00,
      -0.3792482775685616E+00,  -0.3048723680294163E+00,
      -0.3434493977982042E+00,  -0.4348816278906578E+00,
      -0.3292326583568731E+00,  -0.1510147586773290E+00,
      -0.2819547684267144E+00,  -0.2901177668548256E+00,
      -0.2150815207670319E+00,  -0.1439403370753732E+00    };

    vec w0 = {
      0.1557996020289920E-01,
      0.3177233700534134E-02,
      0.1048342663573077E-01,
      0.1320945957774363E-01,
      0.1497500696627150E-01,
      0.1498790444338419E-01,
      0.1333886474102166E-01,
      0.1088917111390201E-01,
      0.8189440660893461E-02,
      0.5575387588607785E-02,
      0.3191216473411976E-02,
      0.1296715144327045E-02,
      0.2982628261349172E-03,
      0.9989056850788964E-03,
      0.4628508491732533E-03,
      0.1234451336382413E-02,
      0.5707198522432062E-03,
      0.1126946125877624E-02,
      0.1747866949407337E-02,
      0.1182818815031657E-02,
      0.1990839294675034E-02,
      0.1900412795035980E-02,
      0.4498365808817451E-02,
      0.3478719460274719E-02,
      0.4102399036723953E-02,
      0.4021761549744162E-02,
      0.6033164660795066E-02,
      0.3946290302129598E-02,
      0.6644044537680268E-02,
      0.8254305856078458E-02,
      0.6496056633406411E-02,
      0.9252778144146602E-02,
      0.9164920726294280E-02,
      0.1156952462809767E-01,
      0.1176111646760917E-01,
      0.1382470218216540E-01 };

    p0.reshape(2,w0.n_elem);
    generateFromOrbit(p0,w0,1,12,23,x,w);

    //x.save("gausspts.dat",raw_ascii);

  }
  else{
    vec x1, w1;
    LineQuadrature Line;
    Line.ComputeQuadrature(p,x1,w1);

    x1 = 2.0*x1-1.0;
    w1 = 2.0*w1;
    const int ng1d = x1.n_elem;

    mat x2(ng1d,ng1d);
    for(int i=0; i<ng1d; i++){
      x2.col(i) = x1;
    }
    mat y2 = x2.t();
    x2 = vectorise(x2);
    y2 = vectorise(y2);
    x.resize(ng1d*ng1d,2);
    x.col(0) = (1.0 + x2 - y2 - x2%y2)/4.0;
    x.col(1) = (1.0 + y2)/2.0;

    w.resize(x.n_rows);
    w = vectorise(w1*w1.t());

    w = w%(1.0-y2)/8.0;
  }
 
}

void QuadQuadrature::ComputeQuadrature(int p, arma::mat &x, arma::vec &w) const{

  arma::vec gp1d, gw1d;
  LineQuadrature Line;
  Line.ComputeQuadrature(p,gp1d,gw1d);
  const int ng1d = gw1d.n_elem;
  const int ngq = pow(ng1d,2);

  w.resize(ngq);
  x.resize(ngq,2);
  for(int i=0,cnt=0; i<ng1d; i++){
    for(int j=0; j<ng1d; j++,cnt++){
      x(cnt,0) = gp1d(j);
      x(cnt,1) = gp1d(i);
      w(cnt) = gw1d(j)*gw1d(i);
    }
  }

}

void TetQuadrature::ComputeQuadrature(int p, arma::mat &x, arma::vec &w) const{

  using namespace arma;
  QUAD quad = getPHGQuad(p);

  if(p <= 14){
    const int NP = quad.npoints;
    w.resize(NP);
    w = arma::vec(quad.weights,NP)/6.0;
  
    arma::mat bary(quad.points,4,NP);
  
  
    arma::mat barynodes = {0,0,0,
			   1,0,0,
			   0,1,0,
			   0,0,1};

    //std::cout << "sum weights: " << arma::accu(w) << std::endl;
    barynodes.reshape(3,4);
    arma::mat temp = barynodes*bary;
    //std::cout << barynodes << std::endl;
    //std::cout << temp << std::endl;

    x = temp.t();
    return;
  }
  else{
    arma::mat x1d;
    arma::vec w1d;

 
    LineQuadrature Line;
    Line.ComputeQuadrature(p,x1d,w1d);
    const int ng = w1d.n_elem;
    const int ng3d = pow(ng,3);

    x.resize(ng3d,3);
    w.resize(ng3d);
    for(int i=0,cnt=0; i<ng; i++){
      for(int j=0; j<ng; j++){
	for(int k=0; k<ng; k++,cnt++){
	  x(cnt,0) = x1d(k)*x1d(j);
	  x(cnt,1) = 1.0-x1d(j);
	  x(cnt,2) = (1.0-x1d(k))*x1d(j)*x1d(i);
	  w(cnt) = w1d(i)*w1d(j)*w1d(k)*(1.-x(cnt,1))*(1.-x(cnt,1)-x(cnt,0));
	}
      }
    }
  }
  //bary.save("bary",arma::raw_ascii);
  //x.save("gpts_tet",arma::raw_ascii);

  /*
  using namespace arma;
  mat x1d;
  vec w1d;

 
  LineQuadrature Line;
  Line.ComputeQuadrature(p,x1d,w1d);
  const int ng = w1d.n_elem;
  const int ng3d = pow(ng,3);
 

  std::cout << "quad n points: " << quad.npoints << std::endl;

  std::cout << "ng Tet: " << ng3d << std::endl;

  x.resize(ng3d,3);
  w.resize(ng3d);
  for(int i=0,cnt=0; i<ng; i++){
    for(int j=0; j<ng; j++){
      for(int k=0; k<ng; k++,cnt++){
	x(cnt,0) = x1d(k)*x1d(j);
	x(cnt,1) = 1.0-x1d(j);
	x(cnt,2) = (1.0-x1d(k))*x1d(j)*x1d(i);
	w(cnt) = w1d(i)*w1d(j)*w1d(k)*(1.-x(cnt,1))*(1.-x(cnt,1)-x(cnt,0));
      }
    }
  }
  */

    /*
    shards::CellTopology tet_4(shards::getCellTopologyData<shards::Tetrahedron<4> >() );

  DefaultCubatureFactory<double>  cubFactory;
  int cubDegree = p;
  Teuchos::RCP<Cubature<double> > tetCub = 
    cubFactory.create(tet_4, cubDegree); 

  
  int cubDim       = tetCub->getDimension();
  int numCubPoints = tetCub->getNumPoints();

  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);

  tetCub->getCubature(cubPoints, cubWeights);
  x.resize(numCubPoints,3);
  w.resize(numCubPoints);
  for(int i=0; i<numCubPoints; i++){
    for(int j=0; j<3; j++) x(i,j) = cubPoints(i,j);
    w(i) = cubWeights(i);

  }

  return 1;
    */

}

void HexQuadrature::ComputeQuadrature(int p, arma::mat &x, arma::vec &w) const{
  LineQuadrature Line;
  

  arma::mat gptsline;
  arma::vec gwline;

  Line.ComputeQuadrature(p,gptsline,gwline);
  const int ngl = gwline.n_elem;
  
  const int nghex = std::pow(ngl,3);
  x.resize(nghex,3);
  w.resize(nghex);

  for(int i = 0, cnt = 0; i < ngl; i++){
    for(int j = 0; j < ngl; j++){
      for(int k = 0; k < ngl; k++, cnt++){
	x(cnt,0) = gptsline(k);
	x(cnt,1) = gptsline(j);
	x(cnt,2) = gptsline(i);
	w(cnt) = gwline(i)*gwline(j)*gwline(k);
      }
    }
  }
}

void PrismQuadrature::ComputeQuadrature(int p, arma::mat &x, arma::vec &w)const{
  using namespace arma;

  LineQuadrature Line;
  TriQuadrature Tri;

  mat gptstri, gptsline;
  vec gwtri, gwline;
  Line.ComputeQuadrature(p,gptsline,gwline);
  Tri.ComputeQuadrature(p,gptstri,gwtri);
  
  const int ngtri = gwtri.n_elem;
  const int ngline = gwline.n_elem;
  const int ngprism = ngtri*ngline;


  x.resize(ngprism,3);
  w.resize(ngprism);
  for(int j=0,cnt=0; j<ngline; j++){
    for(int i=0; i<ngtri; i++,cnt++){
      x(cnt,0) = gptstri(i,0);
      x(cnt,1) = gptstri(i,1);
      x(cnt,2) = gptsline(j);
      w(cnt) = gwtri(i)*gwline(j);
    }
  } 

}


void PyramidQuadrature::ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const{
  std::cout << "Begeinning pyramid quadrature" << std::endl;
  LineQuadrature Line;
  QuadQuadrature Quad;
  arma::mat gptsquad, gptsline;
  arma::vec gwquad, gwline;
  
  Line.ComputeQuadrature(p,gptsline,gwline);
  Quad.ComputeQuadrature(p,gptsquad,gwquad);
  
  int ngquad = gwquad.n_elem;
  int ngline = gwline.n_elem;
  int ngpyramid = ngquad*ngline;

  x.resize(ngpyramid,3);
  w.resize(ngpyramid);
  
  for(int j = 0, cnt = 0; j < ngline; j++){
    double xyscale = 1.0-gptsline(j);
    for(int i = 0; i < ngquad; i++,cnt++){
      x(cnt,0) = xyscale*gptsquad(i,0);
      x(cnt,1) = xyscale*gptsquad(i,1);
      x(cnt,2) = gptsline(j);
      w(cnt) = gwline(j)*gwquad(i)*std::pow(xyscale,2);
    }
  }

}
