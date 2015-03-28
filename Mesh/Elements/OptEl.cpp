//#define ARMA_NO_DEBUG
#include "OptEl.h"
#include "ShapeFunctionMatrices.h"
#include "GeometryContainer.h"


#define eps_fact 1.0e12
const double OptEl::computeGradMeritParam(arma::mat&gradMerit, 
					  double& distortion,
					  double factor,
					  double minDetS) const{

  //std::cout << "compute grad merit param" << std::endl;

  arma::mat grad_temp;
  double detS;
  double mer = computeGradMerit(grad_temp,distortion,detS,factor,minDetS);

  ActiveMEl& active = compel.getActiveElement();
  const MEl* el = active.getMeshElement();
  const gind* nind = el->getNodes();


  const node_map& nodes = active.getGlobalNodes();
  const int nn = active.Ndof();
  
  //const int ncn = el->numCornerNodes();
  //const gind* cn = el->getCornerNodes();

  gradMerit.resize(3*nn,1);
  gradMerit.zeros();
  
  for(int i=0, dofcnt=0; i<nn; i++){
    const MNode* nd = nodes.at(nind[i]).get();
    const int dim = nd->getType();
    const short int entity = nd->getGeoEntity();
 
    if(nd->getType() == 1){
      std::vector<myEdge>& edges = geometry->getEdgesNC();
      myEdge& edge = edges[entity];
      
      const double* u = nd->getParametricCoords();
      arma::rowvec3 deriv;
      edge.D1(u,deriv);
      const double uderiv = arma::as_scalar(deriv*grad_temp.unsafe_col(i));
      gradMerit(dofcnt) = uderiv;
    }
    else if(nd->getType() == 2){
      std::vector<myFace>& faces = geometry->getFacesNC();
      myFace& face = faces[entity];
      if(1){
	//if(ace.areParamsvalid(){

	const double* uv = nd->getParametricCoords();
	arma::mat::fixed<2,3> deriv;
	face.D1(uv,deriv);
	arma::vec2 uvderiv = deriv*grad_temp.unsafe_col(i);
	for(int j=0; j<2; j++){
	  gradMerit(dofcnt+j) = uvderiv[j];
	}
	//cout << "face: " << entity << endl;
      }
    }
    else if(nd->getType() == 3){
      for(int j=0; j<3; j++){
	gradMerit(dofcnt+j) = grad_temp.at(j,i);
      }
    }

    dofcnt+= dim;
    
    //if(i == ncn-1) gradMerit.resize(dofcnt);
  }
  
  //gradMerit.resize(dofcnt);
  return mer;
}


const double OptEl::computeMinDetS() const{

  if(ideal.is_empty()) return 1.0;

  const ActiveMEl& active = compel.getActiveElement();

  const arma::mat nodes = active.getNodesMatrix();

  const double eps = eps_fact*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::vec& gw = 
      (active.getShapeFunctionMatrices())->getQuadratureWeights();


 
  const int dim = dsf.n_slices;

  //cout << "dim: " << dim << endl;

  arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal*dsf.slice(i);
  }
  const int ngp = dsf.n_cols;

  double element_area = 0.0;
  double eta_shap = 0.0;
  double deb = 0.0;
  double minDetS = 1.0e15;

  arma::mat J, JI, S;
    arma::vec3 derivpt[3], derivptI[3];
  for(int i=0; i<ngp; i++){

    for(int j=0; j<dim; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(i);
      derivptI[j] = xyzderivI[j].unsafe_col(i);
    }
    J = compel.computeJacobianPoint(derivpt);
    JI = compel.computeJacobianPoint(derivptI);
 
    //double dj = det(J);
 
    S = inv(JI)*J;

    minDetS = std::min(minDetS,det(S));

  }
 
  return minDetS;

};

double OptEl::computeDistortion() const{


  using std::cout;
  using std::endl;

  if(ideal.is_empty()) return 1.0;


  //cout << "arma mat prealloc: " << ARMA_MAT_PREALLOC << endl;
 const ActiveMEl& active = compel.getActiveElement();

  const arma::mat nodes = active.getNodesMatrix();

  const double eps = eps_fact*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::cube& dsf_ideal = sf_ideal->getQuadratureSFDeriv();

  const arma::vec& gw = 
      (active.getShapeFunctionMatrices())->getQuadratureWeights();


 

  const int dim = dsf.n_slices;

  //cout << "dim: " << dim << endl;

  arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal*dsf_ideal.slice(i);
  }
  const int ngp = dsf.n_cols;

  double element_area = 0.0;
  double eta_shap = 0.0;
  double deb = 0.0;

  arma::mat J, JI, S;
    arma::vec3 derivpt[3], derivptI[3];
  for(int i=0; i<ngp; i++){

    for(int j=0; j<dim; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(i);
      derivptI[j] = xyzderivI[j].unsafe_col(i);
    }
    J = compel.computeJacobianPoint(derivpt);
    JI = compel.computeJacobianPoint(derivptI);
 
    double dj = det(J);
 
   

    S = inv(JI)*J;

    double detS = det(S);
    double fronorm = norm(S,"fro");

    double delta;
    //if(detS < eps) delta = eps;
    if(detS < eps) delta = sqrt(eps*std::abs(eps-detS));
    else delta = 0.0;
    //delta = eps;
    //delta = 0.0;

    double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));
    //double sigma = 
    //  0.5*(detS + sqrt(detS*detS + dim*dim*pow(delta,4.0/(double)dim)));
    //double sigma = 0.5*(detS + std::abs(detS) + delta);
    //double sigma = 0.5*(detS + sqrt(detS*detS) + delta);
    //double sigma = 0.5*(detS + sqrt(detS*detS) + delta);
    
  

 
    double denom = double(dim)*pow(sigma,2.0/double(dim));
    double eta = fronorm*(fronorm/denom);
    element_area+= std::abs(dj)*gw(i);
    eta_shap+= eta*eta*std::abs(dj)*gw(i);
    deb+= sigma;

  }
  //if(eta_shap!= eta_shap) cout << "Eta shap is NaN!" << endl;
  double distortion = sqrt(eta_shap/std::abs(element_area));

  
  /*
  if(distortion != distortion){
    cout << "distortion is NaN!" << endl;
    cout << eta_shap << " " << element_area << " " << deb <<  endl;
    //cout << "dj: " << dj << endl;
    cout << nodes << endl;
    cout << ideal << endl;
    cout << S << endl;
  }
  
  
  test(0) = distortion;
  if(!test.is_finite()){
    cout << "Not finite distortion: " << distortion << endl;
    cout << "eta_shap: " << eta_shap << " " << element_area << endl;
    cout << deb << endl;
  }
  */
  //if(distortion > 1e6) cout << "distortion: " << distortion << endl;
  //if(distortion < 1.0-1.0e-10){
  //  cout << "distortion: " << distortion << endl;
  //  cout << eta_shap << " " << element_area << endl;
  //}

  return distortion;

};

const double OptEl::computeMerit() const{
  return pow(computeDistortion()-1.0,2.0);
  //return (computeDistortion() - 1.0);
}



const double OptEl::computeGradMerit(arma::mat& gradMerit,
				     double& dist,
				     double& DetS,
				     double factor,
				     double minDetS) const{

  //std::cout << "begin computeGradMerit" << std::endl;

  const ActiveMEl& active = compel.getActiveElement();

  const arma::mat nodes = active.getNodesMatrix();

  const double eps = eps_fact*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::cube& dsf_ideal = sf_ideal->getQuadratureSFDeriv();

  const arma::vec& gw = 
      (active.getShapeFunctionMatrices())->getQuadratureWeights();


  const MEl* el = active.getMeshElement();

  const int dim = dsf.n_slices;
  const int ngp = dsf.n_cols;
  const int ndof = dsf.n_rows;

  //std::cout << "element type: " << el->getElementType() << std::endl;

  if(el->getElementType() == 4){
    //cout << "tet gp: " << ngp << endl;
  }
  else if(el->getElementType() == 5){
    //cout << "quad gp: " << ngp << endl;
  }
  
  //double xyzderiv_mem[3*ngp][3];
  double xyzderiv_mem[3][3*ngp];

  arma::mat xyzderiv[3] = {arma::mat(xyzderiv_mem[0],3,ngp,false),
			   arma::mat(xyzderiv_mem[1],3,ngp,false),
			   arma::mat(xyzderiv_mem[2],3,ngp,false)};
 
  
  //arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 

  
  double xyzderivI_mem[3][3*ngp];
  arma::mat xyzderivI[3] = {arma::mat(xyzderivI_mem[0],3,ngp,false),
			   arma::mat(xyzderivI_mem[1],3,ngp,false),
			   arma::mat(xyzderivI_mem[2],3,ngp,false)};
  

  //arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal*dsf_ideal.slice(i);
  }
  
  //std::cout << "here" << std::endl;

  assert(xyzderivI[0].n_cols == xyzderiv[0].n_cols);

  //std::cout << ideal << std::endl;

  double element_area = 0.0;
  double eta_shap = 0.0;
  double eta_actual = 0.0;
  double ideal_area = 0.0;
  //arma::mat33 J, JI, invJI, S;
  
  
  double Jmem[dim*dim], JImem[dim*dim], invJImem[dim*dim], Smem[dim*dim];
  arma::mat J(Jmem,dim,dim,false);
  arma::mat JI(JImem,dim,dim,false);
  arma::mat invJI(invJImem,dim,dim,false);
  arma::mat S(Smem,dim,dim,false);
  
  double gradDetJ_mem[3*dim], gradEtaShape_mem[3*dim], gradsf_mem[dim*ndof];
  double gradFroNorm_mem[3*dim], gradDetS_mem[3*dim];
  arma::mat gradDetS(gradDetS_mem,3,dim,false);
  arma::mat gradDetJ(gradDetJ_mem,3,dim,false);
  arma::mat gradEtaShape(gradEtaShape_mem,3,dim,false);
  arma::mat gradFroNorm(gradFroNorm_mem,3,dim,false);
  arma::mat gradsf(gradsf_mem,dim,ndof,false);

  double ESA_mem[3*ndof], AA_mem[3*ndof];
  
  arma::mat gradEtaShapeAccum(ESA_mem,3,ndof,false);
  arma::mat gradAreaAccum(AA_mem,3,ndof,false);


  gradEtaShapeAccum.zeros();
  gradAreaAccum.zeros();

  //arma::mat gradEtaShapeAccum(3,ndof,arma::fill::zeros);
  //arma::mat gradAreaAccum(3,ndof,arma::fill::zeros);
 
  

  double delta, delta_detS, stab, stab_detS;
  
  //double deriv[dim][3];
  //double derivI[dim][3];
  const double *deriv[dim], *derivI[dim];
  /*
  double min_detS = 0.0;
  for(int i=0; i<ngp; i++){
 
    for(int j=0; j<dim; j++){
      derivI[j] = xyzderivI[j].colptr(i);
    }
    computeJacobianPoint(derivI,JI);

    min_detS = std::min(min_detS,det(JI));
  }
  */

  //std::cout << "before GP loop" << std::endl;

  for(int i=0; i<ngp; i++){
 
    for(int j=0; j<dim; j++){
      deriv[j] = xyzderiv[j].colptr(i);
      derivI[j] = xyzderivI[j].colptr(i);
      //J.unsafe_col(j) = xyzderiv[j].unsafe_col(i);
      //JI.unsafe_col(j) = xyzderivI[j].unsafe_col(i);
    }
    computeJacobianPoint(deriv,J);
    computeJacobianPoint(derivI,JI);

    
    const double dj = det(J);
    const double absDJ = std::abs(dj);


    invJI = inv(JI);

    S = invJI*J;
    //S = pseudoInv*invJI*J;
    
 

    const double detS = det(S);
    const double fronorm = norm(S,"fro");


    //DetS = detS;
    

    computeGradDetS(J,gradDetJ,1.0);
    //computeGradDetS(J,gradDetS, detS);
    computeGradFrobNorm(S,invJI,gradFroNorm);

    //gradDetJ = dj*inv(J).t();

    //gradDetS = detS*inv(J).t();

    //computeGradDetS(J,invJI, gradDetJ);

    //computeGradDetJ(J,gradDetJ);

    //computeGradFroNorm(invJI, S, gradFroNorm);

    //const double detS_gradDetJ = det(invJI);
    //const double alpha = 1.0e-5;
    
    double alpha = 1.0e-4;
    if(detS < 0.0){

      //delta = min_detS*sqrt(alpha*alpha+alpha);
      stab_detS = 0.0;
      // cout << "minDetS: " << minDetS << endl;
      delta = std::abs(detS)*sqrt(alpha*alpha+alpha);
      //delta = std::abs(arma::det(arma::inv(JI)))*sqrt(alpha*alpha+alpha);

      //delta = detS*sqrt(alpha*alpha+alpha);
      stab =  sqrt(detS*detS + 4*delta*delta);
      //stab_detS = 2*detS/sqrt(detS*detS + delta*delta);
      //stab_detS = (2*detS*(1+4*(alpha*alpha+alpha)))/sqrt(detS*detS + delta*delta);
    }
    else{
      stab = detS;
      stab_detS = 1;
    }
    //stab_detS = detS/stab;
    //stab_detS = 0.0;

    /*
    if(detS < eps){
      delta = sqrt(eps*(eps-minDetS));
      stab =  sqrt(detS*detS + dim*dim*pow(delta,4.0/(double)dim));
      stab_detS = detS/stab;
    }
    else{
      stab = detS;
      stab_detS = 1.0;
    }
    */
    
    /*
    if(detS < eps && DetS != -1){
      delta = sqrt(eps*(eps-detS));
      const double delta_detS = -eps/(2.0*delta);
      stab =  sqrt(detS*detS + dim*dim*pow(delta,4.0/(double)dim));
      stab_detS = (detS + 2*dim*pow(delta,4.0/dim-1)*delta_detS)/stab;

      //stab_detS = detS/stab;
    }
    else{
      stab = detS;
      stab_detS = 1.0;
    }
    */
    /*
    double alpha = 1.0e-3;
    if(detS < eps && DetS != -1){
      delta = detS*sqrt(alpha*alpha+alpha);
      const double delta_detS = sqrt(alpha*alpha+alpha);
      //const double delta_detS = -eps/(2.0*delta);
      
      stab =  sqrt(detS*detS + dim*dim*pow(delta,4.0/(double)dim));
      stab_detS = (detS + 2*dim*pow(delta,4.0/dim-1)*delta_detS)/stab;

      //stab_detS = detS/stab;
    }
    else{
      stab = detS;
      stab_detS = 1.0;
    }
    */


    const double sigma = 0.5*(detS + stab);

    const double sigma_detS = 0.5*(1.0+stab_detS);
    



    double denom, denom_sigma;
    //double eta;
    if(dim == 3){
      const double cb = cbrt(sigma);
      denom = dim*cb*cb;
      denom_sigma = 2.0/cb;
      //const double temp = pow(dim,double(dim)/2.0);
      //denom = temp*sigma;
      //denom_sigma = temp;
      //eta = pow(fronorm,dim)/denom;

    }
    else if(dim == 2){
      denom = dim*sigma;
      denom_sigma = dim;
    }
   
    const double eta = fronorm*fronorm/denom;

    const double eta2 = eta*eta;
    

    const double FroOverDenom = fronorm/denom;
    eta_shap+= eta2*absDJ*gw.at(i);
    element_area+= absDJ*gw.at(i);
    ideal_area+= std::abs(det(JI))*gw.at(i);

    //eta_actual+= eta2act*absDJ*gw.at(i);
    
    gradEtaShape= (-2*eta*absDJ*FroOverDenom*FroOverDenom*denom_sigma*
		   sigma_detS*detS + dj*mySign(dj)*eta2)*gradDetJ;

    //gradEtaShape+= mySign(dj)*eta2*gradDetJ;

    gradEtaShape+= eta*absDJ*(4.0/denom)*gradFroNorm;
    

    /*
    gradEtaShape= (-2*eta*absDJ*FroOverDenom*FroOverDenom*denom_sigma*sigma_detS*detS_gradDetJ + mySign(dj)*eta2)*gradDetJ;
    gradEtaShape+= eta*absDJ*(4.0/denom)*(invJI.t()*S);
    */
    
    
    for(int j=0; j<dim; j++){
      for(int k=0; k<ndof; k++) gradsf.at(j,k) = dsf.at(k,i,j)*gw.at(i);
    }



    gradEtaShapeAccum+= gradEtaShape*gradsf;
    gradAreaAccum+= mySign(dj)*dj*gradDetJ*gradsf; 
 

    
  }
 

  //std::cout << "After GP loop" << std::endl;

  double distortion = sqrt(eta_shap/std::abs(element_area));
  arma::mat grad_distortion = (gradEtaShapeAccum/std::abs(element_area) -
			       eta_shap/pow(element_area,2.0)*
			       mySign(element_area)*gradAreaAccum)/(2.0*distortion);
  bool is_nan = false;
  dist = distortion;
  //std::cout << dist << std::endl;
  if(dist == 1.0/0.0){
    std::cout << "is INF" << std::endl;
    std::cout << nodes << std::endl;
  }
  if(distortion != distortion){

    for(int i=0; i<grad_distortion.n_elem; i++){
      if(grad_distortion[i] != grad_distortion[i]) is_nan = true;
      break;
    }
    if(is_nan){
      //cout << grad_distortion << endl;
    }
    //cout << "distortion: " << distortion << DetS <<  endl;
    //cout << nodes << endl;
  }

  //double factor;
 
  if(el->getElementType() != 6 && el->getElementType() != 3){
    factor = 1.0;
  }
  else{
    //cout << "factor " << factor << endl;
  }
  //cout << "factor: " << factor << endl;

  //gradMerit = 2.0*factor*element_area*(distortion-1.0)*grad_distortion +
  //factor*pow(distortion-1.0,2.0)*gradAreaAccum;
  //return element_area*factor*pow(distortion-1.0,2.0);
  
  //gradMerit.zeros(3,ndof);

  //gradMerit = grad_distortion;
  //return distortion;

  //gradMerit = factor*grad_distortion;
  //return factor*(distortion-1);

  //gradMerit = factor*gradEtaShapeAccum/(2.0*sqrt(eta_shap));
  //return factor*sqrt(eta_shap);
  


  //gradMerit = sqrt(element_area)*factor*grad_distortion + 
  //  factor*(distortion-1)/(2.0*sqrt(element_area))*gradAreaAccum;
  //return sqrt(element_area)*factor*(distortion-1);

  //gradMerit = sqrt(ideal_area)*grad_distortion;
  //return sqrt(ideal_area)*(distortion - 1.0);
 
  //gradMerit = 2.0*(factor*(ideal_area*(distortion-1.0)))*factor*grad_distortion;
  // return factor*pow(ideal_area*(distortion-1.0),2);
  
  gradMerit = 2.0*(factor*(distortion-1.0))*factor*grad_distortion;
  const double fdm1 = factor*(distortion-1.0);
  //if(is_nan == true) return 1.0e50;

  //std::cout << "end computeGradMerit" << std::endl;

  return fdm1*fdm1;
  
  //return pow(factor*(distortion-1.0),2);




  //return pow(factor*distortion-1.0,2);
  //return 1.0;

  //return (factor*distortion-1.0)*(factor*distortion-1.0);

 

  gradMerit = 2.0*factor*element_area*(distortion-1.0)*grad_distortion +
      factor*pow(distortion-1.0,2.0)*gradAreaAccum;
  return element_area*factor*pow(distortion-1.0,2.0);
  

}


const arma::mat OptEl::computeGradMeritFD() const{
  
  const ActiveMEl& active = compel.getActiveElement();
  const node_map& nodes = active.getGlobalNodes();

  double merit0 = computeMerit();
  const MEl* el = active.getMeshElement();
  const int nn = el->NumNodes();
  const gind* elnodes = el->getNodes();
  double step = 1.0e-8;
  arma::mat gradmerit(nn,3);
  for(int i=0; i<nn; i++){
    
    arma::vec3 xyz = nodes.at(elnodes[i])->getXYZ();

    for(int j=0; j<3; j++){
 
      arma::vec3 xyzfd = xyz;
      xyzfd(j)+= step;
      nodes.at(elnodes[i])->setXYZ(xyzfd.memptr());
      double meritFD = computeMerit();
      gradmerit(i,j) = (meritFD-merit0)/step;
      nodes.at(elnodes[i])->setXYZ(xyz.memptr());
    }
    
    
  }
  return gradmerit.t();
  
}
void OptEl2D::computeJacobianPoint(const double* deriv[], arma::mat& J) const{
  arma::mat::fixed<3,2> Jtemp;
  for(int i=0; i<2; i++){
    for(int j=0; j<3; j++){
      Jtemp.at(j,i) = deriv[i][j];
    }
  }
  const arma::mat::fixed<3,2>& pseudoInverse = compel2d.getPseudoInverse();
  J = pseudoInverse.t()*Jtemp;
}

void OptEl2D::computeGradDetS(const arma::mat& J, 
			      arma::mat& grad,
			      const double coeff) const{
  const arma::mat::fixed<3,2>& pseudoInv = compel2d.getPseudoInverse();
  grad = coeff*pseudoInv*inv(J).t();

 
}

void OptEl2D::computeGradFrobNorm(const arma::mat& S,
				  const arma::mat& invJI,
				  arma::mat& grad) const{
  
  const arma::mat::fixed<3,2>& pseudoInv = compel2d.getPseudoInverse(); 
  grad = pseudoInv*invJI.t()*S;
}

void OptEl3D::computeJacobianPoint(const double* deriv[], arma::mat& J) const{
  
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      J.at(i,j) = deriv[j][i];
    }
  }
  
}

void OptEl3D::computeGradDetS(const arma::mat& J,
			      arma::mat& grad,
			      const double coeff) const{

  grad = coeff*inv(J).t();

//grad = det(J)*inv(J).t();

}


void OptEl3D::computeGradFrobNorm(const arma::mat& S,
				  const arma::mat& invJI,
				  arma::mat& grad) const{

  grad = invJI.t()*S;
}


const double OptEl::computeDebug() const{}
/*
const double OptEl::computeDebug() const{
  const arma::mat nodes = active.getNodesMatrix();

  const double eps = 2.0*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::vec& gw = 
      (active.getShapeFunctionMatrices())->getQuadratureWeights();


 
  const int dim = dsf.n_slices;


  arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal*dsf.slice(i);
  }
  const int ngp = dsf.n_cols;
  const int ndof = dsf.n_rows;
  double element_area = 0.0;
  double eta_shap = 0.0;
  arma::mat J, JI, invJI, S;
 
  double deb = 0.0;


  for(int i=0; i<ngp; i++){
    arma::vec3 derivpt[3], derivptI[3];
    for(int j=0; j<dim; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(i);
      derivptI[j] = xyzderivI[j].unsafe_col(i);
    }
    J = compel->computeJacobianPoint(derivpt);
    JI = compel->computeJacobianPoint(derivptI);
 
    double dj = std::abs(det(J));
    element_area+= dj*gw(i);
    invJI = inv(JI);
    S = invJI*J;
    double detS = det(S);
    double delta;
    if(detS < eps) delta = sqrt(eps*(eps-detS));
    else delta = 0.0;
    
    double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));

    double fronorm = norm(S,"fro");

    double denom = double(dim)*pow(sigma,double(dim)/double(3));
    

   double eta = fronorm*fronorm/denom;

   eta_shap+= eta*eta*dj*gw(i);
    
  }

  double distortion = sqrt(eta_shap/std::abs(element_area));
  //deb = element_area;
  deb = pow(distortion-1.,2);

  return deb;

}
*/


const arma::mat OptEl::debugGradMeritFD() const{}
/*
const arma::mat OptEl::debugGradMeritFD() const{
  node_map& nodes = active.getGlobalNodes();

  double merit0 = computeMerit();
  const MEl* el = active.getMeshElement();
  const int ncn = el->numCornerNodes();
  const gind* cn = el->getCornerNodes();
  double step = 1.0e-8;
  arma::mat graddebug(ncn,3);
  for(int i=0; i<ncn; i++){
    
    arma::vec3 xyz = nodes[cn[i]]->getXYZ();

    for(int j=0; j<3; j++){
 
      arma::vec3 xyzfd = xyz;
      xyzfd(j)+= step;
      nodes[cn[i]]->setXYZ(xyzfd.memptr());
      double meritFD = computeMerit();
      graddebug(i,j) = (meritFD-merit0)/step;
      nodes[cn[i]]->setXYZ(xyz.memptr());
    }
    
    
  }
  return graddebug;

}
*/
