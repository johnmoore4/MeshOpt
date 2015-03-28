#include "ElementAnalyzer.h"
#include "CompEl.h"
#include "ShapeFunctionMatrices.h"
#include <iomanip>

ElementAnalyzer::ElementAnalyzer(std::shared_ptr<CompEl>& compel, 
				 const ShapeFunctionMatrices* sf_lin,
				 const arma::mat* ideal):
  compel(compel), sf_lin(sf_lin), ideal(ideal){
  
}

double ElementAnalyzer::computeDistortion(){
  using std::cout;
  using std::endl;

  const ActiveMEl& active = compel->getActiveElement();

  const arma::mat nodes = active.getNodesMatrix();


  arma::mat ideal_mat;
  if(!ideal) ideal_mat = active.getLinearNodes();
  else ideal_mat = *ideal;

  //const double eps = eps_fact*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::cube& dsf_ideal = sf_lin->getQuadratureSFDeriv();

  const arma::vec& gw = 
    (active.getShapeFunctionMatrices())->getQuadratureWeights();


  assert(dsf_ideal.slice(0).n_cols == dsf.slice(0).n_cols);

  const int dim = dsf.n_slices;

 
  arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal_mat*dsf_ideal.slice(i);
  }
  const int ngp = dsf.slice(0).n_cols;

  for(int d = 0; d < dim; d++){
    if(xyzderiv[d].n_cols != ngp) std::cout << "ncols/ngp: " << xyzderiv[d].n_cols << " " << ngp << std::endl;
    if(xyzderivI[d].n_cols != ngp) std::cout << "ncols/ngp Ideal: " << 
				     xyzderivI[d].n_cols << " " << ngp << std::endl;
  }

  /*
 if(dim == 1){
    arma::vec avnode = (nodes.col(0) + nodes.col(2))/2;
    double dist = arma::norm(nodes.col(2)-nodes.col(0),2);
    double diff_nodes = arma::norm(avnode-nodes.col(1),2);
    //std::cout << "node rel difference: " << diff_nodes/dist << std::endl;
    //std::cout << "ngp: " << ngp << std::endl;
    std::cout << xyzderiv[0] << std::endl;
    std::cout << xyzderivI[0] << std::endl;
  }
  */

  double element_area = 0.0;
  double eta_shap = 0.0;
  double deb = 0.0;

  //std::cout << "starting gp loop" << std::endl;
  if(dim == 1){
    arma::mat J, JI, S;
    arma::vec3 derivpt[3], derivptI[3];
    double area = 0.0;
    double dist_sum = 0.0;
    for(int i=0; i<ngp; i++){
      //std::cout << "h0" << std::endl;
      for(int j=0; j<dim; j++){
	derivpt[j] = xyzderiv[j].unsafe_col(i);
	derivptI[j] = xyzderivI[j].unsafe_col(i);
      }
      J = compel->computeJacobianPoint(derivpt);
      JI = compel->computeJacobianPoint(derivptI);

      double tc = arma::norm(derivpt[0],2);
      double ti = arma::norm(derivptI[0],2);
      double dj = det(J);
      area+= dj*gw(i);
      dist_sum+= dj*tc/ti*gw(i);
      
    }
    double Distortion = dist_sum/area;
    return Distortion;
  }

  arma::mat J, JI, S;
  arma::vec3 derivpt[3], derivptI[3];
  for(int i=0; i<ngp; i++){
    //std::cout << "h0" << std::endl;
    for(int j=0; j<dim; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(i);
      derivptI[j] = xyzderivI[j].unsafe_col(i);
    }
  
    J = compel->computeJacobianPoint(derivpt);
    JI = compel->computeJacobianPoint(derivptI);


    double dj = det(J);

    S = inv(JI)*J;

    double detS = det(S);
    double fronorm = norm(S,"fro");

    
    double eps = 0.0;

    double delta;
    if(detS < eps) delta = sqrt(eps*std::abs(eps-detS));
    else delta = 0.0;

    delta = 0.0;

    double sigma = detS;
    //double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));
  
    //double denom = double(dim)*pow(sigma,double(dim)/2.0);
    double denom = double(dim)*pow(sigma,2.0/double(dim));
    double eta = fronorm*(fronorm/denom);
    element_area+= std::abs(dj)*gw(i);
    eta_shap+= eta*eta*std::abs(dj)*gw(i);
    deb+= sigma;

  }
  
  double Distortion = sqrt(eta_shap/std::abs(element_area));


  distortion = distortion;


  return Distortion;
}
