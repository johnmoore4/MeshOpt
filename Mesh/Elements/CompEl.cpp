#include "CompEl.h"
#include "ActiveMEl.h"
#include "ShapeFunctionMatrices.h"
#include <iomanip>

double CompEl::computeVolume() const {
  const arma::vec detJ = computeDetJacobianGP();

  const arma::vec& gw = (el.getShapeFunctionMatrices())->getQuadratureWeights();


  //cout << detJ << endl;

  return arma::sum(detJ%gw);

}

CompEl1D::CompEl1D(ActiveMEl el): CompEl(el){
  const arma::mat& lin_nodes = el.getLinearNodes();
  pseudoInverse = (lin_nodes.col(1)-lin_nodes.col(0)).t();
  pseudoInverse/=arma::norm(pseudoInverse,2);
}

const arma::vec CompEl1D::computeGaussNormals() const {
  std::cout << "computeGaussNormals not implemented for CompEl1D!" << std::endl;
}
const arma::vec CompEl1D::computeDetJacobianGP() const{
  std::cout << "computeDetJacobianGP not implemented for CompEl1D!" <<
    std::endl;
}
const arma::mat CompEl1D::computeJacobianPoint(const arma::vec3 grad[]) const{
  //std::cout << "in CompEl1D computeJacobianPoint" << std::endl;
  arma::mat J = pseudoInverse*grad[0];
  //std::cout << "J: " << J << std::endl;
  //std::cout << "grad[0]: " << grad[0] << std::endl;
  //std::cout << "pseudo inv: " << std::endl;
  //std::cout << pseudoInverse << std::endl;

  //std::cout << std::setprecision(15) << J << grad[0] <<  std::endl;
  return pseudoInverse*grad[0];
}

const arma::mat CompEl2D::computeJacobianPoint(const arma::vec3 grad[]) const{

  //std::cout << pseudoInverse << std::endl;
  //std::cout << grad[0] << std::endl;

  arma::mat22 J;
 
  J.unsafe_col(0) = pseudoInverse.t()*grad[0];
  J.unsafe_col(1) = pseudoInverse.t()*grad[1];
 
  
  return J;

}

const arma::mat::fixed<3,2> CompEl2D::computePseudoInverse() const{
  arma::mat::fixed<3,2> pseudoInv;
  const arma::mat nodes = el.getNodesMatrix();
  const std::vector<indtype>& cornern = el.getCornerNodeIndex();
 
  arma::vec3 v0 = nodes.unsafe_col(cornern[0]);
  arma::vec3 l = nodes.unsafe_col(cornern[1]) - v0;
  arma::vec3 r2 = nodes.unsafe_col(cornern[2]) - v0;
  arma::vec3 n = cross(l,r2);
  arma::vec3 m = cross(n,l);
  
  pseudoInv.col(0) = l/norm(l,2.0);
  pseudoInv.col(1) = m/norm(m,2.0);
  return pseudoInv;
}


const arma::vec CompEl::computeDetJacobianGP() const{
 

  const arma::mat nodes = el.getNodesMatrix();

  const arma::cube& dsf = 
    (el.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  
  arma::mat xyzderiv[3];

  for(int i=0; i<dsf.n_slices; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  const int ngp = dsf.n_cols;
  arma::vec detJ(ngp);

  arma::mat J;
  for(int i=0; i<ngp; i++){
    arma::vec3 derivpt[3];
    for(int j=0; j<dsf.n_slices; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(i);
    }
      J = computeJacobianPoint(derivpt);
      //cout << J << endl;
      detJ(i) = det(J);
      //if(detJ(i) < 0.0) cout << detJ(i) << endl;
  }
 
  return detJ;
}

const arma::mat CompEl3D::computeJacobianPoint(const arma::vec3 grad[]) const{

  arma::mat33 J;
  J.unsafe_col(0) = grad[0];
  J.unsafe_col(1) = grad[1];
  J.unsafe_col(2) = grad[2];
  return J;

}

/*
const arma::vec CompEl3D::computeDetJacobianGP() const{

  cout << "IN COMPUTE DET JACOBIAN GP 3D" << endl;
  const arma::mat nodes = el.getNodesMatrix();

  const arma::cube& dsf = 
    (el.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  
  const arma::mat dxi = nodes*dsf.slice(0);
  const arma::mat deta = nodes*dsf.slice(1);
  const arma::mat dzeta = nodes*dsf.slice(2);

  arma::vec detJ(dxi.n_cols);

  arma::mat33 J;
  for(int i=0; i<dxi.n_cols; i++){
    J.unsafe_col(0) = dxi.unsafe_col(i);
    J.unsafe_col(1) = deta.unsafe_col(i);
    J.unsafe_col(2) = dzeta.unsafe_col(i);
    detJ(i) = det(J);
    if(det(J) < 0.0) cout << "det(J) smaller than 0" << endl;
  }
  cout << "detJ n elem: " << detJ.n_elem << endl;
  cout << detJ << endl;

  return detJ;
}
*/
