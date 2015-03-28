#include "El1D.h"
#include "Mnode.h"
//#include "shapeFunctionHandler.h"
#include "myEdge.h"
#include "myFace.h"



//void El1D::swapInParent(element_ptr newel, signed char orientation){
  //parent->setChild(parent_edge,newel,orientation);
//}

/*
arma::mat El1D::nodeParamsOnEdge() const {
  
  const int nn = _nodes.size();
  arma::mat s(1,nn);

  for(int i=0; i<nn; i++){
    if(_nodes[i]->dim() == 0){
      MVertNode* temp = static_cast<MVertNode*>(_nodes[i].get());
      s(0,i) = _edge->paramsOnEdge(temp->getVertex()); 
      // Fix periodic case
      if(i > 0 && std::abs(s(0,i)-_edge->umin) < 1.0e-10) s(0,i) = _edge->umax;
    }
    else if(_nodes[i]->dim() == 1){
      MEdgeNode* temp = static_cast<MEdgeNode*>(_nodes[i].get());

      s(0,i) = temp->getParam();
    }
    else assert(1);
  }

  if(!is_finite(s)){
    cout << "s is not finite in node params on edge!" << endl;
    assert(1);
  }
  

  return s;
  
}
*/

/*
arma::vec3 El1D::elementNormal(arma::vec3& proj) const{
  
  arma::mat nds = nodes();
  arma::mat xyz_xi = nds*_sf->sfderiv.slice(0);
  xyz_xi = xyz_xi.t();
  arma::vec I = 1.0/sqrt(pow(xyz_xi.col(0),2.0) + pow(xyz_xi.col(1),2.0));
  arma::vec3 normal;
  normal(0) = sum(I%xyz_xi.col(1));
  normal(1) = -sum(I%xyz_xi.col(0));
  normal(2) = 0.0;
  normal/= norm(normal,2.0);

  return normal;
  
}
*/

std::vector<node_ptr> projectPointsOnEdge(const arma::mat& xyz0, 
					  const arma::mat& s0,
					  const myEdge* edge){
  /*
  const int nn = xyz0.n_cols;
  std::vector<node_ptr> nodes(nn);

  for(int i=0; i<nn; i++){
    double s = edge->edgeParamsFromPoint(xyz0.colptr(i),s0(0,i));
    arma::vec3 xyz = edge->param2xyz(s);
    if(!is_finite(xyz)){
      cout << "xyz is not finite in projectPointsOnEdge!" << endl;
      assert(1);
    }
    nodes[i] = std::make_shared<MEdgeNode>
      (MEdgeNode(0,xyz.memptr(),1,1,s,edge));
  }

  return nodes;
  */
}

/*
std::vector<node_ptr> El1D::projectPhysicalNodes(){
  
  arma::mat nodes0 = nodes();
  arma::mat map_matrix = _sf->sfmap.cols(_sf->interior_nodes);
  arma::mat xyz_new = nodes0*map_matrix;


  if(_edge){
 
   // Project on edge
    arma::mat s = nodeParamsOnEdge();
    arma::mat s_new = s*map_matrix;  
    if(_porder > 1) _curved = true;
    return projectPointsOnEdge(xyz_new,s_new,_edge);
  }
  else if(_face){
    if(_porder > 1) _curved = true;
    arma::mat uv = nodeParamsOnFace(_nodes,_face);
    arma::mat uv_new = uv*map_matrix;
    return projectPointsOnFace(xyz_new,uv_new,_face);
  }
  else{
    return createFreeNodes(xyz_new);
 
  }
  
}
*/

/*
const arma::mat El1D::computeJacobian(arma::vec3* grad) const{
  arma::mat::fixed<1,1> J;
  J[0] = sqrt(dot(grad[0],grad[0]));

  return J;
}

const arma::vec El1D::computeGradDetJacobian(arma::vec3* grad) const{
  
  return grad[0]/computeJacobian(grad);
}

const arma::mat El1D::computeHessDetJacobian(arma::vec3* grad) const{
  arma::mat33 Hess;
  arma::vec3 gr = grad[0];
  double detJ = sqrt(dot(gr,gr));
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      Hess(i,j) = -0.5*gr[i]*gr[j]/pow(detJ,3);
      if(i==j) Hess(i,j)+= 1.0/detJ;
    }
  }
  return Hess;
}

const arma::vec El1D::computeGradFroNorm(arma::vec3* grad) const{
  return computeGradDetJacobian(grad);
}

const arma::mat El1D::computeHessFroNorm(arma::vec3* grad) const{
  return computeHessDetJacobian(grad);
}

*/
