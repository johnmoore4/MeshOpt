#include "El2D.h"
#include "Mnode.h"
//#include "shapeFunctionHandler.h"
#include "myEdge.h"
#include "myFace.h"

using namespace std;

//void El2D::swapInParent(element_ptr newel, signed char orientation){
  // parent->setChild(parent_face,newel,orientation);
//}


//void El2D::setChild(unsigned char ind, element_ptr el, 
//		    signed char orientation){
  //element1D_ptr el1d = static_pointer_cast<El1D>(el);
  //_children[ind] = el1d;
  //_child_orientation[ind] = orientation;
		      //}

/*
std::vector<node_ptr> El2D::projectPhysicalNodes(){
  
  const arma::uvec& exterior_nodes = _sf->getExteriorNodes();

  arma::mat nodes0 = nodes();
  //nodes0 = nodes0.cols(exterior_nodes);
  arma::mat map_matrix = _sf->sfmap.cols(_sf->interior_nodes);
  //map_matrix = map_matrix.rows(exterior_nodes);

  arma::mat xyz_new = nodes0*map_matrix;

  
  if(_face){
 
    arma::mat uv = nodeParamsOnFace(_nodes,_face);
    //uv = uv.cols(exterior_nodes);
    arma::mat uv_new = uv*map_matrix;

    if(!is_finite(uv) || !is_finite(uv_new)){
      cout << "uv is not finite in El2D project Physical nodes" << endl;
      cout << uv << endl;
      cout << uv_new << endl;
    }
    if(_porder > 1) _curved = true;

    return projectPointsOnFace(xyz_new,uv_new,_face);

  }
  else{
    return createFreeNodes(xyz_new);
  }
  
  
}
*/

//const bool El2D::is_curved() const{
 
  /*
  for(int i=0; i<_children.size(); i++){
    if(_children[i]->is_curved()) return true;
  }
  return _curved;
  */

//}

/*
void El2D::getChildNodes(){
  //for(int i=0; i<_children.size(); i++){
    //const std::vector<node_ptr>& child_nodes = _children[i]->nodes_vector();
    //const std::vector<node_ptr>& int_child_nodes = 
    //  _children[i]->interior_nodes();
    
  // }
}
*/

/*
arma::vec3 El2D::elementNormal() const{
  
  arma::mat nds = nodes();
  arma::mat xyz_xi = nds*_sf->sfderiv.slice(0);
  arma::mat xyz_et = nds*_sf->sfderiv.slice(1);
  const int ng = xyz_xi.n_cols;
  arma::vec3 normal;
  normal.zeros();
  for(int i=0; i<ng; i++){
    normal+= cross(xyz_xi.col(i),xyz_et.col(i));
  }
  normal/=norm(normal,2.0);
  return normal;
  
}
*/

