#include "MeshContainer.h"
#include "Mel.h"
#include "El1D.h"
#include "El2D.h"
#include "El3D.h"


#include "Mnode.h"
//#include <utility>
#include <algorithm>

using namespace std;



arma::vec3 normalFromElnodes(const gind* elnodes, const node_map& nodes){

  arma::vec3 n0 = nodes.at(elnodes[0])->xyzvec3();
  arma::vec3 n1 = nodes.at(elnodes[1])->xyzvec3() - n0;
  arma::vec3 n2 = nodes.at(elnodes[2])->xyzvec3() - n0; 
  arma::vec3 n = cross(n1,n2);
  n/=arma::norm(n,2.0);
  return n;
}

element_set& MeshContainer::getElementsOfDim(int dim){
  if(dim == mesh_dim) return elements;
  else if(dim == mesh_dim-1) return subelements;
  else if(dim == mesh_dim-2) return subsubelements;
  else throw std::runtime_error("mesh dim request does not make sense!");
}
void MeshContainer::computeIsMeshPlanar(){
  // Determine if mesh is planar

  if(mesh_dim != 2){
    is_planar = false;
    return;
  }
  is_planar = false;
  return;
  /*
  const gind* el1nodes = (*elements.begin())->getCornerNodes();
  arma::vec3 n1 = normalFromElnodes(el1nodes,nodes);
  mesh_normal = n1;
  is_planar = true;
  cout << "h1" << endl;
  for(auto el=elements.begin(); el != elements.end(); ++el){
    const gind* elnodes = (*el)->getCornerNodes();
    arma::vec3 n = normalFromElnodes(elnodes,nodes);
    if(norm(n-n1,2.0) > 1.0e-8){
      is_planar = false;
      break;
    }
  }
  */
}

/*
const int MeshContainer::MeshDimension() const{
  
  int dim = 0;
  for(auto el=elements.begin(); el!= elements.end(); ++el){
    dim = std::max(dim,(*el)->getDim());
    //dim = std::max(dim,(*el).first->getDim());
  }
  return dim;
}
*/
