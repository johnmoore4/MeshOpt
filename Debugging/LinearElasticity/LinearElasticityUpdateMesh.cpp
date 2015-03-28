#include "LinearElasticityEvaluator.h"
#include "MeshContainer.h"

int LinearElasticityEvaluator::UpdateMesh(){
  std::cout << "in update mesh " << std::endl;
  /*
 Teuchos::ArrayRCP<ST> DU_data = DU->get1dViewNonConst();
 

 node_map& nodes = mesh.getNodesNC();

 
 const int ndof = nodes.size();
 for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){
   arma::vec3 du;
   du.zeros();
   du(0) = DU_data[nd->first];
   du(1) = DU_data[nd->first+ndof];
   arma::vec3 xnew = nd->second->getXYZ() + du;
   nd->second->setXYZ(xnew.memptr());
 }
  */
  return 1;
} 
