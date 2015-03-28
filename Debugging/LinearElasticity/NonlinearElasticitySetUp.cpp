#include "NonlinearElasticitySolver.h"
#include "MeshContainer.h"

int NonlinearElasticitySolver::
SetUpProblem(const double h, const double lambda_l, const double mu_l){
  lambda = lambda_l;
  mu = mu_l;

  const node_map& nodes = mesh.getNodes();
  
  for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){
    arma::vec3 temp = nd->second->getXYZ();
    x0[nd->first] = {temp[0], temp[1]};
  }

  bcdef.clear();
  nodemap.clear();

  const element_set& subelements = mesh.getSubElements();
 
  arma::vec3 np = {0.0, 0.0, 1.0};
  for(auto el = subelements.begin(); el != subelements.end(); ++el){
    const int ncn = (*el)->numCornerNodes();
    const gind* cn = (*el)->getCornerNodes();
    arma::vec3 r = nodes.at(cn[1])->getXYZ() - nodes.at(cn[0])->getXYZ();
    arma::vec3 n = cross(r,np);
    arma::vec2 temp;
    temp(0) = n(0);
    temp(1) = n(1);
    for(int i=0; i<2; i++){
      auto it = bcdef.find(cn[i]);
      if(it != bcdef.end()){
	bcdef[cn[i]]+= temp;
      }
      else{
	bcdef[cn[i]] = temp;
      }
    }
  }
  for(auto bc = bcdef.begin(); bc != bcdef.end(); ++bc){
    bc->second/=-(norm(bc->second,2.0)/h);
  }
  
  int cnt = 0;
  for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){

    if(bcdef.find(nd->first) == bcdef.end()){
      nodemap[nd->first] = cnt;
      cnt++;
    }
  }
 
  const int ndof = nodes.size()*2;

  Teuchos::RCP<const Teuchos::Comm<int> > comm =     
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  Teuchos::RCP<node_type> node = rcp(new node_type());
 
 
  Teuchos::RCP<const map_type> data_map = 
    Tpetra::createContigMapWithNode<LO,GO>
    (ndof,ndof,comm,node);

  F = rcp(new multivector_type(data_map,1,true) );
  DU = rcp(new multivector_type(data_map,1,true) );

  KGLOB = rcp(new sparse_mat_type(data_map,0));


  return 1;

}
