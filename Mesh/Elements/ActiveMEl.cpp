#include "ActiveMEl.h"
#include <fstream>


using namespace std;
const arma::mat ActiveMEl::getNodesMatrix() const{
 
  //std::cout << "B" << std::endl;
  assert(getMeshElement());
  //std::cout << "ndof: " << Ndof() << std::endl;
  arma::mat allnodes(3,Ndof());
  const MEl* el = getMeshElement();
  assert(el);
  const gind* nodeind = el->getNodes();
  //std::cout << "h0 " << el->NumNodes() << std::endl;

  for(int i=0; i<Ndof(); i++){
    allnodes.unsafe_col(i) = globalnodes.at(nodeind[i])->xyzvec3();
  }
  //std::cout << "A" << std::endl;

  return allnodes;

  /*
  arma::mat allnodes(3,Ndof());
  // Corner nodes
  const std::vector<indtype>& cornerplace = getCornerNodeIndex();
  const gind* cornernodes = getCornerNodes();
  for(auto i=0; i<cornerplace.size(); i++){
    allnodes.unsafe_col(cornerplace[i]) = 
      globalnodes.at(cornernodes[i])->xyzvec3();
  }

  const MEl* el = getMeshElement();
 
  if(el->getOrder() == 1) return allnodes;

  // Child Nodes
 for(auto i=0; i<NumChildren(); i++){
   cout << "Getting child nodes!" << endl;
    const MEl* _child = el->getChild(i);

    ActiveMEl child(_child,index_factory,globalnodes,NULL);
    
    const std::vector<indtype>& interior = getInteriorNodes();

    if(interior.size() > 0){
      const int Nint = interior.size();

      const int child_orientation = el->getChildOrientation(i);

      const std::vector<indtype>& childplace = getChildNodes(i);

      const arma::mat childnodes = child.getNodesMatrix();

      const std::vector<indtype>& orient = 
	child.getOrientedNodes(std::abs(child_orientation));
  

      if(child_orientation < 0){
	const std::vector<indtype>& reversed = child.getReversedNodes();
      
	for(int j=0; j<Nint; j++){
	  allnodes.unsafe_col(childplace[interior[j]]) = 
	    childnodes.unsafe_col(orient[reversed[interior[j]]]);
	}  
      }
      else{
	for(int j=0; j<Nint; j++){
	  allnodes.unsafe_col(childplace[interior[j]]) = 
	    childnodes.unsafe_col(orient[interior[j]]);
	} 
      }
    }
 
  }

  // Interior nodes
  const std::vector<indtype>& interiorplace = getInteriorNodes();
  const std::unique_ptr<std::unique_ptr<MNode>[]>& interior_nodes = 
    el->getInteriorNodes();

  for(auto i=0; i<interiorplace.size(); i++){
    allnodes.unsafe_col(interiorplace[i]) = interior_nodes[i]->xyzvec3();
  }



  return allnodes;

  */
}

/*
const std::vector<gind> ActiveMEl::getAllNodes() const {
  std::vector<int> allnodes(Ndof());

  // Corner nodes
  const std::vector<indtype>& cornerplace = getCornerNodeIndex();
  const gind* cornernodes = getCornerNodes();
  for(auto i=0; i<cornerplace.size(); i++){
    allnodes[cornerplace[i]] = cornernodes[i];
  }

  if(el->getOrder() == 1) return allnodes;

  // Child Nodes
 for(auto i=0; i<NumChildren(); i++){
    const MEl* _child = el->getChild(i);

    ActiveMEl child(_child,index_factory,globalnodes,NULL);
    
    const std::vector<indtype>& interior = getInteriorNodes();

    if(interior.size() > 0){
      const int Nint = interior.size();

      const int child_orientation = el->getChildOrientation(i);

      const std::vector<indtype>& childplace = getChildNodes(i);

      const std::vector<indtype> childnodes = child.getAllNodes();
      const std::vector<indtype>& orient = 
	child.getOrientedNodes(std::abs(child_orientation));
  

      if(child_orientation < 0){
	const std::vector<indtype>& reversed = child.getReversedNodes();
      
	for(int j=0; j<Nint; j++){
	  allnodes[childplace[interior[j]]] = 
	    childnodes[orient[reversed[interior[j]]]];
	}  
      }
      else{
	for(int j=0; j<Nint; j++){
	  allnodes[childplace[interior[j]]] = 
	    childnodes[orient[interior[j]]];
	} 
      }
    }
 
  }

  // Interior nodes
  const std::vector<indtype>& interiorplace = getInteriorNodes();
  const std::unique_ptr<std::unique_ptr<MNode>[]>& interior_nodes = 
    el->getInteriorNodes();

  for(auto i=0; i<interiorplace.size(); i++){
    allnodes[interiorplace[i]] = interior_nodes[i]->getND();
  }



  return allnodes;

}
*/

void ActiveMEl::writeGMSH(ofstream &out, int num){

  
    out << num+1 << " " <<  ndind->getGMSHType() << 
      " " << 2 << " " << el->getBCTag() << " " << el->getGeoEntity()+1;
  


    const std::vector<indtype>& gmshindex = ndind->getGMSHIndex();


    const MEl* el = getMeshElement();
    const gind* nodes = el->getNodes();
    
    //const std::vector<gind> nodes = getAllNodes();
 
    //assert(gmshindex.size() == nodes.size());

    for(int i=0; i<Ndof(); i++){

      out << " " << nodes[gmshindex[i]]+1;
    }

    out << endl;

  
}

arma::mat ActiveMEl::getLinearNodes() const{
  const std::vector<indtype>& cnodes = ndind->getCornerNodes();
  const MEl* el = getMeshElement();
  const gind* nodeind = el->getNodes();

  const int ncn = cnodes.size();
  arma::mat linear(3,ncn);
  for(int i = 0; i < ncn; i++){
    linear.unsafe_col(i) = globalnodes.at(nodeind[cnodes[i]])->xyzvec3();
  }
  return linear;
}
