#include "BLMesh.h"
#include "MeshContainer.h"

BLMesh::BLMesh(const MeshContainer& mesh){
  const element_set& elements = mesh.getElements();
  

  int bl_ndcnt = 0;
  int eptr_cnt = 0;

  for(auto it = elements.begin(); it != elements.end(); ++it){
    const MEl* el = it->get();
    const MEl* symmetry_child_el;
    bool is_bl = false;
    eptrs.push_back(0);
    int symmetry_face=-1;
    for(int ch = 0; ch < el->NumChildren(); ch ++){
      const MEl* child = el->getChild(ch);
      if(child->getBCTag() == 7){
	is_bl = true;

	for(int nd = 0; nd < child->NumNodes(); nd++){
	  int node = el->getNodes()[nd];
	  auto it = blnodemap.find(node);
	  if(it == blnodemap.end()){
	    it = blnodemap.insert(it,std::make_pair(node,bl_ndcnt++));
	    blnodes.push_back(node);
	  }
	  eind.push_back(it->second);
	}
	//bl_faces.push_back(child);
	eptrs.push_back(eptr_cnt+=child->NumNodes());
	//break;
      }
      else if(child->getBCTag() == 20){
	symmetry_face = ch;
	symmetry_child_el = child;
	//std::cout << "Has symmetry face!" << std::endl;
      }
    }
    if(is_bl){
      if(symmetry_face != -1) std::cout << "Has symmetry face" << std::endl;
      blels.push_back(*it);
      symmetry_faces.push_back(symmetry_face);
      if(symmetry_face != -1) symmetry_face_elements.insert(symmetry_child_el);
    }
    //if(is_bl) std::cout << "found bl" << std::endl;
 
  }
}
