#include "MeshOptimizer.h"
#include "MeshContainer.h"
#include "assert.h"
#include "OptElManager.h"


MeshOptimizer::MeshOptimizer(MeshContainer& mesh,
			     OptElManager& optel_manager,
			     const std::map<const MEl*,arma::mat>& ideal_map,
			     int el_dim_to_opt): 
  mesh(mesh), optel_manager(optel_manager), idealElements(ideal_map),
  el_dim_to_opt(el_dim_to_opt){

  //for(auto el = idealElements.begin(); el != idealElements.end(); ++el){
  //  std::cout << el->first->getDim() << std::endl;
  //}
  arma::mat gradMerit;
  std::cout << "initial ideal size: " << idealElements.size() << std::endl;
  
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*, double>& all_merits = mesh.getAllMeritsNC();
  if(all_merits.size() < idealElements.size()){
    for(auto el = idealElements.begin(); el != idealElements.end(); ++el){
      std::unique_ptr<OptEl> optel = 
	optel_manager.CreateOptEl(el->first,el->second);
      double detS;
      optel->computeGradMerit(gradMerit,all_merits[el->first],detS);
    }
  }
  else if(all_merits.size() > idealElements.size()){
    assert(all_merits.size() > idealElements.size());
  }

  std::cout << "All merits size: " << all_merits.size() << std::endl;
}

std::map<gind,std::vector<const MEl*> > MeshOptimizer::
CreateNode2ElementList(const std::unordered_set<const MEl*>& watch,
		       const std::map<gind,gind>& activenodes){

  std::map<gind,std::vector<const MEl*> > Node2ElementList;


  for(auto el = watch.begin(); el != watch.end(); ++el){
    const MEl* element = *el;
    const gind* elnodes = element->getNodes();
    const int nn = element->NumNodes();

    //const gind* cn = (*el)->getCornerNodes();
    //const int ncn = (*el)->numCornerNodes();
    for(int i=0; i<nn; i++){
      auto it = activenodes.find(elnodes[i]);
      if(it != activenodes.end()){
	Node2ElementList[it->first].push_back(*el);
      }
    }
  }

  return Node2ElementList;
}


void MeshOptimizer::InsertActiveNodes(const MEl* el, 
		       std::map<gind,gind>& activenodes,
		       const node_map& nodes){
  
  //only_opt_node_with_bctag;

  const gind* elnodes = el->getNodes();
  const int nn = el->NumNodes();
 
  
  for(int i=0; i<nn; i++){
    int type = nodes.at(elnodes[i])->getType();
    if(dims_to_opt[type]){
      bool opt_face = true;
      if(geo_faces_to_opt.find(nodes.at(elnodes[i])->getGeoEntity()) ==
	 geo_faces_to_opt.end() && geo_faces_to_opt.size() != 0){
	opt_face = false;
      }
      if(opt_face || type == 3){
	if(clamped_nodes.find(elnodes[i]) == clamped_nodes.end())
	  activenodes[elnodes[i]] = 0;
      }

    }
  }
}
bool HasActiveNode(const MEl* el, 
		   std::map<gind,gind>& activenodes){
  //const gind* cn = el->getCornerNodes();
  //const int ncn = el->numCornerNodes();
  const gind* elnodes = el->getNodes();
  const int nn = el->NumNodes();

  for(int i=0; i<nn; i++){
    if(activenodes.find(elnodes[i]) != activenodes.end()) return true;
  }
  return false;
}

void NumberActiveNodes(std::map<gind,gind>& activenodes){
  int cnt=0;
  for(auto it = activenodes.begin(); it != activenodes.end(); ++it, ++cnt){
    it->second = cnt;
  }
}

void MeshOptimizer::
FindActiveElements(std::unordered_set<const MEl*>& watch,
		   std::map<gind,gind>& activenodes,
		   const double threshold,
		   const int Nnearest,
		   const int type_to_find,
		   const int type_to_exclude){
  watch.clear();
  activenodes.clear();


  const node_map& nodes = mesh.getNodes();
  //const element_set& elements = mesh.getElements();
  //const element_set& subelements = mesh.getSubElements();
  //const ElDouble_map& all_merits = mesh.getAllMerits();

  //ideal_map& idealElements = mesh.getIdealElementsNC();
  
  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    //if(el->first->IsBL()) std::cout << "We have a BL element!" << std::endl;
  }
  //std::cout << "threshold: " << threshold << std::endl;
  //std::cout << "type to find: " << type_to_find << std::endl;
  //std::cout << "type to exclude: " << type_to_exclude << std::endl;
  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    if(el->first->getDim() == el_dim_to_opt){
      //std::cout << el->second << " " << threshold << std::endl;
      if((el->second > threshold || 
	  el->first->getElementType() == type_to_find) 
	 && el->first->getElementType() != type_to_exclude){
	watch.insert(el->first);
      }
      else if(el->first->IsBL() && el->second > (1.0 + 1.0e-10 )){
	//watch.insert(el->first);
      }
    }
  }

  //std::cout << "before layer loop" << std::endl;
  //std::cout << "Nnear: " << Nnearest << " " << watch.size() << std::endl;

  for(int lay=0; lay<=Nnearest; lay++){
    for(auto el = watch.begin(); el != watch.end(); ++el){
      //std::cout << "In watch" << std::endl;
      InsertActiveNodes(*el,activenodes,nodes);
    }

    const element_set& elements = mesh.getElementsOfDim(el_dim_to_opt);
    for(auto el = elements.begin(); el != elements.end(); ++el){
      MEl* element = el->get();
      int el_dim = element->getDim();
      if((el_dim == 2 && element->hasGeoEntity()) || element->getDim() == 3){
	if(HasActiveNode(element,activenodes)){
	  watch.insert(element);
	  auto it = idealElements.find(element);
	  if(it == idealElements.end()){
	    idealElements[element] = 
	      optel_manager.CreateIdealMatrix(element);
	  }
	}
      }
    }

      /*
    if(mesh.MeshDimension() == el_dim_to_opt){
      for(auto el = elements.begin(); el != elements.end(); ++el){
	if(HasActiveNode(el->get(),activenodes)){
	  watch.insert(el->get());
	  auto it = idealElements.find(el->get());
	  if(it == idealElements.end()){
	    idealElements[el->get()] = 
	      optel_manager.CreateIdealMatrix(el->get());
	  }
	}
      }
    }
    
    if(mesh.MeshDimension()-1 == el_dim_to_opt){
      for(auto el = subelements.begin(); el != subelements.end(); ++el){
	//cout << (*el)->hasGeoEntity() << endl;
	if((*el)->hasGeoEntity()){
	  if(HasActiveNode(el->get(),activenodes)){
	    watch.insert(el->get());
	    auto it = idealElements.find(el->get());
	    if(it == idealElements.end()){
	      idealElements[el->get()] = 
		optel_manager.CreateIdealMatrix(el->get());
	    }
	  }
	}
      }
    }
      */
  }
 
  NumberActiveNodes(activenodes);

}
int MeshOptimizer::
UpdateMeshMerits(const std::unordered_set<MEl*>& watch){
  //ElDouble_map& all_merits = mesh.getAllMeritsNC();
  //ideal_map& idealElements = mesh.getIdealElementsNC();

  arma::mat gradMerit;

  for(auto el = watch.begin(); el != watch.end(); ++el){
    std::unique_ptr<OptEl> optel = 
      optel_manager.CreateOptEl(*el,idealElements.at(*el));
    double detS;
    optel->computeGradMerit(gradMerit,all_merits[*el], detS);
    //all_merits[*el] = optel->computeDistortion();
    
  }

}
int MeshOptimizer::UpdateMeshMerits(){
  //ElDouble_map& all_merits = mesh.getAllMeritsNC();
  //ideal_map& idealElements = mesh.getIdealElementsNC();
  arma::mat gradMerit;
  for(auto el = idealElements.begin(); el != idealElements.end(); ++el){
    std::unique_ptr<OptEl> optel = 
      optel_manager.CreateOptEl(el->first,el->second);
    double detS=-1;
    optel->computeGradMerit(gradMerit,all_merits[el->first], detS);

  }

  return 1;
}
const double MeshOptimizer::minMeshQuality() const{
  double min_quality = 1.0;
  //const std::map<MEl*,double>& all_merits = mesh.getAllMerits();

  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    if(el->first->getDim() == el_dim_to_opt){
      //if(dims_to_opt[el->first->getDim()]){
      min_quality = std::min(min_quality,1.0/el->second);
    }
  }
  return min_quality;
}
