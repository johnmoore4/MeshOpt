#include "HighOrderMeshGenerator.h"
#include "MeshContainer.h"
#include "GeometryContainer.h"
#include "Mnode.h"
#include "Mel.h"
#include "NodeIndexer.h"
#include "NodeInterpolator.h"
#include "MEl_to_ActiveEl.h"
#include "ActiveMEl.h"
#include "ParametricCoordinateEvaluator.h"
#include "nodeFactory.h"
#include "elementFactory.h"

HighOrderMeshGenerator::
HighOrderMeshGenerator(MeshContainer& mesh, GeometryContainer&geometry):
  mesh(mesh), geometry(geometry){
  meshcurr = std::make_shared<MeshContainer>();
  meshold = std::make_shared<MeshContainer>();
  index_factory = std::make_shared<NodeIndexFactory>();
  node_interpolator = std::make_shared<NodeInterpolator>(0);
  //parametric_coordinate_evaluator = 
  //  std::make_shared<ParametricCoordinateEvaluator>();
}

std::shared_ptr<MeshContainer>
HighOrderMeshGenerator::generateHighOrderMeshRecursive(int order){
  int curr_order = (*mesh.getElements().begin())->getOrder();

  //order = 2;

  for(int ord = curr_order+1; ord <= order; ord++){
    //  for(int ord = order; ord <= order; ord++){
    *meshcurr = mesh;
    meshcurr->getElementsNC().clear();
    meshcurr->getSubElementsNC().clear();
    meshcurr->getSubSubElementsNC().clear();


    for(int dim = 1; dim <= mesh.Dimension(); dim++){
      std::cout << "starting to mesh: " << dim << " " << ord << std::endl;
      meshDimension(dim,ord);
      std::cout << "finished meshing dim: " << dim << std::endl;
    }
    setNewChildrenAndOrientation();
    *meshold = *meshcurr;
    mesh = *meshold;
    
    for(int d = 0; d <= mesh.Dimension(); d++) elmap[d].clear();

    //mesh = *meshcurr;
  }

  return meshcurr;

}

int HighOrderMeshGenerator::setNewChildrenAndOrientation(){
  for(int dim = 2; dim <= mesh.Dimension(); dim++){
    const element_set& elements = mesh.getElementsOfDim(dim);
    //const element_set& elements_new = meshcurr->getElementOfDim(dim);
    for(auto el = elements.begin(); el != elements.end(); ++el){
      auto elit = elmap[dim].find(el->get());
      assert(elit != elmap[dim].end());
      MEl* elnew = elit->second;
      const MEl* elold = el->get();
      for(int ch = 0; ch < elold->NumChildren(); ch++){
	MEl* newchild = const_cast<MEl*>(elmap[dim-1][elold->getChild(ch)]);
	assert(newchild);
	signed char child_or = elold->getChildOrientation(ch);
	elnew->setChild(ch,newchild,child_or);
      }
    }
  }

}

int HighOrderMeshGenerator::meshDimension(int dim, int order){
  element_set& elements = mesh.getElementsOfDim(dim);
  for(auto it = elements.begin(); it != elements.end(); ++it){
    const unique_element_ptr& el = *it;

    MEl_to_ActiveEl m2a(mesh.getNodes(),index_factory);
    const arma::mat& elnodes = m2a.get(el).getNodesMatrix();
    const arma::mat& interpolation_matrix = 
      node_interpolator->getNodeInterpolationMatrix(el->getElementType(),
						    el->getOrder(),order);

    arma::mat newelnodes = elnodes*interpolation_matrix;
    arma::mat new_parametric_coords;

    //std::cout << "before has geo entity" << std::endl;
    //if(el->getGeoEntity() >= 0){
    if(el->hasGeoEntity()){
      //std::cout << "in has geo entity" << std::endl;
      ParametricCoordinateEvaluator 
	parametric_evaluator(geometry,mesh.getNodes());

      const arma::mat& parametric_coords = 
	parametric_evaluator.getParametricCoordinates(el);


      const arma::mat& newparametric = 
	parametric_coords*interpolation_matrix;
      

      newelnodes = 
	parametric_evaluator.projectPointsOnGeometry(el,
						     newelnodes,
						     newparametric,
						     new_parametric_coords);
  
    }

    //std::cout << "After has geo entity" << std::endl;
    
    std::vector<int> nodeindices = 
      insertNewNodes(el,newelnodes,new_parametric_coords,order);


    /*
    for(int i = 0; i < nodeindices.size(); i++){
      std::cout << nodeindices[i] << " ";
    }
    std::cout << std::endl;
    */

    insertNewElement(el,nodeindices,order);

  }

};

std::vector<double> HighOrderMeshGenerator::
getElementInteriorNodesCoordinates(unique_element_ptr& el){
  int nc = el->NumChildren();
  for(int ch = 0; ch < nc; ch++){
    const MEl* child = el->getChild(ch);
    if(child){
      
    }
    else{

    }
  }
}

arma::mat getHighOrderNodeCoordinates(unique_element_ptr& el, int order){
  

}

std::vector<int> HighOrderMeshGenerator::
insertNewNodes(const unique_element_ptr& el,
	       const arma::mat& newelnodes,
	       const arma::mat& parametric_coords,
	       int order){

  const int type = el->getElementType();

  const NodeIndexer* ni_old = 
    index_factory->getNodeIndexer(type,el->getOrder());

  const NodeIndexer* ni_new = 
    index_factory->getNodeIndexer(type,order);
  
  std::vector<int> newnode_indices(ni_new->Ndof(),-1);
  const std::vector<indtype>& corner_new = ni_new->getCornerNodes();
  
  const gind* old_nodes = el->getNodes();
  const std::vector<indtype>& corner_old = ni_old->getCornerNodes();
  
  for(int i = 0; i < el->numCornerNodes(); i++){
    //std::cout << "corner node: " << corner_new[i] << " " << 
    //  corner_old[i] << std::endl;
    newnode_indices[corner_new[i]] = old_nodes[corner_old[i]];
  }
 


    //std::cout << el->getDim() << std::endl;
  for(int ch = 0; ch < el->NumChildren(); ch++){
    const MEl* child_old = el->getChild(ch);
    if(el->getDim() > 1) assert(child_old);
    if(child_old){
      //int child_dim = child_old->getDim();
      auto child_it = elmap[el->getDim()-1].find(child_old);
      assert(child_it != elmap[el->getDim()-1].end());
   
      const MEl* child_new = child_it->second;
      const gind* child_nodes = child_new->getNodes();
      int orient = el->getChildOrientation(ch);
      
      const std::vector<indtype>& child_node_indices = 
	ni_new->getChildNodes(ch);

      const NodeIndexer* ni_child =
	index_factory->getNodeIndexer(child_new->getElementType(),
				      child_new->getOrder());

      const std::vector<indtype>& orient_indices = 
	ni_child->getOrientedNodes(orient);
      assert(child_node_indices.size() == child_new->NumNodes());
      for(int nd = 0; nd < child_new->NumNodes(); nd++){
	newnode_indices[child_node_indices[nd]] = 
	  child_nodes[orient_indices[nd]];
      }
    }
  } // end child loop

  
  // delete the old interior nodes from map
  
  const std::vector<indtype>& old_interior_indices = ni_old->getInteriorNodes();
  for(int i = 0; i < old_interior_indices.size(); i++){
    meshcurr->getNodesNC().erase(old_nodes[old_interior_indices[i]]);
  }
  

  // insert new nodes into node map
  //std::cout << "ni_new type: " << ni_new->getType() << std::endl;
  //arma::uvec new_interior_indices;
  //std::cout << "el type: " << el->getElementType() << std::endl;
  arma::uvec new_interior_indices(ni_new->getInteriorNodes());
  //std::cout << "h3.1" << std::endl;
  //std::cout << new_interior_indices << std::endl;
  //std::cout << newelnodes << std::endl;
  arma::mat new_interior_nodes = newelnodes.cols(new_interior_indices);
  //std::cout << "h3.2" << std::endl;
  //std::cout << new_interior_indices << std::endl;
  //std::cout << parametric_coords << std::endl;
  arma::mat new_interior_parametric_coords;
  if(!parametric_coords.is_empty()){
    new_interior_parametric_coords = 
      parametric_coords.cols(new_interior_indices);
  }

  //std::cout << "h4" << std::endl;

  nodeFactory* node_factory = nodeFactory::Instance();

  int node_type = 3;
  if(el->hasGeoEntity()) node_type = el->getGeoType()+1;
     
  //std::cout << new_interior_nodes << std::endl;

  for(int i = 0; i < new_interior_indices.size(); i++){
    double uv[2];
    for(int j = 0; j < parametric_coords.n_rows; j++){
      uv[j] = new_interior_parametric_coords.at(j,i);
    }
    int nd_count = node_factory->GetNodeCount();
    newnode_indices[new_interior_indices[i]] = nd_count;
    auto newnode = 
      node_factory->CreateNode(node_type,
			       new_interior_nodes.colptr(i),
			       el->getGeoEntity(),
			       uv[0],uv[1]);
    meshcurr->getNodesNC()[nd_count] = std::move(newnode);

  }

  //std::cout << "h5" << std::endl;

  return newnode_indices;
  
}

int HighOrderMeshGenerator::
insertNewElement(const unique_element_ptr& el,
		 std::vector<int>& nodeindices,
		 int order){

  const NodeIndexer* ni_new = 
    index_factory->getNodeIndexer(el->getElementType(),order);

  elementFactory* element_factory = elementFactory::Instance();
  
  element_set& elements = meshcurr->getElementsOfDim(el->getDim());

  auto newel = element_factory->CreateElement(el->getElementType(),
					      nodeindices,
					      ni_new,
					      order,
					      el->getBCTag(),
					      el->getGeoEntity(),
					      el->getGeoType());
  

  elmap[el->getDim()][el.get()] = newel.get();
  

  elements.insert(std::move(newel));


}

int HighOrderMeshGenerator::setMeshIsCurved(){
  /*
  element_set& elements = mesh.getElementsOfDim(1);
    for(auto el = elements.begin(); el != elements.end(); ++el){
      MEl* element = el->get();
      if(element->hasGeoEntity) element->setIsCurved(true);
    }

  for(int dim = 1; dim <= mesh.Dimension(); dim++){


  }  
  */
}
