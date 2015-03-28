#include "BoundaryLayerSubdivider.h"
#include "BLMesh.h"
#include "BLParameterList.h"
#include "DefaultNodeSpacingEvaluator.h"

#include "InterpolationPoints.h"
#include "PolynomialBasis.h"
#include "NodeInterpolator.h"

#include "MeshContainer.h"
#include "nodeFactory.h"
#include "elementFactory.h"
#include "ChildGenerator.h"

#include <iostream>

int BoundaryLayerSubdivider::GenerateElements(){
  BLMesh bl_mesh(mesh);
  std::cout << "N blnodes: " << bl_mesh.blnodes.size() << std::endl;

  int Nblnodes = bl_mesh.blnodes.size();
  
  double minSpacing = bl_parameters->minSpacing;
  double growthRatio = bl_parameters->growthRatio;
  int NLayers = bl_parameters->NLayers;
  
  
  DefaultNodeSpacingEvaluator node_spacing_evaluator(bl_mesh,minSpacing,
						     growthRatio,NLayers);

  std::vector<double> spacing = node_spacing_evaluator.computeSpacing();
    

  InterpolationPointFactory interpolation_factory;
  PolynomialBasisFactory basis_factory;

  
  const PolynomialBasis& line_basis = basis_factory.getPolynomialBasis(1);

  const InterpolationPoints& line_points_computer = 
    interpolation_factory.getInterpolationPoints(1);

  int order = (*mesh.getElements().begin())->getOrder();
  
  arma::mat interp_curr = line_points_computer.ComputePoints(order,0);

  std::cout << interp_curr << std::endl;

  int Nblpts_ho = NLayers+1 + NLayers*(order-1);
  
  arma::mat blpts(Nblpts_ho,1);
  
  NodeInterpolator node_interpolator(0);
  const arma::mat& lin_interp = 
    node_interpolator.getNodeInterpolationMatrix(1,1,order);

  blpts[0] = 0;
  for(int lay = 0, cnt = 1; lay < NLayers; lay++,cnt++){
    arma::vec2 currbl_pts = {spacing[lay],spacing[lay+1]};
    arma::vec ho_currbl_pts = lin_interp.t()*currbl_pts;
    for(int i = 0; i < (order-1); i++,cnt++) blpts[cnt] = ho_currbl_pts[i+1];
    blpts[cnt] = spacing[lay+1];
  }
  //for(int i = 0; i <=NLayers; i++) blpts[i] = spacing[i];
  
  std::cout << blpts << std::endl;

  arma::mat A, f;
  arma::cube dA, df;
  arma::cube dA2, df2;
  line_basis.EvalBasis(order,interp_curr,A,dA,dA2);
  line_basis.EvalBasis(order,blpts,f,df,df2);

  arma::mat interp = arma::solve(A,f);
  
  const std::vector<unique_element_ptr>& blels = bl_mesh.blels;
  
  int ndof;
  if(mesh.Dimension() == 2) ndof = order+1;
  else ndof = (order+1)*(order+2)/2;
  
  arma::mat linenodes(3,order+1);

  int nadd_line = blpts.size()-2;
  int Nnew_nodes = nadd_line*bl_mesh.blnodes.size();
  arma::mat newnode_coords(3,Nnew_nodes);
  
  std::vector<int> node_geo_entities(bl_mesh.blnodes.size());
  std::vector<int> node_types(bl_mesh.blnodes.size());
  std::vector<int> top_bl_nodes(bl_mesh.blnodes.size());

  for(auto it = blels.begin(); it != blels.end(); ++it){
    const MEl* el = it->get();
    const gind* elnodes = el->getNodes();

    for(int nd = 0; nd < ndof; nd++){
      int currnd = elnodes[nd];
      int local_bl_node = bl_mesh.blnodemap[currnd];
      for(int i = 0; i <=order; i++){
	int linend = elnodes[nd+i*ndof];
	const auto& node_ = mesh.getNodes().at(linend);
	linenodes.col(i) = node_->xyzvec3();
	if(i == order){
	  node_geo_entities[local_bl_node] = node_->getGeoEntity();
	  node_types[local_bl_node] = node_->getType();
	  top_bl_nodes[local_bl_node] = linend;
	}
      }
      //const auto& endnode = mesh.getNodes().at(elnodes
      arma::mat newnodes = linenodes*interp;
      newnodes = newnodes.cols(arma::span(1,newnodes.n_cols-2));
      const arma::span& node_span = 
	arma::span(local_bl_node*nadd_line,(local_bl_node+1)*nadd_line-1);
      newnode_coords.cols(node_span) = newnodes;
    }
  }

  std::vector<int> bl_node_ids(Nblnodes*blpts.size());
  
  node_map& nodes = mesh.getNodesNC();

  for(int i = 0; i < bl_mesh.blnodes.size(); i++){
    bl_node_ids[i] = bl_mesh.blnodes[i];
    for(int j = 0; j < nadd_line; j++){
      int node_id = nodeFactory::Instance()->GetNodeCount();
      bl_node_ids[i + Nblnodes*(j+1)] = node_id;
      
      //std::cout << newnode_coords.col(i*nadd_line+j) << std::endl;

      nodes[node_id] = 
	//auto newnode = 
	nodeFactory::Instance()->
	CreateNode(node_types[i],
		  newnode_coords.colptr(i*nadd_line+j),
		  node_geo_entities[i]);
      //nodes[node_id] = newnode;

    }
    bl_node_ids[Nblnodes*(blpts.size()-1)+i] = top_bl_nodes[i];
  }
  

  std::vector<unique_element_ptr> newelements;
  for(auto it = blels.begin(); it != blels.end(); ++it){
    const MEl* el = it->get();
    const gind* elnodes = el->getNodes();
    
    for(int lay = 0, lay_offset=0; lay < NLayers; lay++){
      std::vector<int> node_ids(ndof*(order+1));
      for(int nd = 0; nd < ndof; nd++){
	int currnd = elnodes[nd];
	int local_bl_node = bl_mesh.blnodemap[currnd];
	for(int i = 0; i <=order; i++){
	  node_ids[nd+i*ndof] = 
	    bl_node_ids[local_bl_node+Nblnodes*(lay_offset+i)];

	  if(i > 0 && i < order){
	    auto it = mesh.getNodes().find(elnodes[nd+i*ndof]);
	    if(it != mesh.getNodes().end()){
	      mesh.getNodesNC().erase(it);
	    }
	  }
	}
      }
      elementFactory* element_factory = elementFactory::Instance();
      int eltype = el->getElementType();
      
      unique_element_ptr newel = 
	element_factory->CreateElement(eltype,node_ids,0,order,el->getBCTag(),
				       el->getGeoEntity(),el->getGeoType());
      newelements.push_back(newel);
    
      lay_offset+= order;
    }
    mesh.getElementsNC().erase(*it);
  }
  
  if(mesh.Dimension() == 3){
    ChildGenerator child_generator;
    for(int i = 0; i < blels.size(); i++){
      if(bl_mesh.symmetry_faces[i] != -1){
	for(int lay = 0; lay < NLayers; lay++){
	  int symmetry_face = bl_mesh.symmetry_faces[i];
	  std::vector<element_ptr>& children = 
	    child_generator.GenerateChildren(newelements[i*NLayers+lay]);
	  children[symmetry_face]->setBCTag(20);
	  children[symmetry_face]->setGeoEntity(mesh.SymmetrySurface());
	  mesh.getSubElementsNC().insert
	    (std::move(children[symmetry_face]));
	}
      }
    }

    // Delete the old symmetry face elements
    for(auto it = mesh.getSubElementsNC().begin(); 
	it != mesh.getSubElementsNC().end(); ){
      auto itchild = bl_mesh.symmetry_face_elements.find(it->get());
      if(itchild != bl_mesh.symmetry_face_elements.end()){
	it = mesh.getSubElementsNC().erase(it);
      }
      else{
	++it;
      }
    }
  }
  //auto it = mesh.getElements().find(std::weak_ptr<MEl>(bl_mesh.blels[0].get()));

  for(auto it = newelements.begin(); it != newelements.end(); ++it){
    mesh.getElementsNC().insert(std::move(*it));
  }

 
}
