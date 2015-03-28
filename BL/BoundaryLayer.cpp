#include "BoundaryLayer.h"
#include "MeshContainer.h"
#include "GeometryContainer.h"
#include "ChildGenerator.h"

#include "nodeFactory.h"
#include "elementFactory.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"
#include "NodeReParameterizer.h"
#include "BLParameterList.h"
#include "BoundaryLayerSubdivider.h"

#include "OptElManager.h"

#include "ActiveMEl.h"
#include "OptEl.h"

#include "Mel.h"
#include "Mnode.h"
#include <algorithm>

int BoundaryLayerGenerator::GenerateBL(){
 
  double thickness = bl_parameters->thickness;
  //int Nlayers = bl_parameters->NLayers;
  int Nlayers = 1;
  double factor = 0.01;
  const double safety_factor = 2;
  const double max_normal_thickness_factor = 10;

  findSymmetryFaces();

  std::vector<MEl*> bl_elements = FindBLElements(mesh.getSubElements());

  if(bl_elements.size() == 0) return 1;

  

  // Get node normals
  //normal_map_type normal_map = CreateNodeInfoMap(bl_elements);
 
  normal_map = CreateNodeInfoMap(bl_elements,thickness);

  // Map original nodes to new nodes
  gindmap newnode_map = CreateNewNodeMap(normal_map);

  std::vector<double> max_extrude = 
    getMaximumSafeExtrusionDistance(safety_factor);

  std::vector<double> extrusion_dist = max_extrude;
  for(auto it = extrusion_dist.begin(); it != extrusion_dist.end(); ++it){
    double max_normal_dist = *it*safety_factor;
    double extrude_dist = std::min(max_normal_dist*max_normal_thickness_factor,
				   thickness);
    *it = extrude_dist;
    //double extrude_dist = std::max(*it*safety_factor,m
    //(*it)*= safety_factor;
  }

  // Make Line Symmetry elements
  std::vector<MEl*> symmetry_elements = MakeSymmetryLineElements();

 // Replace elements containing moved boundary node
  ReplaceNeighboringElements(mesh.getElementsNC(), Nlayers);
  ReplaceNeighboringElements(mesh.getSubElementsNC(), Nlayers,20);

 
 
 // Insert the BL nodes
  InsertBLNodes(1,extrusion_dist);
  std::cout << "after insert BL nodes" << std::endl;

  //ResetBLNodes(1,extrusion_dist);

 //ReplaceSymmetryElements(newnode_map,normal_map,Nlayers);

 cout << "symmetry elements size: " << symmetry_elements.size() << endl;

 // Insert the BL elements
 InsertBLElements(mesh.getElementsNC(),bl_elements,Nlayers,thickness);
 InsertBLElements(mesh.getSubElementsNC(),symmetry_elements,Nlayers,
 		  thickness,true);



 // Reset the BL nodes
 InsertBLNodes(1.0/safety_factor,max_extrude);
 //InsertBLNodes(normal_map,Nlayers,thickness,factor,max_extrude);
 //ResetBLNodes(normal_map,newnode_map,Nlayers,thickness,factor);
 
 for(auto nd = mesh.getNodes().begin(); nd != mesh.getNodes().end(); ++nd){
   //std::cout << nd->second->getGeoEntity() << std::endl;
 }
 std::cout << "At end of genereate BL" << std::endl;

  return 1;

}

const arma::vec3 El1DNormal(const gind* cn, const node_map& nodes, 
		      const arma::vec3& normal){

  arma::vec3 r1 = nodes.at(cn[1])->xyzvec3() - nodes.at(cn[0])->xyzvec3();
  arma::vec3 elnorm = cross(r1,normal);
  return elnorm/=norm(elnorm,2.0);
}

const arma::vec3 El2DNormal(const gind* cn, const node_map& nodes){
  arma::vec3 r1 = nodes.at(cn[1])->getXYZ() - nodes.at(cn[0])->getXYZ();
  arma::vec3 r2 =nodes.at(cn[2])->getXYZ() - nodes.at(cn[0])->getXYZ();
  arma::vec3 elnorm = cross(r1,r2);
  return elnorm/=norm(elnorm,2.0);

}		      

int BoundaryLayerGenerator::findSymmetryFaces(){
  if(mesh.Dimension() == 3){
    //std::cout << "finding symmetry faces..." << std::endl;
    auto& subelements = mesh.getSubElements();
    for(auto it = subelements.begin(); it != subelements.end(); ++it){
      const MEl* el = it->get();
      if(el->getBCTag() == 20) symmetry_faces.insert(el->getGeoEntity());
    }
  }
  return 0;
}

int BoundaryLayerGenerator::
SmoothField(double Nsmooth, const std::vector<MEl*>& bl_elements,
	    std::vector<double>& field){
  std::vector<std::set<int> > node_neighbors(normal_map.size());
  for(auto it = bl_elements.begin(); it != bl_elements.end(); ++it){
    MEl* el = *it;
    const gind* elnodes = el->getNodes();
    for(int i = 0; i < el->NumNodes(); i++){
      int ni = elnodes[i];
      for(int j = 0; j < el->NumNodes(); j++){
	int nj = elnodes[j];
	node_neighbors[ni].insert(nj);
      }
    }
  }
  
  for(int smooth = 0; smooth < Nsmooth; smooth++){
    for(int nd = 0; nd < node_neighbors.size(); nd++){
      
    }
  }
}
std::vector<MEl*> BoundaryLayerGenerator::
FindBLElements(const element_set& elements){
  // const element_set& elements = mesh.getSubElements();
  std::vector<MEl*> bl_elements;
  for(auto el = elements.begin(); el != elements.end(); ++el){
    const int bc_tag = (*el)->getBCTag();
    if(bc_tag == 7){
      //if(bc_tag >= 10 && bc_tag < 20){
      bl_elements.push_back(el->get());
    }
  }
  return bl_elements;
}
std::vector<MEl*> BoundaryLayerGenerator::MakeSymmetryLineElements(){
  const element_set& elements = mesh.getSubElements();
  std::vector<MEl*> line_elements;
  ChildGenerator child_generator;
  for(auto el = elements.begin(); el != elements.end(); ++el){
    if((*el)->getBCTag() == 20){
      std::vector<unique_element_ptr>& childels = 
	child_generator.GenerateChildren(*el);
      for(auto ch = childels.begin(); ch != childels.end(); ++ch){
	const int ncn = (*ch)->numCornerNodes();

	const gind* cn = (*ch)->getCornerNodes();
	for(int i=0, cnt=0; i<ncn; i++){
	  if(newnode_map.find(cn[i]) != newnode_map.end()) cnt++;
	  if(cnt == 2){
	    store_symmetry.push_back(std::move(*ch));
	    line_elements.push_back(store_symmetry.back().get());
	    //std::cout << "setting symmetry bl entity to: " << 
	    //  (*el)->getGeoEntity() << std::endl;
	    line_elements.back()->setGeoEntity((*el)->getGeoEntity());
	    line_elements.back()->setBCTag(20);
	    //cout << "found a BL line element!" << endl;
	  }
	}
      }
    }
  }


  return line_elements;
}

normal_map_type BoundaryLayerGenerator::
CreateNodeInfoMap(const std::vector<MEl*>& bl_elements,
		  double thickness){

  arma::vec3 n1 = {0.0, 0.0, 1.0};

  const node_map& nodes = mesh.getNodes();
  const int mesh_dim = mesh.MeshDimension();

  typedef std::pair<normal_map_type::iterator,bool> retpair;
  const int dim = (*bl_elements.begin())->getDim();

  std::set<int> bl_geo_entities;
  for(auto it = bl_elements.begin(); it != bl_elements.end(); ++it){
    auto el = *it;
    //int geo_entity = el->getGeoEntity();
    bl_geo_entities.insert(el->getGeoEntity());
  }

  normal_map_type normal_map;
  for(auto el = bl_elements.begin(); el != bl_elements.end(); ++el){
    const gind* cn = (*el)->getCornerNodes();
    //const gind elnodes = (*el)->getNodes();
    int ncn = (*el)->numCornerNodes();
    std::vector<arma::vec3> node_normals(ncn);
    arma::vec3 n;
    arma::vec3 weight;
    if(dim == 1){
      n = El1DNormal(cn,nodes,n1);
      for(int i = 0; i < ncn; i++) node_normals[i] = n;
      //weight = {1.0, 1.0, 1.0};
    }
    else if(dim == 2){
      n = normalFromElnodes(cn,nodes);
 
      const unsigned int ind[3][2] = {{1,2},{2,0},{0,1}};
      for(int i=0; i<ncn; i++){
	int type = nodes.at(cn[i])->getType();
	
	node_normals[i] = n;
	/*
	if(type == 0){
	  const std::map<int,std::set<int> >& Vertex2FacesMap = 
	    geometry.getVertex2FacesMap();
	  int vertex = nodes.at(cn[i])->getGeoEntity();
	  myVertex& my_vertex = geometry.getVerticesNC()[vertex];
	  
	  auto face_set = Vertex2FacesMap.at(vertex);
	  arma::vec3 average_normal;
	  average_normal.zeros();
	  for(auto fc = face_set.begin(); fc != face_set.end(); ++fc){
	    if(bl_geo_entities.find(*fc) != bl_geo_entities.end()){
	      myFace& my_face = geometry.getFacesNC()[*fc];
	      arma::vec2 params = my_face.paramsOnFace(&my_vertex);
	      const arma::mat& deriv_temp = my_face.D1(params.memptr());
	      arma::mat deriv = deriv_temp.t();
	      arma::vec3 normal = arma::cross(deriv.col(0),deriv.col(1));
	      normal/= arma::norm(normal,2);
	      average_normal+= normal;
	    }
	  }
	  average_normal/= arma::norm(average_normal,2);
	  node_normals[i] = average_normal;
	}
	else if(type == 1){
	  const std::map<int,std::set<int> >& Edge2FacesMap = 
	    geometry.getEdge2FacesMap();


	  int edge = nodes.at(cn[i])->getGeoEntity();
	  double u = *nodes.at(cn[i])->getParametricCoords();
	  myEdge& my_edge = geometry.getEdgesNC()[edge];
	  auto face_set = Edge2FacesMap.at(edge);
	  arma::vec3 average_normal;
	  average_normal.zeros();
	  for(auto fc = face_set.begin(); fc != face_set.end(); ++fc){
	    if(bl_geo_entities.find(*fc) != bl_geo_entities.end()){
	      myFace& my_face = geometry.getFacesNC()[*fc];
	      arma::vec2 params = my_face.paramsOnFace(&my_edge,u);
	      //my_face.D1(params.memptr(),deriv);
	      const arma::mat& deriv_temp = my_face.D1(params.memptr());
	      arma::mat deriv = deriv_temp.t();
	      arma::vec3 normal = arma::cross(deriv.col(0),deriv.col(1));
	      normal/= arma::norm(normal,2);
	      average_normal+= normal;
	    }
	    //if(dot(normal,n) > 0) node_normals[i] = normal;
	    //else node_normals[i] = -normal;
	    //node_normals[i] = -normal;

	    //std::cout << "doing surface evaluation of normal!" << std::endl;
	  }
	  average_normal/= arma::norm(average_normal,2);
	  node_normals[i] = average_normal;
	}
	*/
	arma::vec3 r1 = nodes.at(cn[ind[i][0]])->xyzvec3() 
	  - nodes.at(cn[i])->xyzvec3();
	arma::vec3 r2 = nodes.at(cn[ind[i][1]])->xyzvec3() 
	  - nodes.at(cn[i])->xyzvec3();
	double theta = acos(dot(r1,r2)/(norm(r1,2.0)*norm(r2,2.0)));
	if((*el)->getDim() == 1) theta = 1.0;
	weight(i) = theta/(2.0*arma::datum::pi);
      }
    }
    else{
      cout << "dimension is > 2!" << endl;
    }

    for(int i=0; i<(*el)->numCornerNodes(); i++){
      const int node_dim = nodes.at(cn[i])->getType();
      retpair ret = 
	normal_map.insert(std::make_pair(cn[i],BLNodeInfo(node_normals[i],
							  node_dim,
							  thickness)));
      if(ret.second == false){
	
	//ret.first->second.normal+= n*weight(i);
      }
    }

  }
  for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
    nd->second.normal/=arma::norm(nd->second.normal,2.0);
  }

  if(mesh_dim != 3) return normal_map;

  // Project to symmetry plane
  const std::map<int,std::set<int> >& Edge2FacesMap = 
    geometry.getEdge2FacesMap();


  const std::map<int,std::set<int> >& Vertex2FacesMap = 
    geometry.getVertex2FacesMap();

  //const int symmetry_surf = mesh.SymmetrySurface();
  arma::vec3 symmetry_normal = {0.0, 1.0, 0.0};

  for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
    nd->second.on_symmetry = false;
    for(auto sym = symmetry_faces.begin(); sym != symmetry_faces.end(); ++sym){
      //std::cout << "symmetry faces: " << *sym << std::endl;
      const short int geo_entity = 
	static_cast<MGeoNode*>(nodes.at(nd->first).get())->getGeoEntity();
      //int symmetry_eneity = -1;
      if(nd->second.node_dim == 1){
	auto it = Edge2FacesMap.at(geo_entity).find(*sym);
	if(it != Edge2FacesMap.at(geo_entity).end()){
	  nd->second.on_symmetry = true;
	  nd->second.symmetry_entity = *sym;
	}
      }
      else if(nd->second.node_dim == 0){
	auto it = Vertex2FacesMap.at(geo_entity).find(*sym);
	if(it != Vertex2FacesMap.at(geo_entity).end()){
	  nd->second.on_symmetry = true;
	  nd->second.symmetry_entity = *sym;
	}
      }
    }
    if(nd->second.on_symmetry){
      //std::cout << "node is on symmetry face: " << nd->second.symmetry_entity
      //	<< std::endl;
      nd->second.node_dim = 2;
      arma::vec3 n = nd->second.normal;
      n-= arma::dot(n,symmetry_normal);
      n/=arma::norm(n,2.0);
      nd->second.normal = n;
    } 
  }

  return normal_map;
}

gindmap BoundaryLayerGenerator::
CreateNewNodeMap(const normal_map_type& normal_map){
  const node_map& nodes = mesh.getNodes();
  NNorig = nodes.size();

 
  int cnt=0;
  for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd,cnt++){
    newnode_map[nd->first] = NNorig+cnt;
    int type = nodes.at(nd->first)->getType();
    int entity = nodes.at(nd->first)->getGeoEntity();
    if(type == mesh.Dimension()-2 && entity == 13 ||
       type == mesh.Dimension() - 3 && (entity == 9 || entity == 10)){
      //clamped_nodes.insert(NNorig+cnt);
    }
  }
  return newnode_map;
}

void BoundaryLayerGenerator::
ReplaceSymmetryElements(const std::map<gind,gind>& newnode_map,
			const normal_map_type& normal_map,
			const int Nlayers){
  const int ind[3][2] = {{0,1},{1,2},{2,0}};

  element_set& subelements = mesh.getSubElementsNC();
  node_map& nodes = mesh.getNodesNC();
  //ideal_map& idealElements = mesh.getIdealElementsNC();
  auto& idealElements = ideal_elements;

  const gind Nbln = newnode_map.size();
  std::vector<unique_element_ptr> newels;

  for(auto el = subelements.begin(); el != subelements.end(); ){
    bool found = false;
    if((*el)->getBCTag() == 20){
      const int ncn = (*el)->numCornerNodes();
      const gind* cn = (*el)->getCornerNodes();
      gind newnodes[ncn];
      int nf = 0;
      for(int i=0; i<ncn; i++){
	newnodes[i] = cn[i];
	auto it = normal_map.find(cn[i]);
	if(it != normal_map.end()){
	  found = true;
	  newnodes[i] = newnode_map.at(cn[i]) + (Nlayers-1)*Nbln;
	  nf++;
	}
      }
      if(found){
	
	int type = (*el)->getElementType();
	int geo_entity = (*el)->getGeoEntity();
	int bc_tag = (*el)->getBCTag();
	bool geo_type = (*el)->getGeoType();
	/*
	newels.push_back(elementFactory::Instance()->CreateElement
			 (type,newnodes,0,bc_tag,geo_entity,geo_type));
	
	ActiveMEl activeEl(newels.back().get(),index_factory,nodes,NULL);
	arma::mat ideal = activeEl.getNodesMatrix();
	idealElements[newels.back().get()] = ideal;
	*/
	if(nf == 2){
	  //gind bn[4];
	  std::vector<gind> bn(4);
	  for(int i=0; i<ncn; i++){
	    int cnt=0;
	    for(int j=0; j<2; j++){
	      auto it = normal_map.find(cn[ind[i][j]]);
	      if(it != normal_map.end()){
		bn[j] = it->first;
		bn[2+j] = newnode_map.at(it->first) + (Nlayers-1)*Nbln;
		cnt++;
	      }
	    }
	    if(cnt == 2) break;
	  }


	  
	  newels.push_back(elementFactory::Instance()->CreateElement
			   (3,bn,0,1,bc_tag,geo_entity,geo_type));
	  
	}

	el = subelements.erase(el);	
      }
    }
    if(found == false) ++el;
  }
  
  //node_map& nodes = mesh.getNodesNC();

  for(auto el=newels.begin(); el != newels.end(); ++el){
    ActiveMEl activeEl(el->get(),index_factory,nodes,NULL);
    MEl* currel = el->get();
    subelements.insert(std::move(*el));
    arma::mat ideal = activeEl.getNodesMatrix();
    idealElements[currel] = ideal;
    currel->setIsBL(true);
    
  }
  
}
void BoundaryLayerGenerator::
ReplaceNeighboringElements(element_set& elements,
			   const int Nlayers,
			   const int surf_to_restrict){

  
  //element_set& elements = mesh.getElementsNC();
  const node_map& nodes = mesh.getNodes();
  //ideal_map& idealElements = mesh.getIdealElementsNC();
  auto& idealElements = ideal_elements;


  const gind Nbln = newnode_map.size();
  
  std::vector<unique_element_ptr> newels;
  for(auto el = elements.begin(); el != elements.end();){
    const int bc_tag = (*el)->getBCTag();
    if(bc_tag == surf_to_restrict || surf_to_restrict == -1){

    const gind* cn = (*el)->getCornerNodes();
    bool found=false;
    //gind nds[(*el)->numCornerNodes()];
    std::vector<gind> nds((*el)->numCornerNodes());
    for(int i=0; i<(*el)->numCornerNodes(); i++){
      auto it = newnode_map.find(cn[i]);
      
      if(it != newnode_map.end()){
	nds[i] = it->second+(Nlayers-1)*Nbln;
	found = true;
      }
      else{
	nds[i] = cn[i];
      }
      
    }
    if(found){
      int type = (*el)->getElementType();

      // Save ideal shape
      const gind* cn = (*el)->getCornerNodes();
      const int ncn = (*el)->numCornerNodes();
      arma::mat ideal(3,ncn);
      for(int i=0; i<ncn; i++){
	ideal.col(i) = nodes.at(cn[i])->getXYZ();
      }
      const ShapeFunctionMatrices* sf = 
	sf_factory.getShapeFunction((*el)->getElementType(),
				    (*el)->getOrder(),0);
      
      
      int geo_entity = (*el)->getGeoEntity();
      int bc_tag = (*el)->getBCTag();
      bool geo_type = (*el)->getGeoType();


      el = elements.erase(el);
      newels.push_back(elementFactory::Instance()->CreateElement
		       (type,nds,0,1,bc_tag,geo_entity,geo_type));
      idealElements[newels.back().get()] = ideal;
    }
    else{
      ++el;
    }
    }
    else{
      ++el;
    }

  }
  
  for(auto el=newels.begin(); el != newels.end(); ++el){
    elements.insert(std::move(*el));
  }

}

void BoundaryLayerGenerator::
ResetBLNodes(const double ratio, 
	      const std::vector<double>& extrusion_dist){

}

void BoundaryLayerGenerator::
InsertBLNodes(const double ratio, 
	      const std::vector<double>& extrusion_dist){

  node_map& nodes = mesh.getNodesNC();

  int nodecnt = 0;
  for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd,++nodecnt){
    //std::cout << nodecnt << std::endl;
    MNode* currnode = nodes.at(nd->first).get();
    const int dim = nd->second.node_dim;

    const double h = extrusion_dist[nodecnt];

    
    arma::vec3 xyz = currnode->getXYZ();
    arma::vec3 coord = xyz - ratio*h*nd->second.normal;
    const gind nc = nodeFactory::Instance()->GetNodeCount();

    int entity = -1;

    arma::vec2 params;
    arma::vec3 newpts = coord;

    //std::cout << newpts << std::endl;

    if(mesh.Dimension() < 3) 
      entity = NodeParamsOnGeometry(currnode,nd->second.normal,
				    nd->second.dist,params,xyz);

    //std::cout <<  "h0" << std::endl;
    if(nd->second.on_symmetry){

      entity = nd->second.symmetry_entity;

      const myFace& face = geometry.getFaces()[entity];
	    
      //const arma::vec3 xyznew = currnode->xyzvec3()-
      // 	nd->second.dist*nd->second.normal;
      params = 
	face.faceParamsFromPoint(newpts.memptr(),
				 currnode->getParametricCoords());

      face.param2xyz(params.memptr(),newpts.memptr());

      if(ratio == 1){
	nodes[nc] = nodeFactory::Instance()->
	  CreateNode(2,newpts.memptr(),nd->second.symmetry_entity,
		     params[0],params[1]);
	    
      }
      //std::cout << "end nd->onSymmetry" << std::endl;
    }
    else if(ratio == 1){
      //std::cout << "begin ratio == 1" << std::endl;
      int type = mesh.MeshDimension();
      nodes[nc] = nodeFactory::Instance()->
	CreateNode(type,newpts.memptr(),entity,params[0],params[1]);
      //std::cout << "end ratio == 1" << std::endl;
    }

    //std::cout << "before ratio != 1" << std::endl;
    if(ratio != 1){
      auto newnd = nodes.at(newnode_map[nd->first]);
      //auto& nd = nodes.at(NNorig+nodecnt);
      //std::cout << "setting newpts" << std::endl;
      //std::cout << newpts << std::endl;
      newnd->setXYZ(newpts.memptr());
      if(newnd->getType() != 3) newnd->setParametricCoords(params.memptr());
      else newnd->setParametricCoords(newpts.memptr());
      //newnode_map[nd->first]-NNorig
      //nodes[newnode_map[nd->first]]->setXYZ(xyz.memptr());
      //nodes[newnode_map[nd->first]]->setParametricCoords(params.memptr());

    }
  }



}

/*
void BoundaryLayerGenerator::
InsertBLNodes(const normal_map_type& normal_map,
	      const int Nlayers, const double thickness, const double ratio,
	      const std::vector<double>& max_extrude){

  node_map& nodes = mesh.getNodesNC();

  const double h = double(thickness/Nlayers);
  for(int lay=0; lay<Nlayers; lay++){
    for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
      MNode* currnode = nodes.at(nd->first).get();
      const int dim = nd->second.node_dim;

      
      arma::vec3 xyz = currnode->getXYZ();
      arma::vec3 coord = xyz - double(lay+1)*h*nd->second.normal;
      const gind nc = nodeFactory::Instance()->GetNodeCount();
      int type = mesh.MeshDimension();
      if(nd->second.on_symmetry){
	type = nd->second.node_dim;
	type = 2; // maybe a fix
      }
      
      arma::vec params;
  
      //cout << "par[0]: " << par[0] << endl;
      if(nd->second.on_symmetry || (mesh.MeshDimension() == 2 && dim < 2)){
	arma::vec2 params;

	if(ratio == -1){
	  

	  //node_reparameterizer.ReParametrizeOnFace(entity,
	  short int entity = NodeParamsOnGeometry(currnode,nd->second.normal,
					      nd->second.dist,params,xyz);

	  if(nd->second.on_symmetry){
	    entity = nd->second.symmetry_entity;
	  }

	  if(nd->second.on_symmetry){
	    const myFace& face = geometry.getFaces()[entity];
	    
	    const arma::vec3 xyznew = currnode->xyzvec3()-
	      nd->second.dist*nd->second.normal;
	    arma::vec2 newparams = 
	      face.faceParamsFromPoint(xyznew.memptr(),
				       currnode->getParametricCoords());
	    arma::vec3 newpts;
	    face.param2xyz(newparams.memptr(),newpts.memptr());

	    nodes[nc] = nodeFactory::Instance()->
	      CreateNode(2,newpts.memptr(),nd->second.symmetry_entity,
			 newparams[0],newparams[1]);
	    
	  }
	  else{
	    nodes[nc] = nodeFactory::Instance()->
	      CreateNode(type,xyz.memptr(),entity,params[0],params[1]);
	  }
	}
	else{
	  //std::cout << "in here" << std::endl;

	  if(nd->second.on_symmetry){
	    const myFace& face = 
	      geometry.getFaces()[nd->second.symmetry_entity];
	 
	    assert(newnode_map[nd->first]-NNorig < max_extrude.size());

	    double dist = max_extrude[newnode_map[nd->first]-NNorig];
	    //std::cout << dist << std::endl;
	    const arma::vec3 xyznew = currnode->xyzvec3()-
	      dist*nd->second.normal;
   
	    //const arma::vec3 xyznew = currnode->xyzvec3()-
	    //  ratio*nd->second.dist*nd->second.normal;

	    arma::vec2 newparams = 
	      face.faceParamsFromPoint(xyznew.memptr(),
	    			       currnode->getParametricCoords());
	    arma::vec3 newpts;
	    face.param2xyz(newparams.memptr(),newpts.memptr());

	    nodes[newnode_map[nd->first]]->setXYZ(newpts.memptr());
	    nodes[newnode_map[nd->first]]->
	      setParametricCoords(newparams.memptr());
	  }
	  else{
	    double dist = max_extrude[newnode_map[nd->first]-NNorig];
	    short int entity = 
	      NodeParamsOnGeometry(currnode,nd->second.normal,
	    			   dist,params,xyz);

	    
	    //short int entity = 
	    //  NodeParamsOnGeometry(currnode,nd->second.normal,
	    //			   ratio*nd->second.dist,params,xyz);
	    

	    nodes[newnode_map[nd->first]]->setXYZ(xyz.memptr());
	    nodes[newnode_map[nd->first]]->setParametricCoords(params.memptr());
	  }
	}
      }
      else{
	if(ratio == -1){
	  //std::cout << arma::norm(coord) << std::endl;
	  //if(dim == 1) std::cout << "For some reason in here" << std::endl;
	  nodes[nc] = nodeFactory::Instance()->
	    CreateNode(type,coord.memptr(),-1);
	}
	else{
	  arma::vec3 temp = currnode->getXYZ() - 
	    ratio*double(lay+1)*h*nd->second.normal;
	    nodes[newnode_map[nd->first]]->setXYZ(temp.memptr());
	}
      }
    }
  }
}
*/

/*
void BoundaryLayerGenerator::
ResetBLNodes(const normal_map_type& normal_map,
	     const std::map<gind,gind>& newnode_map,
	     const int Nlayers, const double thickness, const double factor){

  node_map& nodes = mesh.getNodesNC();

  const int Nbln = normal_map.size();
  const double h = double(thickness/Nlayers);
  for(int lay=0; lay<Nlayers; lay++){
    for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
      arma::vec3 xyz = nodes.at(nd->first)->getXYZ();
      arma::vec3 xyz2 = nodes.at(newnode_map.at(nd->first))->getXYZ();
      arma::vec3 newxyz = xyz + factor*(xyz2-xyz);
      nodes.at(newnode_map.at(nd->first))->setXYZ(newxyz.memptr());
    }
    
  }

}
*/

void BoundaryLayerGenerator::
InsertBLElements(element_set& elements,
		 const std::vector<MEl*>& bl_elements,
		 const int Nlayers,
		 const double thickness,
		 const bool copy_geo_entity){

  //OptElManager optel_manager(mesh.getNodesNC(),sf_factory,index_factory);

  //element_set& elements = mesh.getElementsNC();
  
  //ideal_map& idealElements = mesh.getIdealElementsNC();
  auto& idealElements = ideal_elements;

  node_map& nodes = mesh.getNodesNC();

  const arma::vec3 n_planar = {0.0, 0.0, 1.0};

  const int Nbln = newnode_map.size();

  typedef std::pair<element_set::iterator,bool> retpair;

  for(int lay=0; lay< Nlayers; lay++){
    for(auto el = bl_elements.begin(); el != bl_elements.end(); ++el){
      //gind nds[8];
      const int type = (*el)->getElementType();
      const int dim = (*el)->getDim();
      const int ncn = (*el)->numCornerNodes();
      const gind* cn = (*el)->getCornerNodes();

      std::vector<gind> nds(2*ncn);
      arma::vec3 elnorm;
      if(dim == 1){
	elnorm = El1DNormal(cn,nodes,n_planar);
      }
      else{
	elnorm = El2DNormal(cn,nodes);
      }
      arma::mat ideal(3,2*ncn);
      arma::mat idealExtrude(3,2*ncn);
      for(int i=0; i<ncn; i++){
	int rn;
	if(type == 1) rn = i;
	else rn = ncn-i-1;
	if(lay==0) nds[i] = cn[rn];
	else nds[i] = newnode_map.at(cn[rn]) + (lay-1)*Nbln;
	nds[i+ncn] = newnode_map.at(cn[rn]) + lay*Nbln;
	ideal.col(i) = nodes.at(nds[i])->getXYZ();
	idealExtrude.col(i) = nodes.at(nds[i])->getXYZ();
	idealExtrude.col(i+ncn) = ideal.col(i) - thickness*elnorm;
	//ideal.col(i+ncn) = ideal.col(i) - thickness*elnorm;

	ideal.col(i+ncn) = nodes.at(nds[i+ncn])->getXYZ();
      }
      if(type == 1){
	//std::cout << "ideal: " << std::endl;
	//std::cout << ideal << std::endl;
	//ideal.save("ideal",arma::raw_ascii);
      }
      int entity=-1;
      if(copy_geo_entity) entity = (*el)->getGeoEntity();
      //if((*el)->getDim() == 2) entity = 0; 
      //if(type == 1) entity = 0; // Careful, this is wrong!
      if(mesh.Dimension() == 2) entity = 0;

      retpair ret;
      if(type == 1){ // Line to quad
	int bc_tag = 0;
	if(mesh.Dimension() == 3){
	  bc_tag = (*el)->getBCTag();
	  //entity = mesh.SymmetrySurface();
	  entity = (*el)->getGeoEntity();
	  //std::cout << "symmetry ent: " << entity << std::endl;
	}
	//int entity = 9;
	ret = elements.insert(elementFactory::Instance()->
			      CreateElement(3,nds,0,1,bc_tag,entity,1)); // entity is a hack!
      }
      else if(type == 2){ // Tri to Prism
	ret = elements.insert(elementFactory::Instance()->
			      CreateElement(6,nds,0,1,0,-1,0));
      } 

      ret.first->get()->setIsBL(true);

      assert(ideal.n_rows == idealExtrude.n_rows);
      assert(ideal.n_cols == idealExtrude.n_cols);
      //std::cout << idealExtrude << std::endl;

      idealElements[(*ret.first).get()] = ideal;

      //if(dim == mesh.Dimension()-2) idealElements[(*ret.first).get()] = ideal;
      //else idealElements[(*ret.first).get()] = idealExtrude;

      /*
      const MEl* newel = *ret.first;
      
      unique_element_ptr ideal_ = 
	elementFactory::Instance()->CreateElement(newel->getElementType(),
						  ideal_nodes,0,1,0,
						  newel->getGeoEntity(),
						  newel->getGeoType());
      ideal_map[*ret.first] = ideal_;
      */

      /*
      const ShapeFunctionMatrices* sf = 
	sf_factory.getShapeFunction((*ret.first)->getElementType(),
				    (*ret.first)->getOrder(),0);

      
      if(type == 1){
	//cout << "idealExtrude: " << endl;
	//cout << idealExtrude << endl;
      }
       std::unique_ptr<OptEl> optel;
       if((*el)->getDim() == 2){
	optel = std::unique_ptr<OptEl3D>
	  (new OptEl3D(ret.first->get(),index_factory,nodes,sf,ideal));
       }
       else if((*el)->getDim() == 1){
	optel = std::unique_ptr<OptEl2D>
	  (new OptEl2D(ret.first->get(),index_factory,nodes,sf,ideal));
       }
       arma::mat gradMerit;
       double dist;
       double detS;
       optel->computeGradMerit(gradMerit,dist,detS);
       //double dist = optel->computeDistortion();
      
      if(dist > 100.0){
	//ExtrudedIdeal[(*ret.first).get()] = true;
	idealElements[(*ret.first).get()] = idealExtrude;
	if(type == 1 && mesh.MeshDimension() == 3){
	  idealElements[(*ret.first).get()] = ideal;
	}
      }
      else{
	idealElements[(*ret.first).get()] = ideal;
	//ExtrudedIdeal[(*ret.first).get()] = false;
      }
      
      */
    }
  }
}

std::vector<double> BoundaryLayerGenerator::
getMaximumSafeExtrusionDistance(int safety_factor){
  const element_set& elements = mesh.getElements();
  //std::cout << "sizeof set: " << sizeof(std::set<int>) << std::endl;

  std::vector<std::set<int> > bl_node_neighbors(normal_map.size());

  for(auto el = elements.begin(); el != elements.end(); ++el){
    auto element = el->get();
    const gind* elnodes = element->getNodes();
    const int nn = element->NumNodes();
    for(int nd = 0; nd < nn; nd++){
      int node_id = elnodes[nd];
      auto it = normal_map.find(node_id);
      //auto itbl = bl_node_neighbors.find(node_i
      if(it != normal_map.end()){
	
	const int nn = element->NumNodes();
	for(int nd = 0; nd < nn; nd++){
	  int node_id_j = elnodes[nd];
	  bl_node_neighbors[newnode_map[node_id]-NNorig].insert(node_id_j);
	}
      }    
    }
  }

  const node_map& nodes = mesh.getNodes();

  std::vector<double> max_dist(normal_map.size(),1.0/0.0);

  //int bln = 0;
  for(auto nr = normal_map.begin(); nr != normal_map.end(); ++nr){
    int bln = newnode_map[nr->first] - NNorig;
    double min_dist = 1.0/0.0;
    for(auto it = bl_node_neighbors[bln].begin(); 
	it !=bl_node_neighbors[bln].end(); ++it){
 
      if(newnode_map.find(*it) == newnode_map.end()){
	const arma::vec3& r = nodes.at(*it)->xyzvec3()-
	  nodes.at(nr->first)->xyzvec3();
	double r_dot_n = arma::dot(r,-nr->second.normal);
	//double dist = arma::norm(r,2);
	double dist = std::abs(r_dot_n);

	min_dist = std::min(min_dist,dist);
      }
    }
    max_dist[bln] = 1.0/safety_factor*min_dist;
  }


  return max_dist;
}

int BoundaryLayerGenerator::SubDivideBL(){
 BoundaryLayerSubdivider bl_subdivider(mesh,bl_parameters);
 bl_subdivider.GenerateElements();
  
}
