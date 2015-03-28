#include "MeshReader.h"
#include "MeshContainer.h"
#include "GeometryContainer.h"
#include "ChildGenerator.h"

#include "Mel.h"
#include "Mnode.h"

#include "elementFactory.h"
#include "nodeFactory.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"

#include <iostream>
#include <fstream>
#include <algorithm>



using namespace std;

bool checkMeshHeader(ifstream &fname, const string &header){
  string line;
  getline(fname,line);
  if(line == header) return true;
  else return false;

}

int writePrismElements(MeshContainer& mesh){
  const element_set& elements = mesh.getElements();
  const node_map& nodes_map = mesh.getNodes();




  int elcnt = 0;
  int pyrcnt = 0;
  for(auto it = elements.begin(); it != elements.end(); ++it){
    const MEl* el = it->get();
    if(elcnt == 17493){
      arma::mat elnodes(3,5);
      const gind* nds = el->getNodes();
      for(int n = 0; n < el->NumNodes(); n++){
	elnodes.col(n) = nodes_map.at(nds[n])->xyzvec3();
      }
      elnodes.save("nodes"+std::to_string(elcnt),arma::raw_ascii);
      pyrcnt++;
    }
    ++elcnt;
  }
}

int SplitPyramid2Tet(MeshContainer& mesh){
  element_set& elements = mesh.getElementsNC();
  element_set& subelements = mesh.getSubElementsNC();

  const node_map& nodes = mesh.getNodes();

  NodeIndexFactory IndexFactory;

  ChildGenerator child_generater;

  int elcnt = 0;
  for(auto it = elements.begin(); it != elements.end();++elcnt){
    const MEl* el = it->get();
    const int ncn = el->numCornerNodes();
    const int order = el->getOrder();
    const int type = el->getElementType();
    const int dim = el->getDim();

    if(type == 7){
      const NodeIndexer* indexer = 
	IndexFactory.getNodeIndexer(type,order);

      const int geoEntity = el->getGeoEntity();
      const int bcTag = el->getBCTag();
      const bool geoType = el->getGeoType();

      const NodeIndexer* tet_indexer = IndexFactory.getNodeIndexer(4,order);

      const gind* pyr_nodes = el->getNodes();
      /*
      std::cout << "pyramid nodes: " << std::endl;
      for(int i = 0; i < el->NumNodes(); i++){
	std::cout << pyr_nodes[i] << " ";
      }
      std::cout << std::endl;
      */

      std::vector<int> pyramid_nodes(el->NumNodes());
      arma::mat pyrnodes(el->NumNodes(),3);
      for(int i = 0; i < el->NumNodes(); i++){
	pyramid_nodes[i] = pyr_nodes[i];
	for(int j = 0; j < 3; j++){
	  pyrnodes(i,j) = nodes.at(pyr_nodes[i])->xyzptr()[j];
	}
      }
      //pyrnodes.save("pyramid_nodes"+std::to_string(elcnt),arma::raw_ascii);

      const int ndof_tet = tet_indexer->Ndof();
      std::vector<int> tet1_nodes(ndof_tet), tet2_nodes(ndof_tet);
      
      int pyr_base_cnt = 0;
      for(int k = 0, cnt = 0; k <=order; k++){
	for(int j = 0; j <=order-k; j++){
	  for(int i = 0; i <= order-j-k; i++, cnt++){
	    int pyr_ind1 = i + j*(order-k+1)  + pyr_base_cnt;
	    tet1_nodes[cnt] = pyr_nodes[pyr_ind1];
	    //tet1_nodes[cnt] = pyr_ind1;

	    //int pyr_ind2 = i +(order-k) + (order-k-j+1)*(order-k+1)  + pyr_base_cnt;
	    int pyr_ind2 = order-k-i + (order-k-j)*(order-k+1)  + pyr_base_cnt;
	    tet2_nodes[cnt] = pyr_nodes[pyr_ind2];
	    //tet2_nodes[cnt] = pyr_ind2;

	  }
	}
	pyr_base_cnt+= pow(order-k+1,2);
       
      }
      std::vector<int> added_pyr_line = {pyr_nodes[1], pyr_nodes[2]};

      /*
      std::cout << "tet 1/2 nodes: " << std::endl;
      for(int i = 0; i < ndof_tet; i++){
	std::cout << tet1_nodes[i] << " ";
      }
      std::cout << std::endl;
      for(int i = 0; i < ndof_tet; i++){
	std::cout << tet2_nodes[i] << " ";
      }
      std::cout << std::endl;
      */
      

      std::unique_ptr<MEl> tet1 = 
	elementFactory::Instance()->CreateElement
	(4,tet1_nodes,tet_indexer,order,bcTag,geoEntity,geoType);

      std::unique_ptr<MEl> tet2 = 
	elementFactory::Instance()->CreateElement
	(4,tet2_nodes,tet_indexer,order,bcTag,geoEntity,geoType);

      arma::mat tet1nodes(tet1->NumNodes(),3);
      for(int i = 0; i < tet1->NumNodes(); i++){
	for(int j = 0; j < 3; j++){
	  tet1nodes(i,j) = nodes.at(tet1_nodes[i])->xyzptr()[j];
	}
      }
      //tet1nodes.save("tet1_nodes"+std::to_string(elcnt)+"_1",arma::raw_ascii);

      arma::mat tet2nodes(tet1->NumNodes(),3);
      for(int i = 0; i < tet1->NumNodes(); i++){
	for(int j = 0; j < 3; j++){
	  tet2nodes(i,j) = nodes.at(tet2_nodes[i])->xyzptr()[j];
	}
      }
      //tet2nodes.save("tet2_nodes"+std::to_string(elcnt)+"_2",arma::raw_ascii);

      std::unique_ptr<MEl> newpyr = 
	elementFactory::Instance()->CreateElement
	(type,pyramid_nodes,indexer,order,bcTag,geoEntity,geoType);

      const std::vector<element_ptr>& children = 
	child_generater.GenerateChildren(*it);
      

      for(auto ch = children.begin(); ch != children.end(); ++ch){
	const element_ptr& child = *ch;
	auto se = subelements.find(child);
	if(se != subelements.end()){
	  const element_ptr& subel = *se;
	  const int type = subel->getElementType();
	  assert(type == 3);
	  
	  if(type == 3){
	    //std::cout << "replacing quad with tri!" << std::endl;
	    const NodeIndexer* tri_indexer = 
	      IndexFactory.getNodeIndexer(2,order);

	    const int geoEntity = subel->getGeoEntity();
	    const int bcTag = subel->getBCTag();
	    const bool geoType = subel->getGeoType();

	    const gind* quad_nodes = subel->getNodes();
	    
	    const std::vector<indtype>& reversed = tri_indexer->getReversedNodes();
	    
	    //std::cout << "quad nodes: " << std::endl;
	    //for(int i = 0; i < subel->NumNodes(); i++){
	    //  std::cout << quad_nodes[i] << " ";
	    //}
	    //std::cout << std::endl;
	    

	    const int tri_dof = tri_indexer->Ndof();
	    std::vector<int> tri1_nodes(tri_dof), tri2_nodes(tri_dof);
	    for(int j = 0, cnt = 0; j <=order; j++){
	      for(int i = 0; i <= order-j; i++, cnt++){
		int quad_ind1 = i + j*(order+1);
		//tri1_nodes[cnt] = quad_nodes[quad_ind1];
		tri1_nodes[reversed[cnt]] = pyr_nodes[quad_ind1];
		//tet1_nodes[cnt] = pyr_ind1;

		int quad_ind2 = order-i + (order-j)*(order+1);
		//tri2_nodes[cnt] = quad_nodes[quad_ind2];
		tri2_nodes[reversed[cnt]] = pyr_nodes[quad_ind2];
		//tet2_nodes[cnt] = pyr_ind2;

	      }
	    }
	    std::vector<int> added_quad_line = {quad_nodes[1],quad_nodes[2]};
	    /*
	    std::cout << "tri 1/2 nodes: " << std::endl;
	    for(int i = 0; i < tri_dof; i++){
	      std::cout << tri1_nodes[i] << " ";
	    }
	    std::cout << std::endl;
	    for(int i = 0; i < tri_dof; i++){
	      std::cout << tri2_nodes[i] << " ";
	    }
	    std::cout << std::endl;
	    */

	    //std::sort(added_quad_line.begin(),added_quad_line.end());
	    //std::sort(added_pyr_line.begin(),added_pyr_line.end());
	    //assert(std::equal(added_quad_line.begin(),
	    //	      added_quad_line.begin()+order+1,
	    //		      added_pyr_line.begin()));
	    

	    std::unique_ptr<MEl> tri1 = 
	      elementFactory::Instance()->CreateElement
	      (2,tri1_nodes,tri_indexer,order,bcTag,geoEntity,geoType);

	    std::unique_ptr<MEl> tri2 = 
	      elementFactory::Instance()->CreateElement
	      (2,tri2_nodes,tri_indexer,order,bcTag,geoEntity,geoType);

	    se = subelements.erase(se);
	    se = subelements.insert(se,std::move(tri1));
	    se = subelements.insert(se,std::move(tri2));
	  }
	}
      }
      it = elements.erase(it);
      //it = elements.insert(it,std::move(newpyr));
      it = elements.insert(it,std::move(tet1));
      it = elements.insert(it,std::move(tet2));

      //el = 
    }
    else{
      ++it;
    }
  }
}

int CheckElementJacobianValidity(MeshContainer& mesh){
  element_set& elements = mesh.getElementsNC();
  const node_map& nodes = mesh.getNodes();


  ShapeFunctionMatricesFactory sf_factory;


  double minDetJ_all = 1.0/0.0;
  double maxDetJ_all = -1.0/0.0;

  int elcnt = 0;
  for(auto it = elements.begin(); it != elements.end(); ++it, ++elcnt){
    const MEl* el = it->get();
    const int ncn = el->numCornerNodes();
    const int order = el->getOrder();
    const int type = el->getElementType();
    const int dim = el->getDim();

    if(type == 7) std::cout << "Still have prism elements!" << std::endl;


    const ShapeFunctionMatrices& sf = 
      *sf_factory.getShapeFunction(type,order,0);



    const int nn = el->NumNodes();
    const gind* elnodes = el->getNodes();

    arma::mat nds(nn,3);
    for(int n = 0; n < nn; n++){
      for(int i = 0; i < 3; i++){
	nds(n,i) = nodes.at(elnodes[n])->xyzptr()[i];
      }
    }

    const int ng = sf.getQuadratureSF().n_cols;


    const arma::cube& sfderiv = sf.getQuadratureSFDeriv();

    
    const arma::mat& nodes_xi = sfderiv.slice(0).t()*nds;
    const arma::mat& nodes_eta = sfderiv.slice(1).t()*nds;


    arma::vec detJ;
    if(dim == 2){

      detJ = 
	nodes_xi.unsafe_col(0)%nodes_eta.unsafe_col(1) -
	nodes_xi.unsafe_col(1)%nodes_eta.unsafe_col(0);
    }
    else if(dim == 3){
      arma::mat nodes_zeta = sfderiv.slice(2).t()*nds;
      const arma::vec& x_xi = nodes_xi.unsafe_col(0);
      const arma::vec& y_xi = nodes_xi.unsafe_col(1);
      const arma::vec& z_xi = nodes_xi.unsafe_col(2); 

      const arma::vec& x_eta = nodes_eta.unsafe_col(0);
      const arma::vec& y_eta = nodes_eta.unsafe_col(1);
      const arma::vec& z_eta = nodes_eta.unsafe_col(2);     

      const arma::vec& x_zeta = nodes_zeta.unsafe_col(0);
      const arma::vec& y_zeta = nodes_zeta.unsafe_col(1);
      const arma::vec& z_zeta = nodes_zeta.unsafe_col(2); 
      detJ = 
	x_eta%y_zeta%z_xi - 
	x_eta%y_xi%z_zeta + 
	x_xi%y_eta%z_zeta - 
	x_xi%y_zeta%z_eta - 
	x_zeta%y_eta%z_xi + 
	x_zeta%y_xi%z_eta;
    }

    double minDetJ = detJ.min();
    double maxDetJ = detJ.max();
    if(minDetJ < minDetJ_all) minDetJ_all = minDetJ;
    if(maxDetJ > maxDetJ_all) maxDetJ_all = maxDetJ;
    //std::cout << minDetJ << " " << maxDetJ << std::endl;

    if(minDetJ < 0.0){
      std::cout << "min(DetJ) is " << minDetJ << " on element " << elcnt << 
	" of type " << type << std::endl;
      std::cout << nds << std::endl;
      
    }
  }
  std::cout << "Min/Max det(J): " << minDetJ_all << " " << maxDetJ_all << 
    std::endl;
}

int fixPyramidOrientation(MeshContainer& mesh){
  element_set& elements = mesh.getElementsNC();
  const node_map& nodes = mesh.getNodes();

  NodeIndexFactory IndexFactory;

  for(auto it = elements.begin(); it != elements.end(); ){
    const MEl* el = it->get();
    const int ncn = el->numCornerNodes();
    const int order = el->getOrder();
    const int type = el->getElementType();

    if(type == 7){
      const gind* elnodes = el->getNodes();
      const arma::vec3& n0 = nodes.at(elnodes[0])->xyzvec3();
      const arma::vec3& n1 = nodes.at(elnodes[1])->xyzvec3();
      const arma::vec3& n2 = nodes.at(elnodes[2])->xyzvec3();
      const arma::vec3& n3 = nodes.at(elnodes[3])->xyzvec3();
      const arma::vec3& n4 = nodes.at(elnodes[4])->xyzvec3();

      for(int i = 0; i < el->NumNodes(); i++){
	const arma::vec3& temp = nodes.at(elnodes[i])->xyzvec3();
	for(int j = 0; j < 3; j++){
	  std::cout << std::setprecision(16) << std::setw(20) << temp[j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;

      /*
      std::cout << n0.t() << std::endl;
      std::cout << n1.t() << std::endl;
      std::cout << n2.t() << std::endl;
      std::cout << n3.t() << std::endl;
      std::cout << n4.t() << std::endl;
      std::cout << std::endl;
      */
      const arma::vec3& r1 = n1-n0;
      const arma::vec3& r2 = n2-n0;
      const arma::vec3& normal = arma::cross(r1,r2);
      const arma::vec3& r4 = n4-n0;

      if(arma::dot(r4,normal) < 0.0){
	std::cout << "Prism with tip vertex is inverted!" << std::endl;
	std::cout << n4.t() << std::endl;

	const int geoEntity = el->getGeoEntity();
	const int bcTag = el->getBCTag();
	const bool geoType = el->getGeoType();
      
	const NodeIndexer* indexer = 
	  IndexFactory.getNodeIndexer(type,order);

	std::vector<int> nodesvec(el->NumNodes());
	nodesvec[0] = elnodes[2];
	nodesvec[1] = elnodes[3];
	nodesvec[2] = elnodes[0];
	nodesvec[3] = elnodes[1];
	nodesvec[4] = elnodes[4];

	std::unique_ptr<MEl> newel = 
	  elementFactory::Instance()->CreateElement
	  (type,nodesvec,indexer,order,bcTag,geoEntity,geoType);
      
	//it = elements.erase(it);
      
	//elements.insert(std::move(newel));
	++it;

      }
      else{
	++it;
      }
    }
    else{
      ++it;
    }
  }

};

int fix2DNormals(MeshContainer& mesh){
  std::cout << "in fix 2D normals" << std::endl;

  element_set& elements = mesh.getElementsNC();
  const node_map& nodes_map = mesh.getNodes();
  arma::vec3 nz = {0.0, 0.0, 1.0};

  element_set newels;
  
  NodeIndexFactory IndexFactory;

  for(auto it = elements.begin(); it != elements.end(); ){
    
    //const int dim = el->getDim();
    const MEl* el = it->get();
    const int ncn = el->numCornerNodes();
    const int order = el->getOrder();
    const int type = el->getElementType();

    int cn[ncn], cng[ncn];
    el->getCornerNodes(cn);
    const gind* nodes = el->getNodes();
    for(int i = 0; i < ncn; i++){
      cng[i] = nodes[cn[i]];
    }
    const arma::vec3& n = normalFromElnodes(cng,nodes_map);
    if(arma::dot(n,nz) < 0.0){
      const NodeIndexer* indexer = 
	IndexFactory.getNodeIndexer(type,order);

      const int geoEntity = el->getGeoEntity();
      const int bcTag = el->getBCTag();
      const bool geoType = el->getGeoType();
      
      const std::vector<indtype>& reversed = indexer->getReversedNodes();
      std::vector<int> nodesvec(el->NumNodes());
      for(int i = 0; i < el->NumNodes(); i++){
	nodesvec[i] = nodes[reversed[i]];
      }
      std::unique_ptr<MEl> newel = 
	elementFactory::Instance()->CreateElement
	(type,nodesvec,indexer,order,bcTag,geoEntity,geoType);
      
      it = elements.erase(it);
      
      elements.insert(std::move(newel));
      
    }
    else{
      ++it;
    }
    
  }


}

int readGMSHNodes(MeshContainer& mesh, ifstream &fname){
  node_map& mesh_nodes = mesh.getNodesNC();

  if(checkMeshHeader(fname,"$ParametricNodes")){
  
    int nnodes;
    fname >> nnodes;
    fname.ignore(1,'\n');
    

    for(int nd=0; nd<nnodes; nd++){
      double xyz[3];
      int nodeid, dim, gmshID;

      fname >> nodeid;
      nodeid--;

      for(int i=0; i<3; i++) fname >> xyz[i]; 
      fname >> dim >> gmshID;
      if(dim > 3){
	throw std::runtime_error("Node dimension cannot be larger than 3.");
      }
      double uv[2];
      for(int i=0; i<dim; i++){
	if(dim != 3) fname >> uv[i];
      }
      //cout << "dim: " << dim << endl;
      //gind ndcurr = nodeFactory::Instance()->GetNodeCount();
      //std::cout << gmshID << std::endl;
      mesh_nodes[nodeid] = 
	nodeFactory::Instance()->CreateNode(dim,xyz,gmshID-1,uv[0],uv[1]);
      //std::cout << mesh_nodes[nodeid]->getGeoEntity() << std::endl;
      //cout << "type: " << mesh_nodes[ndcurr]->getType() << endl;
    }
    fname.ignore(1,'\n');


    if(!checkMeshHeader(fname,"$EndParametricNodes")){
      throw std::runtime_error("Could not find matching $EndParametricNodes!");
    }
  }
  else{
    std::string error ="$ParametricNodes must be followed by $ParametricNodes!";
    throw std::runtime_error(error);
  }
  std::cout << "After reading nodes" << std::endl;

  return true;
}

int readGMSHElements(MeshContainer& mesh, ifstream &fname,
		     std::map<int,int>& edge2bc_map,
		     std::map<int,int>& face2bc_map){

  //element_map& mesh_elements = mesh.getElementsNC();
  //element_map& subelements = mesh.getSubElementsNC();
  //element_map& subsubelements = mesh.getSubSubElementsNC();

  //inverse_element_map& inv_elements = mesh.getInverseElementsNC();


  element_set& mesh_elements = mesh.getElementsNC();
  element_set& subelements = mesh.getSubElementsNC();
  element_set& subsubelements = mesh.getSubSubElementsNC();

  NodeIndexFactory IndexFactory;

  arma::wall_clock timer;
  timer.tic();
  std::vector<std::unique_ptr<MEl> > element_vector;

  if(checkMeshHeader(fname,"$Elements")){
  
    //int elcnt=0;
    int nelements;
    fname >> nelements;
    fname.ignore(1,'\n');


    bool my_geo_type[16] = {false};
    //my_geo_type[1] = 0;
    my_geo_type[2] = 1;
    my_geo_type[3] = 1;



    for(int el=0; el<nelements; el++){
      int temp, eltype, ntags;
      fname >> temp >> eltype >> ntags;

      while(fname.get() != '\n'){
	fname.unget();

	int tags[10];
	for(int i=0; i<ntags; i++) fname >> tags[i];

	//for(int i=0; i<ntags; i++) tags[i]-=1;
	tags[1]-= 1;
	//if(el_dim[eltype] > 2) tags[1] = -1;
	if(tags[0] == 20) mesh.SetSymmetrySurface(tags[1]);



	const int mytype = IndexFactory.TypeFromGMSHType(eltype);
	const int myorder = IndexFactory.OrderFromGMSHType(eltype);
	
	
	//std::cout << mytype << " " << myorder << std::endl;
	//std::cout << el/double(nelements) << std::endl;
	//std::cout << "h0" << std::endl;
	const NodeIndexer* indexer = 
	  IndexFactory.getNodeIndexer(mytype,myorder);

	//std::cout << "h1" << std::endl;

	const std::vector<indtype>& gmind = indexer->getGMSHIndex();
	
	//std::cout << eltype << " " << mytype << " " << myorder << std::endl;
	std::vector<gind> nodes(gmind.size());
	std::vector<gind> temp(gmind.size());
	//std::cout << gmind.size() << std::endl;
	//gind nodes[8];
	for(int i=0; i<gmind.size(); i++){
	  //std::cout << gmind[i] << " ";
	  fname >> nodes[gmind[i]];
	  nodes[gmind[i]]-=1;
	  //fname >> nodes[i];
	  //nodes[i]-=1;
	}

	//std::cout << std::endl;

	if(mytype != 0){
	  short int geo_entity = tags[1];
	  if(indexer->Dimension() == 3) geo_entity = -1;
	  
	  bool insert = true;
	  if(mytype == 1){
	    if(nodes[0] == nodes[1]) insert = false;
	  }
	  if(insert){
	    element_vector.push_back(elementFactory::Instance()->CreateElement
				     (mytype,nodes,indexer,myorder,tags[0],
				      geo_entity,my_geo_type[mytype]));
	  }

	  
	  if(element_vector.back()->getDim() == 1){
	    //edge2bc_map[tags[1]] = tags[0];
	  }
	  else if(element_vector.back()->getDim() == 2){
	    //face2bc_map[tags[1]] = tags[0];
	  }
	  
	  
	}
       
      }
    }
    timer.toc();

    timer.tic();
 

  
    int mesh_dim=0;
    for(auto el=element_vector.begin(); el != element_vector.end(); ++el){
      mesh_dim = std::max(mesh_dim,(*el)->getDim());
    }


    for(auto el=element_vector.begin(); el != element_vector.end(); ++el){
      
      if((*el)->getDim() == mesh_dim){
	//std::cout << (*el)->getElementType() << (*el)->getDim() << std::endl;
	const MEl* curr = el->get();
	const int nn = curr->NumNodes();
	const gind* nds = curr->getNodes();

	mesh_elements.insert(std::move(*el));
      }
      else if((*el)->getDim() == mesh_dim-1){
	subelements.insert(std::move(*el));
      }
      else{
	subsubelements.insert(std::move(*el));
      }    
      
    }

    mesh.SetMeshDimension(mesh_dim);
    mesh.computeIsMeshPlanar();

    /*
    if(mesh.IsMeshPlanar()){ // Replace all surface nodes with free nodes
      
      node_map& nodes = mesh.getNodesNC();
      node_map newnodes;
      cout << "before erasing" << endl;
      for(auto nd = nodes.begin(); nd != nodes.end();){
	if(nd->second->getType() == 2){
	  gind number = nd->second->getND();
	  arma::vec3 xyz = nd->second->getXYZ();
	  newnodes[nd->first] = 
	    std::unique_ptr<MNode>(new MNode(number,xyz.memptr()));
	  nd = nodes.erase(nd);
	}
	else{
	  ++nd;
	}
      }
      cout << "after erasing" << endl;
      for(auto nd = newnodes.begin(); nd != newnodes.end(); ++nd){
        nodes[nd->first] = std::move(nd->second);
      }

      for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){
	if(nd->second->getType() == 2) cout << "ERROR!" << endl;
      }
      cout << "at end of is planar" << endl;
 
    }
    
    cout << "element set insert time: " << timer.toc() << endl;
    */
    if(!checkMeshHeader(fname,"$EndElements")){
      throw std::runtime_error("Could not find matching $EndElements!");
    }
  }
  else{
    throw std::runtime_error("$ParametricNodes must be followed by $Elements!");
  }
    
  return 1;
}

void GMSHReader::ReadMesh(MeshContainer& mesh, std::string filename){
  try{
    ifstream fname;
    fname.open(filename.c_str());

    if(!fname.good()){
      throw std::runtime_error("Could not open gmsh file: "+filename);
    }  
  
    //*************** Read mesh format ****************//
    if(checkMeshHeader(fname,"$MeshFormat")){
      double GMSH_version;
      int GMSH_file_type, GMSH_data_size;
      fname >> GMSH_version >> GMSH_file_type >> GMSH_data_size;
      fname.ignore(1,'\n');
      if(GMSH_version != 2.2){
	std::string error = "GMSH file version is: " + to_string(GMSH_version)+
	  "! The only supported GMSH file version is 2.2!";
	throw std::runtime_error(error);
      }

      if(!checkMeshHeader(fname,"$EndMeshFormat")){
	throw std::runtime_error("Could not find matching $EndMeshFormat!");
      }
    }
    else{
      throw std::
	runtime_error("First line of GMSH file should be $MeshFormat!");
    }

    arma::wall_clock timer;
    timer.tic();
    //*************** Read parametric nodes ****************//
    readGMSHNodes(mesh,fname);


    //****************** Read all elements *******************//
    //readElements(fname,mesh_elements,mesh_nodes,occ_handler,sf_handler);
    readGMSHElements(mesh,fname,edge2bc_map, face2bc_map);
    std::cout << "Time to read GMSH file: " << timer.toc() << std::endl;

    fname.close();

    
    if(mesh.MeshDimension() == 2){
      fix2DNormals(mesh);
    }

    SplitPyramid2Tet(mesh);

    std::cout << "before fix byramid orientation" << std::endl;
    CheckElementJacobianValidity(mesh);
    //fixPyramidOrientation(mesh);
    std::cout << "After fix pyramid orientation" << std::endl;



    //throw std::runtime_error("Quitting now!");

    //writePrismElements(mesh);
    
  }
  catch(std::exception& e) {
    throw;
  }
 

}

void GMSHReader::SetGeoTags(GeometryContainer& geometry){
  std::vector<myEdge>& edges = geometry.getEdgesNC();
  std::vector<myFace>& faces = geometry.getFacesNC();

  for(auto it = edge2bc_map.begin(); it != edge2bc_map.end(); ++it){
    edges[it->first].setBCTag(it->second);
  }

  for(auto it = face2bc_map.begin(); it != face2bc_map.end(); ++it){
    faces[it->first].setBCTag(it->second);
  }
}
