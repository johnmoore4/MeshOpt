//#include "HModel.h"
//#include "Mel.h"
//#include "Mnode.h"
//#include "shapeFunctionHandler.h"

//#include <vector>
//#include <memory>
#include "MeshContainer.h"
#include "ChildGenerator.h"
#include "elementFactory.h"
#include "El1D.h"
#include "El2D.h"
#include "El3D.h"
#include "Mnode.h"
#include <algorithm>  
#include <assert.h>
#include <list>
#include <forward_list>
#include <iostream>
#include <array>
using std::cout;
using std::endl;



std::unique_ptr<MEl> getOrientedEl(const unique_element_ptr& el, 
				   NodeIndexFactory& ni_factory,
				   int orient){


  const int ncn = el->numCornerNodes();
  const int order = el->getOrder();
  const int type = el->getElementType();
  const int dim = el->getDim();

  const int geoEntity = el->getGeoEntity();
  const int bcTag = el->getBCTag();
  const bool geoType = el->getGeoType();

  const NodeIndexer* indexer = ni_factory.getNodeIndexer(type,order);

  const std::vector<indtype>& oriented_nodes = 
    indexer->getOrientedNodes(orient);

  const gind* nodes = el->getNodes();
  const int nn = el->NumNodes();

  std::vector<int> newnodes(nn);
  
  for(int i = 0; i < nn; i++){
    newnodes[i] = nodes[oriented_nodes[i]];
  }

  return elementFactory::Instance()->CreateElement
    (type,newnodes,indexer,order,bcTag,geoEntity,geoType);


}

int fixFaceOrientations(MeshContainer& mesh){
  ChildGenerator child_generator;
  NodeIndexFactory index_factory;

  const element_set& elements = mesh.getElements();

  element_set& subelements = mesh.getSubElementsNC();
  
  for(auto el = elements.begin(); el != elements.end(); ++el){
    std::vector<unique_element_ptr>& childels = 
      child_generator.GenerateChildren(*el);
    for(int i=0; i<childels.size(); i++){
      //retpair ret = subelements.insert(std::move(childels[i]));
      auto it = subelements.find(childels[i]);
      if(it != subelements.end()){
	signed char orient = ElementOrientation(childels[i],*it,index_factory);
	assert((*it)->getBCTag() != 0);
	if(orient != 1){
	  unique_element_ptr newel = 
	    getOrientedEl(*it,index_factory,orient);
	  subelements.erase(*it);
	  subelements.insert(std::move(newel));
	}
      }
    }
  }

}

void generateSubelements(const element_set& elements, element_set& subelements){
  ChildGenerator child_generator;
  NodeIndexFactory index_factory;
  int cnt=0;
  int pos_cnt = 0, neg_cnt = 0;
  int pos_count[] = {0,0,0};
  int neg_count[] = {0,0,0};
  
  //std::cout << "Initial size of subelements: " << subelements.size() << std::endl;
  typedef std::pair<element_set::iterator,bool> retpair;
  for(auto el=elements.begin(); el != elements.end(); ++el){
    std::vector<unique_element_ptr>& childels = 
      child_generator.GenerateChildren(*el);
    //std::cout << "childels size: " << childels.size() << std::endl;
    for(int i=0; i<childels.size(); i++){
      //(*el)->setChild(i,childels[i].get(),1);
      retpair ret = subelements.insert(std::move(childels[i]));
      
      if(ret.second == false){
	signed char orient = 
	ElementOrientation(childels[i],*(ret.first),index_factory);
	if(0){
	  //if((*ret.first)->getBCTag() != 0){
	  unique_element_ptr newel = 
	    getOrientedEl(*ret.first,index_factory,orient);
	  subelements.erase(ret.first);
	  auto it = subelements.insert(std::move(newel));
	  (*el)->setChild(i,it.first->get(),1);
	}
	else{
	  (*el)->setChild(i,ret.first->get(),orient);
	}
      }
      else{
	(*el)->setChild(i,ret.first->get(),1);
      }
      int orient = (*el)->getChildOrientation(i);
      if(orient > 0){
	pos_count[orient-1]++;
      }
      else if (orient < 0){
	neg_count[-orient-1]++;
      }
      else{
	std::cout << "orient is 0!" << std::endl;
      }
    }

  }
  /*
  std::cout << "Number of orientations: " << std::endl;
  std::cout << "pos: ";
  for(int i=0; i<3; i++){
    std::cout << " " << pos_count[i];
  }
  std::cout << " neg: ";
  for(int i=0; i<3; i++){
    std::cout << " " << neg_count[i];
  }
  std::cout << std::endl;
  */

  //std:cout << "number of non-inserted els: " << cnt << std::endl;
}

/*
void generateSubelements2(element_set& elements, int dim){
  ChildGenerator child_generator;
  typedef std::pair<element_set::iterator,bool> retpair;
  int els_inserted=0;
  int els_attempted =0;
  for(auto el=elements.begin(); el != elements.end(); ++el){
    if((*el)->getDim() == dim){
      std::vector<unique_element_ptr>& childels = 
	child_generator.GenerateChildren(*el);
 
      for(int i=0; i<childels.size(); i++){
	(*el)->setChild(i,childels[i].get());
	assert(childels[i]);
	
	retpair ret = elements.insert(std::move(childels[i]));
	
	if(ret.second == false){
	  (*el)->setChild(i,(*ret.first).get());
	}
	else els_inserted++;
	
	els_attempted++;
	const MEl* child = (*el)->getChild(i);
      }
    }
  }
  cout << "els_inserted:  " << els_inserted << endl;
  cout << "els_attempted: " << els_attempted << endl;
}
*/
void ComputeMeshConnectivity(MeshContainer& mesh){
  int node_size = 0;
  const node_map& nodes = mesh.getNodes();

  element_set& elements = mesh.getElementsNC();
  element_set& subelements = mesh.getSubElementsNC();
  element_set& subsubelements = mesh.getSubSubElementsNC();

  //subelements.max_load_factor(2.0);

  //const element_map& elements = mesh.getElements();
  //element_map& subelements = mesh.getSubElementsNC();
  //element_map& subsubelements = mesh.getSubSubElementsNC();

  /*
  cout << "number of elements: " << elements.size() << endl;
  cout << "number of subelements: " << subelements.size() << endl;
  cout << "number of subsubelements: " << subsubelements.size() << endl;
  cout << "sizeof(MEl): " << sizeof(MEl) << endl;
  cout << "sizeof(El1D): " << sizeof(El1D) << endl;
  cout << "sizeof(El2D): " << sizeof(El2D) << endl;
  cout << "sizeof(MLine): " << sizeof(MLine) << endl;
  cout << "sizeof(MTri): " << sizeof(MTri) << endl;   
  cout << "sizeof(Mquad): " << sizeof(MQuad) << endl;
  cout << "sizeof(MTet): " << sizeof(MTet) << endl;
  */

  std::cout << "Number of elements: " << elements.size() << std::endl;
  std::cout << "Number of subelements: " << subelements.size() << std::endl;

  {
    int symm_cnt = 0;
    for(auto it = subelements.begin(); it != subelements.end(); ++it){
      auto el = it->get();
      //if(el->getOrder() != 2) std::cout << "Order not 2!" << std::endl;
      if(el->getBCTag() == 5) symm_cnt++;
    }
    std::cout << "Number of symmetry faces: " << symm_cnt << std::endl;
  }


  arma::wall_clock timer;
  timer.tic();

  for(int d = 1; d< mesh.Dimension(); d++){
    for(int se = d; se >= 1; se--){
      generateSubelements(mesh.getElementsOfDim(se+1),mesh.getElementsOfDim(se));
    }
  }

  /*
  generateSubelements(elements,subelements);
  
  if(mesh.MeshDimension() == 3){
    //generateSubelements2(elements,2);
    generateSubelements(subelements,subsubelements);
    //generateSubelements2(subelements);
    //generateSubelements(subelements,subelements);
  }
  */
  //cout << "subelement time: " << timer.toc() << endl;


  int  mesh_bytes = (sizeof(MTet)+8)*elements.size() + 
    (sizeof(MTri)+8)*subelements.size() +
    (sizeof(MLine)+8)*subsubelements.size();

  int pos = 0;
  int neg = 0;
  for(auto el=elements.begin(); el != elements.end(); ++el){
    for(int i=0; i<(*el)->NumChildren(); i++){
      //std::cout << int((*el)->getChildOrientation(i)) << std::endl;
      if(int((*el)->getChildOrientation(i)) > 0) pos++;
      else neg++;
    }
  }
  
  cout << "number of positive orientations: " << pos << endl;
  cout << "number of negative orientations: " << neg << endl;
  std::cout << "number of faces: " << subelements.size() << std::endl;

  
  //cout << "Mesh in GB: " << mesh_bytes/1e9 << endl;
  //cout << "NOdes size: " << (sizeof(MNode)+12)*nodes.size()/1e9 << endl;
  cout << "number of subelements: " << subelements.size() <<  endl;
  cout << "number of subsubelements: " << subsubelements.size() << endl;
  

}

//int HModel::getMeshConnectivity(){

  /*
  int const *array = new int[3];

  cout << "sizeof(MEl): " << sizeof(MEl) << endl;
  cout << "sizeof(El2D): " << sizeof(El2D) << endl;
  //cout << "sizeof(MTRI: " << sizeof(MTri) << endl;
  cout << "sizeof(int[3]): " << sizeof(int[3]) << endl;
  cout << "sizeof(array): " << sizeof(array) << endl;
  cout << "sozeof(unique_ptr): " << sizeof(std::unique_ptr<MEl>) << endl;

  arma::wall_clock timer;
  timer.tic();

  //element_set::iterator it;
  typedef std::pair<element_set::iterator,bool> insertPair;
  for(auto el=mesh_elements.begin(); el != mesh_elements.end(); ++el){
    const std::vector<element_ptr> children = (*el)->generateChildren();
    for(int i=0; i<children.size(); i++){
      insertPair flag = subelements.insert(children[i]);
      if(flag.second == false){
	signed char orient = elOrientation(children[i],(*flag.first));
	MEl* parent = children[i]->getParent();
	unsigned char parent_face = children[i]->getParentFace();
	parent->setChild(parent_face,orient,(*flag.first));
      }
    }
  }

  if(mesh_dim == 3){
    cout << "In mesh dim == 3" << endl;
    for(auto el=subelements.begin(); el != subelements.end(); ++el){
      const std::vector<element_ptr> children = (*el)->generateChildren();
      for(int i=0; i<children.size(); i++){
	insertPair flag = subelements.insert(children[i]);
	if(flag.second == false){
	  signed char orient = elOrientation(children[i],(*flag.first));
	  MEl* parent = (*flag.first)->getParent();
	  unsigned char parent_face = (*flag.first)->getParentFace();
	  parent->setChild(parent_face,orient,children[i]);
	}
      }
    }
  }
  */

  //cout << "mesh connectivity time: " << timer.toc() << endl;


  //return 1;
//}
