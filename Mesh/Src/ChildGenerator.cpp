#include "ChildGenerator.h"
#include "GlobalDefines.h"
#include "elementFactory.h"
#include "Mel.h"



std::vector<element_ptr>& ChildGenerator::
GenerateChildren(const element_ptr& el){
  //const gind* allnodes = el->getCornerNodes();

  const int bc_tag = el->getBCTag();
  const int geo_entity = el->getGeoEntity();
  const bool geo_type = el->getGeoType();
  const int type = el->getElementType();
  const unsigned char order = el->getOrder();
  //std::cout << "order: " << int(order) << std::endl;

  const NodeIndexer* indexer = 
    index_factory.getNodeIndexer(type,order);

  const gind* ndind = el->getNodes();

  for(int i=0; i<indexer->Ndof(); i++){
    //std::cout << ndind[i] << " ";
  }
  //std::cout << std::endl;

  children.resize(indexer->NumChildren());
  //cout << "geo entity: " << geo_entity << " ";

  for(int i=0; i<indexer->NumChildren(); i++){
    //gind nodes[indexer->NumChildCornerNodes(i)];
    //gind nodes[4];
    //const int* cn = indexer->ChildCornerNodes(i); 
    const std::vector<indtype>& cnodes = indexer->getChildNodes(i);
    std::vector<int> cnvec(cnodes.size());
    
    for(int j=0; j<cnodes.size(); j++){
      cnvec[j] = ndind[cnodes[j]];
    }
 

    //for(int j=0; j<indexer->NumChildCornerNodes(i); j++){
    //  nodes[j] = allnodes[cn[j]];
    //}
    //std::cout << "child type: " << indexer->ChildType(i) << std::endl;

    const NodeIndexer* child_indexer = 
      index_factory.getNodeIndexer(indexer->ChildType(i),order);

    children[i] = elementFactory::Instance()->
      CreateElement(indexer->ChildType(i),cnvec,child_indexer,
		    order,bc_tag,geo_entity,geo_type);
    //const gind* temp = children[i]->getNodes();
   
  }
  //cout << endl;
  return children;
 

}
