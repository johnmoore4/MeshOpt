#include "nodeFactory.h"
#include "Mnode.h"
#include "assert.h"

nodeFactory* nodeFactory::_instance=0;
gind nodeFactory::NODE_COUNT=0;

gind nodeFactory::GetNodeCount(){ return NODE_COUNT; }


std::unique_ptr<MNode> nodeFactory::CreateNode
(int ndtype, double xyz[3], 
 short int geo_entity, double u, double v){

  
  switch(ndtype){
  case 0:
    {
      return std::unique_ptr<MGeoNode>(new MGeoNode(NODE_COUNT++,
						    xyz,geo_entity));
    }
  case 1:
    {
      return std::unique_ptr<MEdgeNode>(new MEdgeNode(NODE_COUNT++,
						     xyz,geo_entity,u));
    }
  case 2:
    {
      return std::unique_ptr<MFaceNode>(new MFaceNode(NODE_COUNT++,
						     xyz,geo_entity,u,v));
    }
  case 3:
    {
      return std::unique_ptr<MNode>(new MNode(NODE_COUNT++,
						    xyz));
    }
  default:
    {
      assert("MNode type should be less than 4!");
    }
  }
 

}
nodeFactory* nodeFactory::Instance(){
 
  if(!_instance){
    _instance = new nodeFactory;
  }
  return _instance;
}
