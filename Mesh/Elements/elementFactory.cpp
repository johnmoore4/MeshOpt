#include "elementFactory.h"
#include "Mel.h"
#include "El1D.h"
#include "El2D.h"
#include "El3D.h"

elementFactory* elementFactory::_instance=0;

gind elementFactory::ELEMENT_COUNT=0;

std::unique_ptr<MEl> elementFactory::CreateElement
(int elType, std::vector<gind>& nodes, const NodeIndexer* indexer, 
 unsigned char order, 
 unsigned char bc_tag, short int geo_entity, 
 bool geo_type){

  ELEMENT_COUNT++;
  switch(elType){
  case 1:
    {
      return std::unique_ptr<MLine>(new MLine(nodes,indexer,order,bc_tag,
					      geo_entity,geo_type));
    }
  case 2:
    {
      return std::unique_ptr<MTri>(new MTri(nodes,indexer,order,bc_tag,
					    geo_entity, geo_type));
    }
  case 3:
    {
      return std::unique_ptr<MQuad>(new MQuad(nodes,indexer,order,bc_tag,
      					    geo_entity, geo_type));
    }
  case 4:
    {
      return std::unique_ptr<MTet>(new MTet(nodes,indexer,order,bc_tag,
      					    geo_entity, geo_type));
    }
  case 5:
    {
      return std::unique_ptr<MHex>(new MHex(nodes,indexer,order,bc_tag,
      					    geo_entity, geo_type));
    }
  case 6:
    {
      return std::unique_ptr<MPrism>(new MPrism(nodes,indexer,order,bc_tag,
      					    geo_entity, geo_type));
    }
  case 7:
    {
      return std::unique_ptr<MPyramid>
	(new MPyramid(nodes,indexer,order, bc_tag,geo_entity,geo_type));
    }
  default:
    { 
      throw std::runtime_error("Unsupported element type of " + std::to_string(elType));
    }
  }
  //std::unique_ptr<MEl> elementFactory::CreateElement(int ID) const {
  //return std::unique_ptr<El1D>(El1D(elno,);

}


elementFactory* elementFactory::Instance(){
 
  if(!_instance){
    _instance = new elementFactory;
  }
  return _instance;
}

