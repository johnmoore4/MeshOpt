#pragma once

//#include "Mel.h"
#include "GlobalDefines.h"

#include <memory>
#include <vector>

class NodeIndexer;

class elementFactory;
class MNode;
class MEl;

class elementFactory{
 private:
  static elementFactory* _instance;
  static gind ELEMENT_COUNT;

 protected:
 
  elementFactory(){}
  elementFactory(const elementFactory &){}
  void operator=(const elementFactory&){}

 public:
  static elementFactory* Instance();

  std::unique_ptr<MEl> 
    CreateElement(int elType, 
		  std::vector<gind>& nodes, 
		  const NodeIndexer* indexer,
		  unsigned char order, 
		  unsigned char bc_tag, 
		  short int geo_entity=-1, 
		  bool geo_type=0);

  gind GetElementCount(){ return ELEMENT_COUNT; }

};
