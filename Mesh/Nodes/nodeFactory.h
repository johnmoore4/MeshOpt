#pragma once

//#include "Mel.h"
#include "GlobalDefines.h"

#include <memory>

class nodeFactory;
class MNode;


class nodeFactory{
 private:
  static nodeFactory* _instance;
  static gind NODE_COUNT;
  
 protected:
 
  nodeFactory(){}
  nodeFactory(const nodeFactory &){}
  void operator=(const nodeFactory&){}

 public:
  static nodeFactory* Instance();

  std::unique_ptr<MNode> 
    CreateNode(int ndtype, double xyz[3], 
	       short int geo_entity=-1, double u=0, double v=0);
  gind GetNodeCount();

  

};
