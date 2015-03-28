#pragma once
#include "MeshTypedefs.h"
#include <memory>

class ActiveMEl;
class NodeIndexFactory;

class MEl_to_ActiveEl{
 public:
  MEl_to_ActiveEl(const node_map& nodes,
		  std::shared_ptr<NodeIndexFactory>& index_factory);
  ActiveMEl get(const unique_element_ptr& el);
  
 private:  
  const node_map& nodemap;
  std::shared_ptr<NodeIndexFactory>& index_factory;

};
