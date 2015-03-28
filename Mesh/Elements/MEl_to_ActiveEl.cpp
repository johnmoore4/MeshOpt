#include "MEl_to_ActiveEl.h"
#include "ActiveMEl.h"
#include "NodeIndexer.h"

MEl_to_ActiveEl::
MEl_to_ActiveEl(const node_map& nodes,
		std::shared_ptr<NodeIndexFactory>& index_factory): 
  nodemap(nodes), index_factory(index_factory){
  
}

ActiveMEl MEl_to_ActiveEl::get(const unique_element_ptr& el){
  //const NodeIndexer* ni = 
  //  index_factory->getNodeIndexer(el->getElementType(),el->getOrder());

  return ActiveMEl(el.get(),*index_factory,nodemap,NULL);
}
