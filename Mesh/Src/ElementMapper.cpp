#include "ElementMapper.h"
#include "MeshTypedefs.h"
  
unique_element_ptr& ElementMapper::getMappedElement(const unique_element_ptr& from){

}

int ElementMapper::insertElementPair(const unique_element_ptr& from,
				     const unique_element_ptr& to){

}

int ElementMapper::Erase(){
  for(int i = 0; i < 3; i++) elmap[i].erase();
}
