#include "MyParameterList.h"
//#include <boost/mpi.hpp>

void MyParameterList::Initialize(const parameter_map_type& plist){
 
  FindInvalidParameters(plist,GetValidParameters());
  Populate(plist);
  //Broadcast();
  
}

void MyParameterList::
FindInvalidParameters(const parameter_map_type& plist,
		      const std::set<std::string>& valid) const{
  for(auto it = plist.begin(); it != plist.end(); ++it){
    auto fn = valid.find(it->first);
    if(fn == valid.end()){
      throw std::runtime_error("Invalid parameter: " + *fn + "!");
    }

  }

}
