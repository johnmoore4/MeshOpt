#include "ParameterManip.h"
#include <sstream>
#include <iostream>
#include <stdlib.h>

//#include "stdlib.h"
//#include <boost/algorithm/string.hpp>

void ParameterManip::OverrideString(std::string &str, const std::string& val, 
		   const parameter_map_type& map){

  auto it = map.find(val);
  if(it != map.end()){
    str = it->second;
  }

}

void ParameterManip::OverrideInteger(int &in, const std::string& val,
		    const parameter_map_type& map){
  auto it = map.find(val);
  if(it != map.end()){
    in = atoi(it->second.c_str());
  }

}

void ParameterManip::OverrideDouble(double &db, const std::string& val,
		   const parameter_map_type& map){

  auto it = map.find(val);
  if(it != map.end()){
    db = atof(it->second.c_str());
  }

}

void ParameterManip::OverrideBool(bool & bl, const std::string& val, 
				 const parameter_map_type& map){
  auto it = map.find(val);
  if(it != map.end()){
    if(it->second == "true" || it->second == "True") bl = true;
    else if(it->second == "false" || it->second == "False") bl = false;
    else{
      throw std::runtime_error("val must be \"true\" or \"false\"");
    }
  }
}

void ParameterManip::OverrideVector(std::vector<double>& vec,
				    const std::string& val, 
				    const parameter_map_type& map){


  vec.resize(0);

  auto it = map.find(val);
  if(it != map.end()){
    std::string vs = it->second;
    if(vs.front() == '{' && vs.back() == '}'){
      std::string temp = vs.substr(1,vs.size()-1);
      std::stringstream ss(temp);
      char vl;
      std::string substr;
      while(ss >> vl){
	//vec.push_back(vl);
	substr.push_back(vl);
	if(ss.peek() == ',' || ss.peek() == '}'){
	  ss.ignore();
	  vec.push_back(atof(substr.c_str()));
	  substr.resize(0);
	}
      }
      
    }
    else{
      throw std::runtime_error
	("Vector input must be enclosed with brackets {}!");
    }
  }


}
