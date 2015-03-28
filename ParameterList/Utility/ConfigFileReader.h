#pragma once

#include "MyParameterList.h"

//typedef std::map<std::string,std::string> parameter_map_type;

typedef std::map<std::string, parameter_map_type > all_parameter_type;



class ConfigFileReader{
 private:
  all_parameter_type parameter_list;

  bool IgnoreLine(std::string line);
 public:
  int ReadConfigFile(std::string filename);
  all_parameter_type& GetParameterList();

};

