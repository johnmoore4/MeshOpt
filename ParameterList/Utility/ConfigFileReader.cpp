#include "ConfigFileReader.h"

#include <iostream>

#include <fstream>
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>

bool ConfigFileReader::IgnoreLine(std::string line){
  
  if(line.length() == 0 || line[0] == '/' || line[0] == '#') return true;
  if(line.find_first_not_of(' ') == std::string::npos) return true;
  if(line.find_first_not_of('\t') == std::string::npos) return true;
  else return false;
}

int ConfigFileReader::ReadConfigFile(std::string filename){
  using namespace std;
  ifstream in;
  in.open(filename);

  string line;
  while(std::getline(in,line)){

    if(!IgnoreLine(line)){

      size_t pos_bracket = line.find("{");
      if(pos_bracket != string::npos){
	std::string category = line.substr(0,pos_bracket);
	parameter_map_type temp_map;
	//cout << "category: " << category << endl;

	std::string line2;
	while(getline(in,line2)){
	  if(!IgnoreLine(line2)){
	    if(line2[0] == '}') break;
	  
	    size_t pos_equal = line2.find("=");
	    if(pos_equal == string::npos){
	      throw std::runtime_error("Parameter list assignment must be " 
				       "seperated by an equal sign!");
	    }
	    size_t pos_semicolon = line2.find(";");
	    if(pos_semicolon == string::npos){
	      throw std::runtime_error("Parameter list assignments must be " 
				       "terminated by a semicolon!");
	    }

	    std::string s1 = line2.substr(0,pos_equal);
	    std::string s2 = 
	      line2.substr(pos_equal+1,pos_semicolon-pos_equal-1);
	    boost::algorithm::trim(s1);
	    boost::algorithm::trim(s2);
	    temp_map[s1] = s2;
	  }
	}
	parameter_list[category] = temp_map;
      }
    }

  }
  in.close();

}

all_parameter_type& ConfigFileReader::GetParameterList(){
  return parameter_list;

}
