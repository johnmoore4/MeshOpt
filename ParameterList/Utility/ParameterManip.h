#include <string>
#include <vector>
#include <map>

class ParameterManip{

typedef std::map<std::string,std::string> parameter_map_type;
 public:
 static void OverrideString(std::string &str, const std::string& val, 
			    const parameter_map_type& map);

  static void OverrideInteger(int &in, const std::string& val,
		      const parameter_map_type& map);

  static void OverrideDouble(double &db, const std::string& val,
		     const parameter_map_type& map);
  

  static void OverrideBool(bool & bl, const std::string& val, 
			   const parameter_map_type& map);

 
  static  void OverrideVector(std::vector<double>& vec, const std::string& val, 
			      const parameter_map_type& map);
};
