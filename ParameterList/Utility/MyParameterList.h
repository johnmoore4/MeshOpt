#pragma once
#include <set>
#include <map>
#include <string>


typedef std::map<std::string,std::string> parameter_map_type;

class MyParameterList{
 public:

  void Initialize(const parameter_map_type& pmap);
  virtual void Broadcast() = 0;
 protected:

  
  virtual void Populate(const parameter_map_type& pmap) = 0;

  virtual std::set<std::string> GetValidParameters() const = 0;

 private:
  void FindInvalidParameters(const parameter_map_type& plist,
			     const std::set<std::string>& valid) const;
};
