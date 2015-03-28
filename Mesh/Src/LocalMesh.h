#pragma once
#include "GlobalDefines.h"
#include "CEl.h"

#include <vector>
#include <unordered_map>
#include <map>
#include <memory>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>



class LocalMesh{
 private:

 protected:
  friend class boost::serialization::access;

  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & elements;
      ar & nodes;
      ar & node_numbers;
      ar & bc_tags;
      ar & Dimension;


    }

 public:

virtual int Write(std::string filename) = 0;
virtual int Read(std::string filename) = 0;

  std::vector<CEl> elements;
  std::vector<std::array<double,3> > nodes;
  std::vector<gind> node_numbers;
  std::vector<unsigned char> bc_tags;

  int Dimension;
  

};
