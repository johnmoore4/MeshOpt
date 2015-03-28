#pragma once
#include "GlobalDefines.h"

#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
//#include <boost/serialization/vector.hpp>


class CEl{
 private:
  friend class boost::serialization::access;

  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & nodes;
      ar & children;
      ar & bc_tag;
      ar & el_type;
      ar & porder;
      ar & elno;
      ar & dim;
      ar & is_curved;
    }

 
  
 protected:

 public:
  CEl(){}
  CEl(std::vector<gind> nodes, 
      std::vector<std::pair<int, signed char> >children,
      unsigned char el_type, 
      unsigned char porder, 
      short int bc_tag,
      gind elno,
      unsigned char dim,
      bool is_curved):
  nodes(nodes), children(children), el_type(el_type), porder(porder), 
    bc_tag(bc_tag), elno(elno), dim(dim), is_curved(is_curved){}

  std::vector<int> nodes;
  std::vector<std::pair<int, signed char> > children;
  gind elno;
  short int bc_tag;
  unsigned char el_type, porder, dim;
  bool is_curved;

};
