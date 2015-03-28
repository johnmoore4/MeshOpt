#pragma once
#include "LocalMesh.h"

#include "CEl.h"
#include <array>
#include <boost/serialization/access.hpp>

class DGFaceInfo{
 public:
  int l, r;
  unsigned char orl=0, orr=0;
  unsigned char pfl, pfr;

 private:
  friend class boost::serialization::access;

  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & l;
      ar & r;
      ar & orl;
      ar & orr;
      ar & pfl;
      ar & pfr;
    }
};

class DGMesh: public LocalMesh{
 public:

  int Write(std::string format){}
  int Read(std::string format){}


  std::vector<DGFaceInfo> faces;
  std::vector<unsigned int> numElDofs;
//std::vector<CEl> elements;
//std::vector<std::array<double,3> > nodes;
  std::vector<int> globalElementPtrs;
//std::vector<unsigned char> bc_tags;
  int NElLocal, NDGDofsGlobal;
//int Dimension;

 private:
  friend class boost::serialization::access;

  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & node_numbers;
      ar & faces;
      ar & numElDofs;
      ar & elements;
      ar & nodes;
      ar & globalElementPtrs;
      ar & bc_tags;
      ar & NElLocal;
      ar & NDGDofsGlobal;
      ar & Dimension;
    }
};
