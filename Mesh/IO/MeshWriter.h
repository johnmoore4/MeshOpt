#pragma once
#include "MeshContainer.h"

class MeshWriter{
 protected:
  MeshContainer& mesh;
 public:
 MeshWriter(MeshContainer& mesh_t): mesh(mesh_t){}

  virtual void Save(std::string filename,std::string format)=0;

};

class GMSHWriter: public MeshWriter{
  

 public:
 GMSHWriter(MeshContainer& mesh_t): MeshWriter(mesh_t) {}
  void Save(std::string filename="default",std::string format="");
  

};
