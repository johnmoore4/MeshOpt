#pragma once 
#include <string>

class MeshContainer;

class MatlabMeshWriter{
 public:
 MatlabMeshWriter(const MeshContainer& mesh): mesh(mesh){}
  void Write(std::string filename);

 private:
  const MeshContainer& mesh;
};
