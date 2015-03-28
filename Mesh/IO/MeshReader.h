#pragma once
#include <string>
#include <map>

class MeshContainer;
class GeometryContainer;

class MeshReader{
 private:

 protected:

 public:
  virtual void ReadMesh(MeshContainer& mesh, std::string filename)=0;
  virtual void SetGeoTags(GeometryContainer& geometry) = 0;
};


class GMSHReader: public MeshReader{
 private:

 protected:
  std::map<int,int> edge2bc_map, face2bc_map;
 public:
  
  void ReadMesh(MeshContainer& mesh, std::string filename);
  void SetGeoTags(GeometryContainer& geometry);

};
