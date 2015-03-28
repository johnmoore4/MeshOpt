#pragma once

#include <vector>
#include <map>
#include <memory>
#include <set>

class MeshContainer;
class MEl;

class BLMesh{
 public:
  BLMesh(const MeshContainer& mesh);

  std::map<int,int> blnodemap;
  std::vector<int> blnodes;
  std::vector<int> eptrs;
  std::vector<int> eind;
  std::vector<std::shared_ptr<MEl> > blels;
  std::set<const MEl*> symmetry_face_elements;
  std::vector<int> symmetry_faces;

  //std::vector<const MEl*> bl_faces;

 private:
   

};
