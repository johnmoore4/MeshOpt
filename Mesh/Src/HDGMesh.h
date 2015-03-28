#pragma once
#include "LocalMesh.h"


class HDGMesh: public LocalMesh{
 public:

  int Nfdof_local, Nfdof_global, StartDof;
  int Nfaces_global, Nfaces_local;
  int StartFace;

  std::vector<bool> is_local;  
  std::vector<int> start_dofs;
  std::vector<unsigned char> num_dofs;

  std::vector<int> start_global_dofs;
  std::vector<int> end_global_dofs;
  std::vector<int> face_nodes;
  std::vector<int> edg_node_map;
  std::vector<int> global_faces;
  std::vector<int> FirstGlobalPointInBlocks;
  std::vector<int> BlockSizes;
  std::vector<int> node_perm;
  std::vector<std::vector<int> > faceNodes;
  std::vector<int> reordered_nodes;

  int Write(std::string format);
  int Read(std::string format);

  
 private:

  friend class boost::serialization::access;

  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & elements;
      ar & nodes;
      ar & node_numbers;
      ar & bc_tags;
      ar & Dimension;

      ar & start_dofs;
      ar & num_dofs;
      ar & face_nodes;
      ar & edg_node_map;
      ar & is_local;
      ar & Nfdof_local;
      ar & Nfdof_global;
      ar & StartDof;
      ar & global_faces;
      ar & Nfaces_global;
      ar & Nfaces_local;
      ar & FirstGlobalPointInBlocks;
      ar & BlockSizes;
      ar & StartFace;
      ar & node_perm;
      ar & faceNodes;
      ar & reordered_nodes;
    }

};
