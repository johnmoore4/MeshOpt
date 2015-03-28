#pragma once

#include <memory>

class MeshContainer;
class BLParameterList;


class BoundaryLayerSubdivider{
 public:
BoundaryLayerSubdivider(MeshContainer& mesh, 
			std::shared_ptr<BLParameterList> bl_parameters):
  mesh(mesh), bl_parameters(bl_parameters){}

  int GenerateElements();

 private:
  MeshContainer& mesh;
  std::shared_ptr<BLParameterList> bl_parameters;
};
