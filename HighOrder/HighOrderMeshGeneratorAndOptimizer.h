#pragma once
#include <memory>
class MeshContainer;
class GeometryContainer;

class HighOrderMeshGeneratorAndOptimizer{
 public:
  HighOrderMeshGeneratorAndOptimizer(MeshContainer& mesh, 
				     GeometryContainer&geometry);

  std::shared_ptr<MeshContainer> GenerateAndOptimize(int order);

 private:
  MeshContainer& mesh;
  GeometryContainer& geometry;
};
