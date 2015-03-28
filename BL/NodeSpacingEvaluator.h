#pragma once
#include <vector>

class BLMesh;

class NodeSpacingEvaluator{
 public:
 NodeSpacingEvaluator(BLMesh& mesh): mesh(mesh){}
  virtual std::vector<double> computeSpacing() = 0;

 protected:
  BLMesh& mesh;
};
