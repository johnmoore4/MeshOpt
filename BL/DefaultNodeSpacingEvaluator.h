#pragma once
#include "NodeSpacingEvaluator.h"

class DefaultNodeSpacingEvaluator: public NodeSpacingEvaluator{
 public:
  DefaultNodeSpacingEvaluator(BLMesh& mesh, double minSpacing, 
			      double growthRatio, int NLayers):
  NodeSpacingEvaluator(mesh), minSpacing(minSpacing), 
    growthRatio(growthRatio), NLayers(NLayers){}

  std::vector<double> computeSpacing();

 private:
  const double minSpacing;
  const double growthRatio;
  const int NLayers;
};
