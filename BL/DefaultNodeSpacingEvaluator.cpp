#include "DefaultNodeSpacingEvaluator.h"
#include "BLMesh.h"
#include <cmath>
#include <iostream>

std::vector<double> DefaultNodeSpacingEvaluator::computeSpacing(){
  int Nbl_nodes = mesh.blnodes.size();
  int Nrows = NLayers+1;
  
  std::vector<double> spacings(Nbl_nodes*Nrows);
  using std::pow;

  double positions[Nrows];
  
  
  for(int row = 0; row < Nrows; row++){
    positions[row] = double(row)/NLayers*pow(growthRatio,row)/pow(growthRatio,NLayers);
    std::cout << "pos: " << positions[row] << std::endl;
  }

  for(int nd = 0; nd < Nbl_nodes; nd++){
    for(int row = 0; row < Nrows; row++){
      spacings[nd*Nrows+row] = positions[row];
    }
  }
  
  return spacings;
}
