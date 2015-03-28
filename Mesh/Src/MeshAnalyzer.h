#pragma once

#include <memory>
#include <vector>

class OptElManager;
class MeshContainer;

class MeshAnalyzer{
 public:
  MeshAnalyzer(MeshContainer& mesh);
  int Analyze();
  int SetMeshCurvedElements();

  //double minQuality(){ return min_quality; }
  //double maxQuality(){ return max_quality; }

  //const std::vector<double>& getMeshQualities(){ return qualities; }

 private:
  int ComputeMeshQualities(int dim);
  
  //int ComputeQuality();
  
  //double min_quality = 0;
  //double max_quality = 1;
  
  MeshContainer& mesh;
  
  std::vector<double> qualities[3];

};
