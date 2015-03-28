#pragma once
#include "MeshOptimizer.h"
#include "OptElManager.h"
#include "BoundaryLayer.h"


class GlobalMeshOptimizer: public MeshOptimizer{

 public:
 GlobalMeshOptimizer(MeshContainer& mesh_t,
		     OptElManager& optel_manager_t,
		     const std::map<const MEl*,arma::mat>& ideal=
		     std::map<const MEl*,arma::mat>(),
		     int el_dim_to_opt=2): 
  MeshOptimizer(mesh_t, optel_manager_t,ideal,el_dim_to_opt){}

  int Optimize(double threshold = 1.0, std::vector<bool> to_opt = {1,1,1});




};
