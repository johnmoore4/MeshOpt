#pragma once
#include "MeshOptimizer.h"
#include "OptElManager.h"



class SubmeshOptimizer: public MeshOptimizer{
 private:
  //const normal_map_type& normal_map;
  //const gindmap& newnode_map;
  //const std::map<MEl*, bool>& ExtrudedIdeal;
 protected:
 
 public:
 SubmeshOptimizer(MeshContainer& mesh,
		  OptElManager& optel_manager,
		  const std::map<const MEl*, arma::mat>& ideal_map =
		  std::map<const MEl*, arma::mat>(),
		  int el_dim_to_opt=2): 
  MeshOptimizer(mesh, optel_manager,ideal_map,el_dim_to_opt){}

  int Optimize(double threshold = 1.0, std::vector<bool> to_opt = {1,1,1},
	       std::set<unsigned short int> geo_faces_to_opt=
	       std::set<unsigned short int>(),
	       std::unordered_set<int> clamped_nodes_t = 
	       std::unordered_set<int>());
  //int SnapToBLExtent();
  double OptimizeOneByOne(double threshold=1.0, int Nnearest = 0, int type = -1,
		       double factor=1.0);

  int OptimizeGlobal(double threshold=1.0, int Nnearest = 0, int type = -1,
		      double factor=1.0);

  int RandomizeNodes(double factor);
};
