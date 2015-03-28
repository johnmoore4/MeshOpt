#pragma once
#include "GlobalDefines.h"
#include <unordered_set>
#include <map>
#include <vector>
#include <armadillo>
#include <memory>
#include "MeshTypedefs.h"

class MeshContainer;
class OptElManager;
class MEl;
class MNode;

//typedef std::map<gind,std::unique_ptr<MNode> > node_map;


class MeshOptimizer{

 public:
  MeshOptimizer(MeshContainer& mesh, OptElManager& optel_manager,
		const std::map<const MEl*,arma::mat>& ideal_map,
		int el_dim_to_opt);

  virtual int Optimize(double threshold, std::vector<bool> to_opt,
		       std::set<unsigned short int> geo_faces_to_opt = 
		       std::set<unsigned short int>(),
		       std::unordered_set<int> clamped_faces =
		       std::unordered_set<int>()) = 0;
  int UpdateMeshMerits(const std::unordered_set<MEl*>& watch);
  int UpdateMeshMerits();
  const double minMeshQuality() const;


 protected:
  MeshContainer& mesh;
  OptElManager& optel_manager;
  std::map<const MEl*,arma::mat> idealElements;
  std::map<const MEl*,double> all_merits;
  //int only_opt_node_with_bctag = -1;
  std::set<unsigned short int> geo_faces_to_opt;
  std::unordered_set<int> clamped_nodes;

  const int el_dim_to_opt;
  //ShapeFunctionMatricesFactory& sf_factory;
  //OptElManager& optel_manager;

  std::vector<bool> dims_to_opt = {false,false,true,false};

  void FindActiveElements(std::unordered_set<const MEl*>& watch,
			  std::map<gind,gind>& activenodes,
			  const double threshold,
			  const int Nnearest,
			  const int type_to_find=-1,
			  const int type_to_exclude = -1);

  std::map<gind, std::vector<const MEl*> >
    CreateNode2ElementList(const std::unordered_set<const MEl*>& watch,
			   const std::map<gind,gind>& activenodes);

 private:
 void InsertActiveNodes(const MEl* el, 
			std::map<gind,gind>& activenodes,
			const node_map& nodes);

};

