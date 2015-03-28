#pragma once
#include "MeshTypedefs.h"

#include <memory>
#include <unordered_map>
#include <vector>
#include <armadillo>

class MeshContainer;
class GeometryContainer;
class IndexFactory;
class NodeInterpolator;

class HighOrderMeshGenerator{
 public:
  HighOrderMeshGenerator(MeshContainer& mesh, GeometryContainer&geometry);
  std::shared_ptr<MeshContainer> generateHighOrderMeshRecursive(int order);
  
 private:
  int meshDimension(int dim, int order);
  
  MeshContainer& mesh;
  GeometryContainer& geometry;
  std::shared_ptr<NodeIndexFactory> index_factory;

  std::shared_ptr<NodeInterpolator> node_interpolator;

  std::shared_ptr<MeshContainer> meshcurr, meshold;

  std::unordered_map<const MEl*,MEl*> elmap[4];
  
  std::vector<double> 
    getElementInteriorNodesCoordinates(unique_element_ptr& el);

  arma::mat getHighOrderNodeCoordinates(unique_element_ptr& el, int order);

  std::vector<int> insertNewNodes(const unique_element_ptr& el,
				  const arma::mat& newelnodes,
				  const arma::mat& parametric_coords,
				  int order);

  int insertNewElement(const unique_element_ptr& el,
		       std::vector<int>& nodeindices,
		       int order);

  int setNewChildrenAndOrientation();

  int setMeshIsCurved();
  
};
