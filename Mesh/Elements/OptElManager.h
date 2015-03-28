#pragma once
#include "OptEl.h"
#include <armadillo>
#include <memory>

class ShapeFunctionMatricesFactory;
class NodeIndexFactory;
class GeometryContainer;

class ElIdealPair{

 public:
  const MEl* el;
  const arma::mat& ideal;
 ElIdealPair(const MEl* el_t, const arma::mat& ideal_t): 
  el(el_t), ideal(ideal_t){}
};

class OptElManager{
 private:
  node_map& nodes;
  GeometryContainer& geometry;
  ShapeFunctionMatricesFactory& sf_factory;
  NodeIndexFactory& node_indexer;


 public:
  OptElManager(node_map& nodes_t,
	       GeometryContainer& geometry_t,
	       ShapeFunctionMatricesFactory& sf_factory_t,
	       NodeIndexFactory& node_indexer_t): 
  nodes(nodes_t), 
    geometry(geometry_t), 
    sf_factory(sf_factory_t), 
    node_indexer(node_indexer_t){}

  std::unique_ptr<OptEl> CreateOptEl(const MEl* el, const arma::mat& ideal);

  arma::mat CreateIdealMatrix(const MEl* el);

  GeometryContainer& getGeometry(){ return geometry; }
  node_map& getNodes() { return nodes; }
};

