#pragma once
#include "MeshTypedefs.h"
#include <armadillo>

class GeometryContainer;

class ParametricCoordinateEvaluator{
 public:
  ParametricCoordinateEvaluator(GeometryContainer& geometry,
				const node_map& nodes);

  arma::mat getParametricCoordinates(const unique_element_ptr& el);
  arma::mat projectPointsOnGeometry(const unique_element_ptr& el,
				    const arma::mat& nodes_guess,
				    const arma::mat& parametric_guess,
				    arma::mat& parametric_coords);

 private:
  GeometryContainer& geometry;
  const node_map& nodes;

};
