#pragma once
#include <armadillo>

class GeometryContainer;
class MNode;

class NodeReParameterizer{
 public:
  NodeReParameterizer(GeometryContainer& geometry);

  int ReParameterizeOnFace(const int face, 
			   const MNode* node,
			   arma::vec3& xyz,
			   arma::vec2& params);

		       
 private:
  
  GeometryContainer& geometry;
};
