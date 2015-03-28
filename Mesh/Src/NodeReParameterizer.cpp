#include "NodeReParameterizer.h"
#include "Mnode.h"
#include "GeometryContainer.h"
#include "myFace.h"

NodeReParameterizer::NodeReParameterizer(GeometryContainer& geometry):
  geometry(geometry){

}

int NodeReParameterizer::ReParameterizeOnFace(const int face_id, 
					      const MNode* node,
					      arma::vec3& xyz,
					      arma::vec2& params){

  const myFace& face = geometry.getFaces()[face_id];
  int node_dim = node->getType();
  int node_geo_entity = node->getGeoEntity();
  if(node_dim == 0){
    const myVertex& vertex = geometry.getVertices()[node_geo_entity]; 
    params = face.paramsOnFace(&vertex);
  }
  else if(node_dim == 1){
    const myEdge& edge = geometry.getEdges()[node_geo_entity];
    params = face.paramsOnFace(&edge,node->getParametricCoords()[0]);
  }
  else{
    throw std::logic_error("Cannot re-parametrize a face node on a "
			   "different face!");
  }

  face.param2xyz(params.memptr(),xyz.memptr());

  return 0;
}
