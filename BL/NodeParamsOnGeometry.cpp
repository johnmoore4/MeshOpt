#include "BoundaryLayer.h"
#include "assert.h"
#include "Mnode.h"
#include "GeometryContainer.h"
#include <assert.h>

short int BoundaryLayerGenerator::
NodeParamsOnGeometry(const MNode* node, const arma::vec3& normal, 
		     const double dist,
		     arma::vec& params,
		     arma::vec3& xyz){

 
  std::vector<myVertex>& vertices = geometry.getVerticesNC();
  std::vector<myEdge>& edges = geometry.getEdgesNC();
  std::vector<myFace>& faces = geometry.getFacesNC();
  const std::map<int,std::set<int> >& Vertex2EdgesMap = 
    geometry.getVertex2EdgesMap();

  const std::map<int,std::set<int> >& Vertex2FacesMap = 
    geometry.getVertex2FacesMap();

  const std::map<int,std::set<int> >& Edge2FacesMap = 
    geometry.getEdge2FacesMap();

  const int dim = node->getType();
  //const short int entity = node->getGeoEntity();
  const short int initial_entity = node->getGeoEntity();
  short int entity=-1;
 

  //arma::vec pcoords = node->getParametricCoordsVec();
  const double* pcoords = node->getParametricCoords();
  double geo_pcoords[2];
  Geo* phys_entity, *to_param;

  
  if(dim == 0){
    to_param = &vertices[initial_entity];
   
    const std::set<int>& vert_edges = Vertex2EdgesMap.at(initial_entity);
    const std::set<int>& vert_faces = Vertex2FacesMap.at(initial_entity);
    bool face_found = false;

    for(auto it = vert_faces.begin(); it != vert_faces.end(); ++it){
      const myFace& face = faces[*it];
      if(face.getBCTag() < 10 || face.getBCTag() >= 20){
	const int fc = *it;
	phys_entity = &faces[fc];
	entity = fc;
	face_found = true;
	break;
      }
    }
    if(!face_found){
      for(auto it = vert_edges.begin(); it != vert_edges.end(); ++it){
	const myEdge& edge = edges[*it];
	if(edge.getBCTag() < 10 || edge.getBCTag() >= 20){
	  phys_entity = &edges[*it];
	  entity = *it;
	  break;
	} 
      } 
    }  
 
  }
  else if(dim == 1){
 
    to_param = &edges[initial_entity];
 
    const std::set<int>& edge_faces = Edge2FacesMap.at(initial_entity);
    for(auto fc = edge_faces.begin(); fc != edge_faces.end(); ++fc){
      const myFace& face = faces[*fc];

      if(face.getBCTag() < 10 || face.getBCTag() >= 20){
	entity = *fc;
	phys_entity = &faces[*fc];
	break;
      } 
    }

 
    //myFace& face = faces[entity];
    
  }
 
  // Make sure entity is found
  if(!phys_entity){
    cout << "Physical entity not found!" << endl;
  }
  assert(phys_entity);
  const int node_dim = phys_entity->getDim();

 
  phys_entity->paramsOnGeo(*to_param,pcoords,geo_pcoords);
 

  const arma::mat& grad = phys_entity->D1(geo_pcoords);
  //cout << "grad*normal: " << grad*normal << endl;
  
  params = arma::vec(geo_pcoords,node_dim) - dist*grad*normal;
  phys_entity->param2xyz(params.memptr(),xyz.memptr());
 


  return entity;

}
