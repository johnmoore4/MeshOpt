#include "ParametricCoordinateEvaluator.h"
#include "GeometryContainer.h"
#include "Mel.h"

ParametricCoordinateEvaluator::
ParametricCoordinateEvaluator(GeometryContainer& geometry,
			      const node_map& nodes):
  geometry(geometry), nodes(nodes){

}

arma::mat ParametricCoordinateEvaluator::
getParametricCoordinates(const unique_element_ptr& el){
  unsigned char geo_entity = el->getGeoEntity();
  int geo_dim = el->getGeoType()+1;

  //int dim = el->getDim();
  assert(el->hasGeoEntity());

  //std::cout << "dim/geo entity: " << geo_dim << " " << int(geo_entity) << std::endl;

  int nn = el->NumNodes();
  const gind* elnodes = el->getNodes();
  
  //std::cout << "elnodes: " << std::endl;
  for(int nd = 0; nd < nn; nd++){
    const node_ptr& node = nodes.at(elnodes[nd]);
    //std::cout << node->xyzvec3().t()  << std::endl;

  }
  arma::mat parametric_coords(geo_dim,nn);
  for(int nd = 0; nd < nn; nd++){
    const node_ptr& node = nodes.at(elnodes[nd]);
    int node_geo_dim = node->getType();
    //std::cout << "node geo dim/entity: " << node_geo_dim << " " << node->getGeoEntity() << std::endl;
    //std::cout << node->xyzvec3().t()  << std::endl;
    if(node_geo_dim == geo_dim){
      double* colptr = parametric_coords.colptr(nd);
      for(int d = 0; d < geo_dim; d++){
	colptr[d] = node->getParametricCoords()[d];
      }
    }
    else if(node_geo_dim == 1){
      // re-paramaterize edge node on face
      if(geo_dim == 2){
	const myFace& face = geometry.getFaces()[el->getGeoEntity()];
	const myEdge& edge = geometry.getEdges()[node->getGeoEntity()];
	parametric_coords.unsafe_col(nd) = 
	  face.paramsOnFace(&edge,node->getParametricCoords()[0]);
      }
      else throw std::logic_error("This doesn't compute");
    }
    else if(node_geo_dim == 0){
      const myVertex& vertex = geometry.getVertices()[node->getGeoEntity()];
      // re-parameterize vertex node on face
      if(geo_dim == 2){

	const myFace& face = geometry.getFaces()[el->getGeoEntity()];
	parametric_coords.unsafe_col(nd) = face.paramsOnFace(&vertex);

      }
      // re-parameterize vertex node on line
      else if(geo_dim == 1){
	const myEdge& edge = geometry.getEdges()[el->getGeoEntity()];
	parametric_coords.unsafe_col(nd)[0] = edge.paramsOnEdge(&vertex);
      }
      else throw std::logic_error("This does not compute");
    }
  }

  return parametric_coords;
}
arma::mat ParametricCoordinateEvaluator::
projectPointsOnGeometry(const unique_element_ptr& el,
			const arma::mat& nodes_guess,
			const arma::mat& parametric_guess,
			arma::mat& parametric_coords){
  //std::cout << "b" << std::endl;
  const int neval = nodes_guess.n_cols;
  assert(neval == parametric_guess.n_cols);

  //const int nn = el->NumNodes();
  //const int dim = el->getDim();

  int el_geo_dim = el->getGeoType()+1;

  parametric_coords.resize(el_geo_dim,neval);
  arma::mat projected_nodes(3,neval);
  for(int nd = 0; nd < neval; nd++){

    if(el_geo_dim == 1){
      //std::cout << "nodes_guess: " << std::endl;
      //std::cout << nodes_guess << std::endl;
      const myEdge& edge = geometry.getEdges()[el->getGeoEntity()];
      parametric_coords[nd] = 
	edge.edgeParamsFromPoint(nodes_guess.colptr(nd),parametric_guess[nd]);
      edge.param2xyz(parametric_coords.colptr(nd),projected_nodes.colptr(nd));
      //std::cout << parametric_coords[nd] << std::endl;
    }
    else if(el_geo_dim == 2){
      const myFace& face = geometry.getFaces()[el->getGeoEntity()];
      parametric_coords.unsafe_col(nd) = 
	face.faceParamsFromPoint(nodes_guess.colptr(nd),
				 parametric_guess.colptr(nd));
      face.param2xyz(parametric_coords.colptr(nd),projected_nodes.colptr(nd));
    }
    else throw std::logic_error("Dim must be 1 or 2");

    
  }
  //std::cout << "end" << std::endl;

  return projected_nodes;
}
