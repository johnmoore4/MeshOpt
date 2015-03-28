#pragma once
#include "GlobalDefines.h"
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include "Mel.h"
//#include "OptEl.h"
#include "MeshTypedefs.h"




arma::vec3 normalFromElnodes(const gind* elnodes, const node_map& nodes);

class MeshContainer{
 private:
  int mesh_dim;
  int symmetry_surface;
  bool is_planar;
  bool has_connectivity;
  
  arma::vec3 mesh_normal;
  //ElementHasher hasher(1);

 protected:
  //element_set elements;
  element_set elements;
  element_set subelements, subsubelements;

  //element_map elements, subelements, subsubelements;

  //inverse_element_map inv_elements, inv_subelements, inv_subsubelements;

  node_map nodes;

  //ideal_map idealElements;

  //std::map<MEl*,double> all_merits;

 public:
  /*
 MeshContainer(): elements(1,ElementHasher(),ElementEqual(1)){
    //ElementEqual equal(1);
  }
  */

  element_set& getElementsOfDim(int dim);

  const element_set& getElements() const {return elements; }
  const element_set& getSubElements() const { return subelements; }
  const element_set& getSubSubElements() const { return subsubelements; }
  const node_map& getNodes() const { return nodes; }
  //const ideal_map& getIdealElements() const{ return idealElements; }
  //const ElDouble_map& getAllMerits(){ return all_merits; }
  element_set& getElementsNC() {return elements; }
  element_set& getSubElementsNC() { return subelements; }
  element_set& getSubSubElementsNC() { return subsubelements; }
  node_map& getNodesNC() { return nodes; }
  //ideal_map& getIdealElementsNC() { return idealElements; }
  //ElDouble_map& getAllMeritsNC(){ return all_merits; }

  int Dimension() const{ return mesh_dim; }
  const int MeshDimension() const { return mesh_dim; } 
  const bool IsMeshPlanar() const { return is_planar; }
  const arma::vec3& getMeshNormal() const{ return mesh_normal; }
  const int SymmetrySurface() const{ return symmetry_surface; }
  bool hasConnectivity() const{ return has_connectivity; }

  void SetMeshDimension(const int dim){ mesh_dim = dim; }
  void computeIsMeshPlanar();
  void SetSymmetrySurface(int surf){ symmetry_surface = surf; }
  void setHasConnectivity(bool has){ has_connectivity = has; }


};
