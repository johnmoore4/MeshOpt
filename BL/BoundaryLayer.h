#pragma once
#include "GlobalDefines.h"
//#include <map>
#include "armadillo"
#include "NodeIndexer.h"
#include <map>
#include <unordered_set>
#include "MeshContainer.h"

#include <vector>
#include <memory>

//class MeshContainer;
class GeometryContainer;
class ShapeFunctionMatricesFactory;
class MEl;
class MNode;
class BLParameterList;

class BLNodeInfo{

public:
 BLNodeInfo(arma::vec3 normal_t, int dim, double dist_t): 
  normal(normal_t), node_dim(dim), dist(dist_t){
    on_symmetry = false; }
  arma::vec3 normal;
  int node_dim;
  bool on_symmetry;
  int symmetry_entity;
  double dist;

};
typedef std::map<gind,BLNodeInfo> normal_map_type;
typedef std::map<gind,gind> gindmap;

class BoundaryLayerGenerator{


 private:
  MeshContainer& mesh;
  GeometryContainer& geometry;
  ShapeFunctionMatricesFactory& sf_factory;
  std::shared_ptr<BLParameterList>& bl_parameters;

  std::map<const MEl*,arma::mat> ideal_elements;
  std::set<unsigned short int> symmetry_faces;

  int NNorig;

  NodeIndexFactory index_factory;
  gindmap newnode_map;
  normal_map_type normal_map;
  std::map<MEl*,bool> ExtrudedIdeal;

  std::vector<unique_element_ptr> store_symmetry;

  std::vector<MEl*> FindBLElements(const element_set& elements);
  std::vector<MEl*> MakeSymmetryLineElements();

  int findSymmetryFaces();

  std::vector<double> getMaximumSafeExtrusionDistance(int safety_factor);
  
  normal_map_type CreateNodeInfoMap(const std::vector<MEl*>& bl_elements,
				    double thickness);

  gindmap CreateNewNodeMap(const normal_map_type& normal_map);

  void ReplaceNeighboringElements(element_set& elements,
				  const int Nlayers,
				  const int surf_to_restrict=-1);

  void ReplaceSymmetryElements(const std::map<gind,gind>& newnode_map,
			       const normal_map_type& normal_map,
			       const int Nlayers);
  /*
  void InsertBLNodes(const normal_map_type& normal_map,
		     const int Nlayers, const double thickness, 
		     const double ratio=-1,
		     const std::vector<double>& max_extrude=
		     std::vector<double>());
  */

  int SmoothField(double Nsmooth, const std::vector<MEl*>& bl_elements,
		  std::vector<double>& field);

  void InsertBLNodes(const double ratio, 
		     const std::vector<double>& extrusion_dist);

  void ResetBLNodes(const double ratio, 
		     const std::vector<double>& extrusion_dist);

  void InsertBLElements(element_set& elements,
			const std::vector<MEl*>& bl_elements,
			const int Nlayers,
			const double thickness,
			const bool copy_geo_entity = false);
  /*
  void ResetBLNodes(const normal_map_type& normal_map,
		    const std::map<gind,gind>& newnode_map,
		    const int Nlayers, const double thickness, 
		    const double factor);
  */

  //short int NodeParamsOnGeometry(const MNode* node, arma::vec& params);
  short int NodeParamsOnGeometry(const MNode* node, const arma::vec3& normal, 
				 const double dist,
				 arma::vec& params,
				 arma::vec3& xyz);
 protected:

 public:
 BoundaryLayerGenerator(MeshContainer& mesh_t, GeometryContainer& geometry_t,
			ShapeFunctionMatricesFactory& sf_factory_t,
			std::shared_ptr<BLParameterList>& bl_parameters): 
  mesh(mesh_t), geometry(geometry_t), sf_factory(sf_factory_t),
    bl_parameters(bl_parameters){}

  int GenerateBL();
  int OptimizeBL();
  int OptimizeBLSubmesh();
  int SubDivideBL();
  
  std::map<const MEl*,arma::mat>& getIdealElements(){ return ideal_elements; }

  std::unordered_set<int> clamped_nodes;

  const normal_map_type& getNormalMap() { return normal_map; }
  

};
