#pragma once

#include <string>
#include <vector>
#include <set>
#include <memory>
//#include <fstream>
//#include "Mel.h"
//#include "Mnode.h"

//#include "Mnode.h"

//#include "shapeFunctionHandler.h"
//#include "elementFactory.h"
#include "ShapeFunctionMatrices.h"
#include "GlobalDefines.h"
#include "MeshContainer.h"
#include "GeometryContainer.h"

#include <unordered_set>
#include <map>
#include <functional>



class sf_factory;
class elementFactory;
class nodeFactory;

//class GMSHReader;
class MeshReader;
class GeometryReader;
class MeshOptConfig;
class BoundaryLayerGenerator;

class HModel{
 private:

 protected:

  int _porder;

  //int _numNodesOrig;

  int mesh_dim;

  // Name of input files
  std::string step_file_name;
  std::string gmsh_file_name;


  
  // Shape function factory
  //std::shared_ptr<sf_factory> sf_handler;

  elementFactory *element_factory;

  nodeFactory *node_factory;

  std::unique_ptr<MeshReader> mesh_reader;
  std::unique_ptr<GeometryReader> geometry_reader;
  std::shared_ptr<MeshOptConfig> mesh_opt_config;
  std::shared_ptr<BoundaryLayerGenerator> bl_generator;

  MeshContainer mesh_container;
  GeometryContainer geometry_container;
  ShapeFunctionMatricesFactory sf_factory;
  std::shared_ptr<MeshContainer> high_order_mesh;
  

 public:

  HModel();
  ~HModel();


  int _optimal;

  int parseCommandLineOptions(int argc, char **arcv);
  

  
  //std::unique_ptr<OCC_Handler> occ_handler;

  int readMesh(void);
  int ReadGeometry(void);
  int WriteMesh(std::string filename, std::string format);
  int PrepareMesh();

  //int readGMSH(void);
  //int readGMSHElements(std::ifstream &fname);
  //int readGMSHNodes(std::ifstream &fname);
  //int readOCC(void);

  //int associateVerticesWithGeom(std::unique_ptr<OCC_Handler> &);

  int MeshHighOrder();

  //int meshHighOrder2(int porder, int pfinal);

  //int writeTec(const std::string filename){}

  //int writeGMSH(const std::string filename);

  //int computeMeshStats();

  //int optimizeMesh();

  //double computeGlobalMerit();

  //int computeUpdatedNodes(double step_size, bool recompute);
  
  //int computeUpdateSubmesh(const double distortion_threshold);

  //int revertNodes();

  //int generateHierarchical(const int porder);

  //int getMeshMerits(double &max_merit, double &global_merit);

  //int getMeshConnectivity();

  int generateHOMesh(int porder);
  
  //int debugSF();

  //int totalVolume();

  //int totalSurfaceArea();

  int generateBL();

  //int DebugLinearElasticity(){};


};
