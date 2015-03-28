#pragma once

#include "OCC_standard_includes.h"
#include <string>
#include <vector>
#include <memory>
#include "armadillo"

class myEdge;
class myFace;
class myVertex;

class OCC_Handler{
 
 protected:
  // The OCC geometry
  TopoDS_Shape _shape;
 
 public:
  typedef std::shared_ptr<myEdge> edge_ptr;
  typedef std::shared_ptr<myFace> face_ptr;
  typedef std::shared_ptr<myVertex> vertex_ptr;

  TopTools_IndexedDataMapOfShapeListOfShape 
    EdgeFaceMap, VertexEdgeMap, VertexFaceMap;

  TopTools_IndexedMapOfShape fmap, emap, vmap, somap, shmap, wmap;
  std::vector< BRepAdaptor_Curve > curves;
  std::vector< BRepAdaptor_Surface > surfaces;
  std::vector< TopoDS_Edge > edges;
  std::vector< TopoDS_Face > faces;
  //std::vector< TopoDS_Vertex > vertices;
  std::vector< face_ptr > myFaces;
  std::vector< edge_ptr> myEdges;
  std::vector< vertex_ptr> myVertices;


  void readSTEP(const std::string name);
  void makeLists();
 
  
};
