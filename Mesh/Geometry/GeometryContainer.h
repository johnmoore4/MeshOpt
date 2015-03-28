#pragma once
#include "myFace.h"
#include "myEdge.h"
#include "myVertex.h"
#include "OCC_standard_includes.h"
#include <set>

class GeometryContainer{
 private:

 protected:
  std::map<int,std::set<int> > Edge2FacesMap;
  std::map<int,std::set<int> > Vertex2FacesMap;
  std::map<int,std::set<int> > Vertex2EdgesMap;

  std::vector<myFace > myFaces;
  std::vector<myEdge> myEdges;
  std::vector<myVertex> myVertices;
 public:  
  TopoDS_Shape shape;
  TopTools_IndexedMapOfShape fmap, emap, vmap;
  const std::vector<myFace>& getFaces() { return myFaces; }
  const std::vector<myEdge>& getEdges() { return myEdges; }
  const std::vector<myVertex>& getVertices() { return myVertices; }
  const std::map<int,std::set<int> >& getEdge2FacesMap() const{ 
    return Edge2FacesMap; 
  }
  const std::map<int,std::set<int> >& getVertex2FacesMap() const{ 
    return Vertex2FacesMap; 
  }
  const std::map<int,std::set<int> >& getVertex2EdgesMap() const{
    return Vertex2EdgesMap;
  }
  std::vector<myFace>& getFacesNC() { return myFaces; }
  std::vector<myEdge>& getEdgesNC() { return myEdges; }
  std::vector<myVertex>& getVerticesNC() { return myVertices; }
  std::map<int,std::set<int> >& getEdge2FacesMapNC() { return Edge2FacesMap; }
  std::map<int,std::set<int> >& getVertex2FacesMapNC() { 
    return Vertex2FacesMap; 
  }
  std::map<int,std::set<int> >& getVertex2EdgesMap(){
    return Vertex2EdgesMap;
  }

};

