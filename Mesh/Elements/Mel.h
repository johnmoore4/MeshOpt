#pragma once
#include <memory>
#include <armadillo>
#include <boost/functional/hash.hpp>
#include <unordered_map>

#include "GlobalDefines.h"

#include "Mnode.h"
#include "NodeIndexer.h"
#include "MeshTypedefs.h"

class myEdge;
class myFace;
typedef std::shared_ptr<myEdge> edge_ptr;
typedef std::shared_ptr<myFace> face_ptr;
typedef std::shared_ptr<MNode> node_ptr;

//class NodeIndexFactory;
//class NodeIndexer;

class MEl;
class El1D;
class El2D;
class El3D;

//typedef std::unique_ptr<MEl> unique_element_ptr;
//typedef std::unique_ptr<MNode> unique_node_ptr;



signed char ElementOrientation(const unique_element_ptr& el1,
			       const unique_element_ptr& el2,
			       NodeIndexFactory& indexfactory);

arma::mat nodeParamsOnFace(const std::vector<node_ptr> &nodes,
			   const myFace* face);

std::vector<node_ptr> projectPointsOnFace(const arma::mat xyz0, 
					  const arma::mat& uv0,
					  const myFace* face);

std::vector<node_ptr> createFreeNodes(arma::mat& xyz);


class Element{

 public:

};

class MEl: public Element{
 
 private:
  
 protected:
  //std::unique_ptr<MEl[]> 
  std::unique_ptr<gind[]> _nodes;
  //std::unique_ptr<std::unique_ptr<MNode>[] > _interior_nodes;
  short int geo_entity;
  unsigned char bc_tag, porder;
  bool curved, geo_type;
  bool is_bl = false;
  bool is_curved = true;

  void sortNodes(unsigned char* sorted_nodes, const NodeIndexer* indexer);
 

 public:
 
  inline const int getGeoEntity() const { return geo_entity; }
  inline const int getBCTag() const { return bc_tag; }
  inline const bool getGeoType() const { return geo_type; }
  inline const bool hasGeoEntity() const{ return (geo_entity != -1); }
  inline const unsigned char getOrder() const{ return porder; }
  inline const bool IsBL() const { return is_bl; }
  inline void setIsBL(bool is) { is_bl = is; }
  inline void setGeoEntity(short int entity) { geo_entity = entity; }
  inline const gind* getNodes() const{ return _nodes.get(); }
  inline gind* getNodesNC(){ return _nodes.get(); }
  inline const gind* getCornerNodes() const{
    return getNodes();
  }
  inline void setIsCurved(bool is){ is_curved = is; }
  inline int isCurved() const { return is_curved; }
  inline void setBCTag(unsigned char bct){ bc_tag = bct; }
  /*
  inline  const std::unique_ptr<std::unique_ptr<MNode>[] >& getInteriorNodes()
  const {
    return _interior_nodes;
  }
  */
 MEl(std::vector<gind>& nodes, unsigned char order, unsigned char bc_tag_t,
     short int geo_entity_t=-1,bool geo_type_t=0);

  
 
 void ComputeSortedCornerNodes(gind* sorted);

  virtual const int getDim() const = 0; 
  //virtual const bool is_curved() const=0;
  virtual const unsigned char* getSortedCornerNodes() const = 0;
  virtual const int getElementType() const = 0;
  virtual void setChild(const int c, MEl* child, signed char orn=1) = 0;
  virtual const MEl* getChild(const int c) const = 0;
  virtual const signed char getChildOrientation(const int c) const = 0;
  //virtual const gind* getCornerNodes() const = 0;
  virtual const int numCornerNodes() const = 0;
  virtual const int NumNodes() const = 0;
  virtual const int NumChildren() const  = 0;

  virtual void getCornerNodes(int* cnodes) const = 0;


  //void setNode(const int n, node_ptr node);
  void setOrder(int p);

  //std::unique_ptr<gind[]> computeSortedNodes()const;



 
};











