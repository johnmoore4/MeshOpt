#pragma once
#include "Mel.h"

class El2D: public MEl{

 protected:
 
 public:


 El2D(std::vector<gind>& nodes, unsigned char order, unsigned char bc_tag_t,
      short int geo_entity_t=-1, bool geo_type_t=0):
  MEl(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){}

  //virtual const gind* getCornerNodes() const = 0;
  virtual const unsigned char* getSortedCornerNodes() const = 0;
  virtual const int numCornerNodes() const = 0;
  virtual const int getElementType() const = 0;
  virtual const MEl* getChild(const int c) const = 0;
  virtual const signed char getChildOrientation(const int c) const = 0;

  virtual void setChild(const int c, MEl* child, signed char orn=1) = 0;
  virtual void getCornerNodes(int* cnodes) const = 0;
  virtual const int NumNodes() const = 0;
  virtual const int NumChildren() const  = 0;

  inline const int getDim() const{ return 2; }
  inline const bool hasPhysical() const {return true; }
  //const bool is_curved() const {return false; }



};

class MTri: public El2D{
 protected:
  MEl* children[3];
  //const gind nodes[3];

  //unsigned char sorted_nodes[3];
  unsigned char child_orientation[3];


 public:
 MTri(std::vector<gind>& nodes, const NodeIndexer* indexer, unsigned char order, 
      unsigned char bc_tag_t, short int geo_entity_t=-1, bool geo_type_t=0): 
  El2D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }

  void getCornerNodes(int* cnodes) const{
    cnodes[0] = 0;
    cnodes[1] = porder;
    cnodes[2] = (porder+1)*(porder+2)/2-1;
  }

  //inline const gind* getCornerNodes() const { return nodes; }
  //inline const unsigned char* getSortedCornerNodes() const{ return sorted_nodes; }
  inline const unsigned char* getSortedCornerNodes() const{}

  //inline const unsigned char* getSortedCornerNodes() const{}
  inline const int NumChildren() const { return 3; }
  inline const int NumNodes() const { return (porder+1)*(porder+2)/2;}
  inline const int numCornerNodes() const {return 3; }
  inline const int getElementType() const{ return 2; }
  inline void setChild(const int c, MEl* child, signed char orn){ 
    children[c] = child; 
    child_orientation[c] = orn;
  }
  //inline void setChild(const int c, MEl* child){}
  inline const MEl* getChild(const int c) const{ return children[c]; }
  inline const signed char getChildOrientation(const int c) const{
    return child_orientation[c];
  }
};

class MQuad: public El2D{
 protected:	
  //gind children[4];
  //unsigned char sorted_nodes[4];
  signed char child_orientation[4];
  MEl* children[4];
  //const gind nodes[4];


 public:
 MQuad(std::vector<gind>& nodes, const NodeIndexer* indexer, 
       unsigned char order, unsigned char bc_tag_t,
       short int geo_entity_t=-1, bool geo_type_t=0): 
  El2D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }
 
  void getCornerNodes(int* cnodes) const{
    const int ndof = pow((porder+1),2);
    cnodes[0] = 0;
    cnodes[1] = porder;
    cnodes[2] = ndof-porder-1;
    cnodes[3] = ndof-1;
  }

  //inline const gind* getCornerNodes() const { return nodes; }
  //inline const unsigned char* getSortedCornerNodes() const{ return sorted_nodes; }
  inline const unsigned char* getSortedCornerNodes() const{}

  //inline const unsigned char* getSortedCornerNodes() const{}
  inline const int NumChildren() const { return 4; }
  inline const int NumNodes() const { return (porder+1)*(porder+1);}
  inline const int numCornerNodes() const {return 4; }
  inline const int getElementType() const{ return 3; }
   inline void setChild(const int c, MEl* child, signed char orn){ 
    children[c] = child; 
    child_orientation[c] = orn;
  }
  //inline void setChild(const int c, MEl* child){}
  inline const MEl* getChild(const int c) const{ return children[c]; }
  inline const signed char getChildOrientation(const int c) const{
    return child_orientation[c];
  }
};
