#pragma once
#include "Mel.h"

class El3D: public MEl{

 protected:

 public:
 El3D(std::vector<gind>& nodes, unsigned char order, unsigned char bc_tag_t,
      short int geo_entity_t=-1,bool geo_type_t=0):
  MEl(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){}



  //virtual const gind* getCornerNodes() const = 0;
  virtual const unsigned char* getSortedCornerNodes() const = 0;
  virtual const int numCornerNodes() const = 0;
  virtual const int getElementType() const = 0;
  virtual const MEl* getChild(const int c) const = 0;
  virtual void setChild(const int c, MEl* child, signed char orn=1) = 0;
  virtual const signed char getChildOrientation(const int c) const = 0;
  virtual void getCornerNodes(int* cnodes) const = 0;
  virtual const int NumNodes() const = 0; 
  virtual const int NumChildren() const  = 0;

  inline const bool hasPhysical() const {return false; }
  inline const int getDim() const{ return 3; }
  //inline const bool is_curved() const{};


};

class MTet: public El3D{
 protected:

  MEl* children[4];
  //const gind nodes[4];
  //unsigned char sorted_nodes[4];
  signed char child_orientation[4];

 public:
 MTet(std::vector<gind>& nodes, const NodeIndexer* indexer, 
      unsigned char order, unsigned char bc_tag_t,
      short int geo_entity_t=-1, bool geo_type_t=0): 
  El3D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }

  void getCornerNodes(int* cnodes) const{
    cnodes[0] = 0;
    cnodes[1] = porder;
    cnodes[2] = (porder+1)*(porder+2)/2-1;
    cnodes[3] = (porder+1)*(porder+2)*(porder+3)/6-1;
  }

  //inline const gind* getCornerNodes() const { return nodes; }
  //inline const unsigned char* getSortedCornerNodes()const {}
  inline const unsigned char* getSortedCornerNodes()const {}
  inline const int NumChildren() const { return 4; }
  inline const int NumNodes() const {return (porder+1)*(porder+2)*(porder+3)/6;}
  inline const int numCornerNodes() const {return 4; }
  inline const int getElementType() const{ return 4; }
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

class MHex: public El3D{
 protected:
  MEl* children[6];
  signed char child_orientation[6];

 public:
 MHex(std::vector<gind>& nodes, const NodeIndexer* indexer, 
      unsigned char order, unsigned char bc_tag_t,
      short int geo_entity_t=-1, bool geo_type_t=0): 
  El3D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }

  void getCornerNodes(int* cnodes) const{
    int ndofq = pow(porder+1,2);
    cnodes[0] = 0;
    cnodes[1] = porder;
    cnodes[2] = ndofq-porder-1;
    cnodes[3] = ndofq-1;
    for(int i = 0; i < 4; i++){
      cnodes[4+i] = cnodes[i] + porder*ndofq;
    }

  }

  //inline const gind* getCornerNodes() const { return nodes; }
  inline const unsigned char* getSortedCornerNodes()const {}
  inline const int NumChildren() const { return 6; }
  inline const int NumNodes() const { return pow(porder+1,3); }
  inline const int numCornerNodes() const {return 8; }
  inline const int getElementType() const{ return 5; }
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

class MPrism: public El3D{
 protected:
  MEl* children[5];
  //const gind nodes[6];
  //unsigned char sorted_nodes[6];
  signed char child_orientation[5];

 public:
 MPrism(std::vector<gind>& nodes, const NodeIndexer* indexer, 
	unsigned char order, unsigned char bc_tag_t,
	short int geo_entity_t=-1, bool geo_type_t=0): 
  El3D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }

  void getCornerNodes(int* cnodes) const{
    const int ndoftri = (porder+1)*(porder+2)/2;
    cnodes[0] = 0;
    cnodes[1] = porder;
    cnodes[2] = ndoftri-1;
    cnodes[3] = porder*ndoftri;
    cnodes[4] = porder*ndoftri+porder;
    cnodes[5] = (porder+1)*ndoftri-1;
  }

  //inline const gind* getCornerNodes() const { return nodes; }
  inline const unsigned char* getSortedCornerNodes()const {}
  inline const int NumChildren() const { return 5; }
  inline const int NumNodes() const {return (porder+1)*(porder+2)/2*(porder+1);}
  inline const int numCornerNodes() const {return 6; }
  inline const int getElementType() const{ return 6; }
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


class MPyramid: public El3D{
 protected:
  MEl* children[5];
  signed char child_orientation[5];

 public:
 MPyramid(std::vector<gind>& nodes, const NodeIndexer* indexer, 
      unsigned char order, unsigned char bc_tag_t,
      short int geo_entity_t=-1, bool geo_type_t=0): 
  El3D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }

  void getCornerNodes(int* cnodes) const{
    int ndofq = pow(porder+1,2);
    cnodes[0] = 0;
    cnodes[1] = porder;
    cnodes[2] = ndofq-porder-1;
    cnodes[3] = ndofq-1;
    cnodes[4] = NumNodes()-1;
    

  }

  //inline const gind* getCornerNodes() const { return nodes; }
  inline const unsigned char* getSortedCornerNodes()const {}
  inline const int NumChildren() const { return 5; }
  inline const int NumNodes() const { return (porder+1)*(porder+2)*(2*(porder+1)+1)/6; }
  inline const int numCornerNodes() const {return 5; }
  inline const int getElementType() const{ return 7; }
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
