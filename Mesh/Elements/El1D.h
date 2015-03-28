#pragma once
#include "Mel.h"

class El1D: public MEl{
 protected:
 
  
 public:
 
 El1D(std::vector<gind>& nodes, unsigned char order, unsigned char bc_tag_t,
      short int geo_entity_t=-1, bool geo_type_t=0):
  MEl(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){}

  //inline const bool is_curved() const { return false; }
  inline const int getDim() const{ return 1; }  

  //virtual const gind* getCornerNodes() const = 0;
  virtual const unsigned char* getSortedCornerNodes() const = 0;
  virtual const int numCornerNodes() const = 0;
  virtual const int getElementType() const = 0;
  virtual const MEl* getChild(const int c) const = 0;
  virtual const signed char getChildOrientation(const int c) const = 0;
  virtual void getCornerNodes(int* cnodes) const = 0;
  virtual const int NumNodes() const = 0;
  virtual const int NumChildren() const  = 0;
  void setChild(const int c, MEl* child, signed char orn){};

};

class MLine: public El1D{
 private:
  //unsigned char sorted_nodes[2];
  //const gind nodes[2];


 public:
 MLine(std::vector<gind>& nodes, const NodeIndexer* indexer, unsigned char order, 
       unsigned char bc_tag_t, short int geo_entity_t=-1, bool geo_type_t=0): 
  El1D(nodes,order,bc_tag_t,geo_entity_t,geo_type_t){
    //sortNodes(sorted_nodes,indexer);
  }

  void getCornerNodes(int* cnodes) const{
    cnodes[0] = 0;
    cnodes[1] = porder;
  }

  //inline const gind* getCornerNodes() const { return nodes; }
  //inline const unsigned char* getSortedCornerNodes() const{ return sorted_nodes; }
  inline const unsigned char* getSortedCornerNodes() const{}

  //inline const unsigned char* getSortedCornerNodes() const{}
  inline const int NumChildren() const { return 2; }
  inline const int NumNodes() const { return porder+1;}
  inline const int numCornerNodes() const {return 2; }
  inline const int getElementType() const{ return 1; }
  //inline void setChild(const int c, MEl* child) {};
  inline const MEl* getChild(const int c) const{ return NULL; }
  const signed char getChildOrientation(const int c) const { return 0; }
  
};
