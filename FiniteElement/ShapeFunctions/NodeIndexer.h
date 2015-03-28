#pragma once
#include <memory>
#include <vector>

typedef unsigned int indtype;
class NodeIndexer{
 private:

 protected:
  int order;
  int gmsh_type;
  int my_type=0;

  int ComputeLocalFaceIndices();

  std::vector<std::vector<indtype> > child_nodes;
  std::vector<std::vector<indtype> > oriented_nodes;
  std::vector<indtype> corner_nodes;
  std::vector<indtype> reversed_nodes;
  std::vector<indtype> interior_nodes;
  std::vector<indtype> gmsh_index;
  std::vector<int> gmsh_types;
  std::vector<indtype> gmsh_corner_nodes;
  std::vector<std::vector<indtype> > local_face_indices;
  std::vector<indtype> rev_local_face_indices;
 public:  
 
  int Initialize(int p);

  virtual const int Ndof(const int p) const = 0;
  virtual const int NumChildren() const = 0;
  virtual const int NumBoundary() const = 0;
  virtual const int* ChildCornerNodes(const int child) const = 0;
  virtual const int NumChildCornerNodes(const int child)const = 0;
  virtual const int ChildType(const int child) const = 0;
  virtual const int Dimension() const = 0;
  virtual void SetOrder(const int p) = 0;

  const int getP() const{ return order; }
  const int getGMSHType() const{ return gmsh_type; }
  const int getMyType() const{ return my_type; }
  const int Ndof() const { return Ndof(getP()); }
  const std::vector<indtype>& getChildNodes(const int child) const{ 
    return child_nodes[child];
  } 

  const std::vector<indtype> getOrientedNodes(const int orient) const{ 
    if(orient > 0){
      return oriented_nodes[orient-1];
    }
    else{
      const std::vector<indtype>& reversed = getReversedNodes();
      std::vector<indtype> temp(reversed.size());
      for(int i=0; i<reversed.size(); i++){
	temp[i] = oriented_nodes[-orient-1][reversed[i]];
      }
      return temp;
    }
  } 

  const std::vector<indtype> getExteriorNodes() const;


  const std::vector<indtype>& getReversedNodes() const{ return reversed_nodes;}
  const std::vector<indtype>& getInteriorNodes() const{ return interior_nodes;}
  const std::vector<indtype>& getCornerNodes() const{ return corner_nodes; }
  const std::vector<indtype>& getGMSHIndex() const{ return gmsh_index; }
  const std::vector<int>& getGMSHTypes() const { return gmsh_types; }
  const std::vector<indtype>& gmshCornerNodes() const{
    return gmsh_corner_nodes;
  }
  const std::vector<std::vector<indtype> >& getLocalFaceIndices() const{
    return local_face_indices;
  } 
  const std::vector<indtype>& getRevLocalFaceIndices() const{
    return rev_local_face_indices;
  }
};


class PointNodeIndexer: public NodeIndexer{
 private:

 public:
  PointNodeIndexer(int p){ 
    //SetOrder(p);
  }
  inline const int Dimension() const{ return 0; }
  inline const int Ndof(const int p) const { return 1; }
  inline const int NumChildren() const { return 0; }
  inline const int NumBoundary() const { return 0; }
  inline const int NumChildCornerNodes(const int child) const { return 0; }
  inline const int ChildType(const int child) const { return 0; }
  inline const int* ChildCornerNodes(const int child) const { }
  void SetOrder(int p);
};

class LineNodeIndexer: public NodeIndexer{
 private:
  const int index[2] = {0,1};
 public:
  LineNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
    //porder = p;
  }
  inline const int Dimension() const{ return 1; }
  inline const int Ndof(const int p) const { return p+1; }
  inline const int NumChildren() const { return 1; }
  inline const int NumBoundary() const { return 1; }
  inline const int NumChildCornerNodes(const int child) const { return 2; }
  inline const int ChildType(const int child) const { return 1; }
  inline const int* ChildCornerNodes(const int child) const { return index; }
 void SetOrder(const int p);
};

class TriNodeIndexer: public NodeIndexer{
 private:
  const int index[3][2] = {{0,1},{1,2},{2,0}};
 public:
  TriNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
  }
  inline const int Dimension() const{ return 2; }
  inline const int Ndof(const int p) const { return (p+1)*(p+2)/2; }
  inline const int NumChildren() const { return 3; }
  inline const int NumBoundary() const { return 3; }
  inline const int NumChildCornerNodes(const int child )const { return 2; } 
  inline const int ChildType(const int child) const { return 1; } 
  inline const int* ChildCornerNodes(const int child) const {
    return index[child];
  }
  void SetOrder(const int p);
};

class QuadNodeIndexer: public NodeIndexer{
 private:
  //const int index[4][2] = {{0,1},{1,2},{2,3},{3,0}};
  const int index[4][2] = {{0,1},{1,3},{3,2},{2,0}};
 public:
  QuadNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
  }
  inline const int Dimension() const{ return 2; }
  inline const int Ndof(const int p) const { return (p+1)*(p+1); }
  inline const int NumChildren() const { return 4; }
  inline const int NumBoundary() const { return 4; }
  inline const int NumChildCornerNodes(const int child ) const { return 2; }
  inline const int ChildType(const int child) const { return 1; } 
  inline const int* ChildCornerNodes(const int child) const { 
    return index[child]; 
  }
  void SetOrder(const int p);


};

class TetNodeIndexer: public NodeIndexer{
 private:
  const int index[4][3] = {{0,2,1},{0,1,3},{1,2,3},{0,3,2}};
  
 public:
  TetNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
  }
  inline const int Dimension() const{ return 3; }
  inline const int Ndof(const int p) const { return (p+1)*(p+2)*(p+3)/6; }
  inline const int NumChildren() const { return 4; };
  inline const int NumBoundary() const { return 0; }
  inline const int NumChildCornerNodes(const int child) const { return 3; }
  inline const int ChildType(const int child) const { return 2; } 
  inline const int* ChildCornerNodes(const int child) const{ 
    return index[child]; 
  }
  void SetOrder(const int p);
};

class PrismNodeIndexer: public NodeIndexer{
 private:
 
  const int indextri [2][3] = {{0,2,1},{3,4,5}};
  const int indexhex [3][4] = {{0,1,4,3},{1,2,5,4},{2,0,3,5}};

 public:
  PrismNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
  }
  inline const int Dimension() const{ return 3; }
  inline const int Ndof(const int p) const { return (p+1)*(p+2)/2*(p+1); }
  inline const int NumChildren() const { return 5; };
  inline const int NumBoundary() const { return 0; }
  inline const int ChildType(const int child) const{
    if(child < 2) return 2;
    else return 3;
  }
  inline const int NumChildCornerNodes(const int child) const { 
    if(child < 2) return 3;
    else return 4;
  }
  inline const int* ChildCornerNodes(const int child) const{ 
    if(child < 2) return indextri[child];
    else return indexhex[child-2];
  }
  void SetOrder(const int p);

};

class HexNodeIndexer: public NodeIndexer{
 private:
  const int index [6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},
  			    {3,0,4,7}};
 public:
  HexNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
  }
  inline const int Dimension() const{ return 3; }
  inline const int Ndof(const int p) const { return (p+1)*(p+1)*(p+1); }
  inline const int NumChildren() const { return 6; };
  inline const int NumBoundary() const { return 0; }
  inline const int NumChildCornerNodes(const int child) const { return 4; }
  inline const int ChildType(const int child) const { return 3; } 
  inline const int* ChildCornerNodes(const int child) const{ 
    return index[child]; 
  }
  void SetOrder(const int p);
};

class PyramidNodeIndexer: public NodeIndexer{
 private:
  //const int index [6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},
  //			    {3,0,4,7}};
 public:
  PyramidNodeIndexer(int p){
    SetOrder(p);
    ComputeLocalFaceIndices();
  }
  inline const int Dimension() const{ return 3; }
  inline const int Ndof(const int p) const {
    return (p+1)*(p+2)*(2*(p+1)+1)/6;
  }
  inline const int NumChildren() const { return 5; };
  inline const int NumBoundary() const { return 0; }
  inline const int NumChildCornerNodes(const int child) const { return 5; }
  inline const int ChildType(const int child) const { 
    if(child == 0) return 3;
    else return 2;
  } 
  inline const int* ChildCornerNodes(const int child) const{ 
    //return index[child]; 
    throw std::runtime_error("ChildCornerNodes not supported for Pyramid!");
  }
  void SetOrder(const int p);
};

class NodeIndexFactory{
 private:
  //std::shared_ptr<NodeIndexer> Point, Line, Tri, Quad, Tet, Prism, Hex, Pyramid;

  std::shared_ptr<NodeIndexer> indexers[8][11];
  
  std::vector<std::pair<int,int> > GMSH2Type;

  const std::shared_ptr<NodeIndexer> NodeIndexerFromType(int type, int order);

 public:
  NodeIndexFactory();

  const NodeIndexer*  
    getNodeIndexer(const int type, const int p); 

 

  const int TypeFromGMSHType(const int type);
  const int OrderFromGMSHType(const int type);
};
