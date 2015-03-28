#pragma once

#include <memory>
#include <map>
#include "GlobalDefines.h"
#include "armadillo"
//#include "OCC_standard_includes.h"


class BRepAdaptor_Curve;
class BRepAdaptor_Surface;
class MEl;
class myVertex;
class myEdge;
class myFace;

class MNode;

class MNode{ // Abstract base class
 private:

 protected:
  double xyz[3];
  const gind  nd;
  //double* temp;



 public:
 MNode(const int nd_t, double xyz_t[3]): 
  nd(nd_t), xyz{xyz_t[0],xyz_t[1],xyz_t[2]}{}


  virtual inline const int getType() const { return 3; }
  virtual const short int getGeoEntity() const {};
  virtual const double* getParametricCoords() const { return xyz; }


  virtual void setParametricCoords(double* par){ 
    for(int i=0; i<3; i++) xyz[i] = par[i];
  }
  void addElement(MEl *element);

  void eraseElement(MEl *element);


  void setXYZ(const double *xyz_t) {
    for(int i=0; i<3; i++) xyz[i] = xyz_t[i];
  }

  arma::vec getXYZ() const {
    arma::vec3 temp;
    for(int i=0; i<3; i++) temp[i] = xyz[i];
    return temp;
    //return arma::vec(xyz,3);
  }
  
  double* xyzptr() { return xyz; }

  const arma::vec3 xyzvec3() const{ return arma::vec3(xyz); }


  void optimizeSubmesh(const double distortion_threshold);
  void revertNode();
  void saveParams();
  //void setIsStatic(bool s){ _isStatic = s;}
  double computeSubmeshMerit(double& minmerit);
  double getSubmeshMerit(double& minmerit);


  int writeTec(std::ofstream &);
  const gind getND() const {return nd;}
 


};

class MGeoNode: public MNode{
 protected:
  const short int geo_entity;
 
 public:
 MGeoNode(const gind nd_t, double xyz_t[3], 
	  const short int geo_entity_t=-1):
  MNode(nd_t,xyz_t), geo_entity(geo_entity_t){}

  inline const int getType() const { return 0; }
  inline const short int getGeoEntity() const{ return geo_entity; }
  virtual const double* getParametricCoords() const {};

  virtual void setParametricCoords(double* par){}

};


class MEdgeNode: public MGeoNode{
 protected:
  double u;
 
 public:
 MEdgeNode(const int nd_t, double xyz_t[3], 
	   const short int geo_entity_t=-1, double u_t=0):
  MGeoNode(nd_t,xyz_t,geo_entity_t), u(u_t) {}

  inline const double* getParametricCoords() const {return &u; }

  inline const int getType() const { return 1; }
  const double GetU() const { return u; }
  void setParametricCoords(double* par){ u = par[0]; }
  



};

class MFaceNode: public MGeoNode{
 protected:
  double uv[2];

 public:
 MFaceNode(const gind nd_t, double xyz_t[3], 
	   const short int geo_entity_t=-1, double u=0,double v=0):
  MGeoNode(nd_t,xyz_t,geo_entity_t), uv{u,v} {}

  inline const int getType() const { return 2; }
  const double* GetUV() const{ return uv; }
  inline const double* getParametricCoords() const {return uv; }

  void setParametricCoords(double* par){ 
    for(int i=0; i<2; i++) uv[i] = par[i];
  }

};

/*
class MVertNode;

class MVertNode: public MNode{

 private:

 protected:

  const myVertex* _vertex;
 public:
 

 MVertNode(int nd, double xyz[3], int dim, int ID, 
	  const myVertex* vertex): MNode(nd,xyz,dim,ID), 
    _vertex(vertex) {
    } 

  const myVertex* getVertex() const{
    return _vertex;
  }

 





};


class MEdgeNode;

class MEdgeNode: public MNode{


 protected:
  
  const myEdge* _edge;
  double _u;
 
 public:
 
  MEdgeNode(int nd, double xyz[3], int dim, int ID, double u, const myEdge* edge):
  MNode(nd,xyz,dim,ID), _edge(edge), _u(u){}

  const myEdge*  getEdge() const{
    return _edge;
  }

  arma::vec3 param2xyz(double s);

  double getParam() const {return _u; }


 

};

class MFaceNode;

class MFaceNode: public MNode{
 protected:
  const myFace* _face;
  double _uv[2];
 public:

 MFaceNode(int nd, double xyz[3], int dim, int ID, double uv[2], const myFace* face):
  MNode(nd,xyz,dim,ID), _face(face){
    for(int i=0; i<2; i++) _uv[i] = uv[i];
  }

  const myFace* getFace() const{
    return _face;
  }

  arma::vec3 param2xyz(double u, double v);

  //void getParams(double &u, double&v) const {u = _uv[0]; v = _uv[1]; }
  arma::vec2 getParams() const {return arma::vec(_uv,2); }


  std::shared_ptr<MNode> copyNode(){
    return std::shared_ptr<MFaceNode>(new MFaceNode(*this));
  }

  
};

class MVolumeNode: public MNode{
 protected:
 public:

 MVolumeNode(int nd, double xyz[3], int ID):
  MNode(nd,xyz,3,ID){}

 

  
};
*/
