#pragma once
#include <vector>
#include <memory>
#include "Geo.h"
#include "OCC_standard_includes.h"
#include "armadillo"

class myEdge;
class myFace;
class myVertex;

typedef std::shared_ptr<myEdge> edge_ptr;
typedef std::shared_ptr<myFace> face_ptr;
typedef std::shared_ptr<myVertex> vertex_ptr;

class myEdge: public Geo{
 private:
  const Handle(Geom_Curve) _geom_curve;
  const TopoDS_Edge edge;
  const BRepAdaptor_Curve curve;
  ShapeAnalysis_Curve shape_analysis_curve;
 

  arma::rowvec3 D1_mat;
  double umin, umax;
 public:  
 


 myEdge(const TopoDS_Edge edge):
  edge(edge), curve(edge), _geom_curve(BRep_Tool::Curve(edge,umin,umax)){
    //std::cout << "umin: " << umin << " umax: " << umax << std::endl;
  }
 
 
  double edgeParamsFromPoint(const double *pt0, const double guess) const;


  double paramsOnEdge(const myVertex* vertex) const;

  //const arma::vec3 param2xyz(double s) const;

  const TopoDS_Edge& getTopoDS_Edge() const {return edge;}

  const BRepAdaptor_Curve& getBRepAdaptor_Curve() const{ return curve;}

  const myFace* findSharedFace(const myEdge* edge) const;

  void D1(const double* u, arma::rowvec3& deriv);
  const arma::mat& D1(const double* u);

  const int getDim() const { return 1; }
  void paramsOnGeo(const Geo& geom, const double* geo_params, 
		   double* eval_params);

  //void param2xyz(const double *params, double* u) const;

  const gp_Pnt Value(const double* par) const{ return curve.Value(*par); }

 bool areParamsValid(const double* params);

};
