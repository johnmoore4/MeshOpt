#pragma once

#include <vector>
#include <memory>
#include "Geo.h"
#include "OCC_standard_includes.h"

class myEdge;
class myFace;

class myVertex: public Geo{

 private:
  const TopoDS_Vertex vertex;
  const gp_Pnt point;

 public:
 myVertex(const TopoDS_Vertex vertex):
  vertex(vertex) , point(BRep_Tool::Pnt(vertex)) {};

  const TopoDS_Vertex& getTopoDS_Vertex() const{ return vertex; }

  const int getDim() const { return 0; }
  void paramsOnGeo(const Geo& geom, const double* geo_params, 
		   double* eval_params);

  const arma::mat& D1(const double* u){}
  //void param2xyz(const double *params, double* u) const{}
  const gp_Pnt Value(const double* par) const{ return point; }

};
