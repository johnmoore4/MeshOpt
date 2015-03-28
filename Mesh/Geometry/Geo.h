#pragma once
#include <armadillo>
#include "OCC_standard_includes.h"

class Geo{
 protected:
  int bc_tag=0;


 public:
  Geo(){};
  //Geo(Handle(Geom_Geometry)& subgeom): geometry(subgeom){}


  virtual const int getDim() const = 0;
  virtual void paramsOnGeo(const Geo& geom, const double* geo_params, double* 
			  eval_params) = 0;

  virtual const arma::mat& D1(const double* u) = 0;
  virtual const gp_Pnt Value(const double* par) const = 0;


  void param2xyz(const double *params, double* xyz) const;
  

  //virtual const arma::vec3 myEdge::param2xyz(double s) const

  void setBCTag(int tag){ bc_tag = tag; }
  const int getBCTag() const{ return bc_tag; }
 
};
