#pragma once

#include "Geo.h"
#include "OCC_standard_includes.h"
#include <memory>
#include <vector>
#include "armadillo"

class myVertex;
class myEdge;
class myFace;
//typedef std::shared_ptr<myFace> face_ptr;

class myFace: public Geo{
  
 private:
  Handle(Geom_Surface) _geom_surface;
  const TopoDS_Face face;
  const BRepAdaptor_Surface surface;


  std::unique_ptr<ShapeAnalysis_Surface> _sas;
  double umin, umax, vmin, vmax;
  arma::mat::fixed<2,3> D1_mat;

 public:

 myFace(const TopoDS_Face face):  
  _geom_surface(BRep_Tool::Surface(face)),
    face(face), 
    surface(face){
      _sas = std::unique_ptr<ShapeAnalysis_Surface>
	(new ShapeAnalysis_Surface(_geom_surface));
      ShapeAnalysis::GetFaceUVBounds(face, umin, umax, vmin, vmax);
      const double du = umax - umin;
      const double dv = vmax - vmin;
      //umin -= fabs(du) / 100.0;
      //vmin -= fabs(dv) / 100.0;
      //umax += fabs(du) / 100.0;
      //vmax += fabs(dv) / 100.0;
  }

  arma::vec2 paramsOnFace(const myVertex* vertex) const;

  arma::vec2 paramsOnFace(const myEdge* edge, double u) const;

  arma::vec2 faceParamsFromPoint(const double *pt, const double *guess) const;

  /*
  arma::vec3 params2xyz(double u, double v) const;
  void param2xyz(const double *params, double* u) const;
  */

  const TopoDS_Face& getTopoDS_Face() const{ return face;}

  const BRepAdaptor_Surface& getBRepAdaptor_Surface() const{ return surface;}

  void D1(const double* uv, arma::mat& deriv);
  const arma::mat& D1(const double* uv);

  const int getDim() const { return 2; }
  void paramsOnGeo(const Geo& geom, const double* geo_params, 
		   double* eval_params);

  const gp_Pnt Value(const double* par) const{ 
    return surface.Value(par[0],par[1]); 
  }

  bool areParamsValid(const double* params);
};
