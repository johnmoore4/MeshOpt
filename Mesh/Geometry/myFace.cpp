#include "myFace.h"
#include "myVertex.h"
#include "myEdge.h"
#include <assert.h>

using namespace std;
using namespace arma;

arma::vec2 myFace::paramsOnFace(const myVertex* vertex) const {
  gp_Pnt2d uvgp = BRep_Tool::Parameters(vertex->getTopoDS_Vertex(),face);
  arma::vec2 uv;
  uvgp.Coord(uv(0),uv(1));
  return uv;
}

arma::vec2 myFace::paramsOnFace(const myEdge* edge, double u) const {
  Handle(Geom2d_Curve) curve_handle;
  double f,l;
  curve_handle = BRep_Tool::CurveOnSurface(edge->getTopoDS_Edge(),face,f,l); 
  gp_Pnt2d point = curve_handle->Value(u);
  arma::vec2 uv;
  point.Coord(uv(0),uv(1));
  return uv;
}


arma::vec2 myFace::faceParamsFromPoint(const double *pt, const double *guess) 
  const{

  gp_Pnt pnt(pt[0],pt[1],pt[2]);
  gp_Pnt2d prev(guess[0],guess[1]);

  arma::vec2 pr;
  gp_Pnt2d proj = _sas->NextValueOfUV(prev,pnt,1.0e-10);
  pr[0] = proj.X();
  pr[1] = proj.Y();

  return pr;
}

/*
arma::vec3 myFace::params2xyz(double u, double v) const {

  vec3 xyz;
 
  gp_Pnt pnt = _surface.Value(u,v);
  xyz[0] = pnt.X();
  xyz[1] = pnt.Y();
  xyz[2] = pnt.Z();

  return xyz;
}

void myFace::param2xyz(const double *params, double* u) const{


}
*/


void myFace::D1(const double* uv, arma::mat& deriv){

  gp_Pnt P;
  gp_Vec D1U, D1V;

 
  surface.D1(uv[0],uv[1],P,D1U,D1V);

  deriv(0,0) = D1U.X();
  deriv(0,1) = D1U.Y();
  deriv(0,2) = D1U.Z();

  deriv(1,0) = D1V.X();
  deriv(1,1) = D1V.Y();
  deriv(1,2) = D1V.Z();

}

const arma::mat& myFace::D1(const double* uv){
  D1(uv,D1_mat);
  return D1_mat;
}

void myFace::paramsOnGeo(const Geo& geom, const double* geo_params, 
		 double* eval_params){

  const int dim = geom.getDim();

  if(dim < 2){
    assert("Trying to re-paramaterize something"
	   "other than a vertex or a line on a face!");
  }

  gp_Pnt2d uvgp;

  if(dim == 0){
    const myVertex& vertex = static_cast<const myVertex&>(geom);
     uvgp = BRep_Tool::Parameters(vertex.getTopoDS_Vertex(),face);
  }
  else if(dim == 1){
    const myEdge& edge = static_cast<const myEdge&>(geom);
    Handle(Geom2d_Curve) curve_handle;
    double f,l;
    curve_handle = BRep_Tool::CurveOnSurface(edge.getTopoDS_Edge(),face,f,l); 
    uvgp = curve_handle->Value(geo_params[0]);
  }

  uvgp.Coord(eval_params[0],eval_params[1]);

}
bool myFace::areParamsValid(const double* params){
  if(params[0] >= umin && params[0] <= umax &&
     params[1] >= vmin && params[1] <= vmax) return true;
  else return false;
}
