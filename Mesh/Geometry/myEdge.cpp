#include "myVertex.h"
#include "myEdge.h"
#include "myFace.h"
#include <assert.h>


using namespace std;
using namespace arma;
double myEdge::edgeParamsFromPoint(const double *pt, const double guess) const{
  
  gp_Pnt point(pt[0],pt[1],pt[2]);
  gp_Pnt proj;

  double pr;
  double ret = shape_analysis_curve.NextProject(guess,_geom_curve,
						point,1.0e-10,proj,
						pr,umin,umax);

  //if(ret > 1.0e-8) cout << "projection falied!" << endl;

  return pr;

}

double myEdge::paramsOnEdge(const myVertex* vertex) const {
  auto vert = vertex->getTopoDS_Vertex();
  return  BRep_Tool::Parameter(vertex->getTopoDS_Vertex(),edge);
}

/*
const arma::vec3 myEdge::param2xyz(double s) const{

  vec3 xyz;
  gp_Pnt pnt = curve.Value(s);

  xyz[0] = pnt.X();
  xyz[1] = pnt.Y();
  xyz[2] = pnt.Z();

  return xyz;
}

void myEdge::param2xyz(const double *params, double* u) const{
  gp_Pnt pnt = curve.Value(*params);
 

  //gp_Pnt pnt = _curve.Value(s);
  //for(int i=0; i<3; i++) u[i] = 
}
*/

const myFace* myEdge::findSharedFace(const myEdge* edge) const {
  
 
}

void myEdge::D1(const double* u, arma::rowvec3& deriv){
 
  gp_Pnt P;
  gp_Vec V;
 
  curve.D1(u[0],P,V);
 
  deriv.at(0) = V.X();
  deriv.at(1) = V.Y();
  deriv.at(2) = V.Z();


}
const arma::mat& myEdge::D1(const double* u){
  D1(u,D1_mat);
  return D1_mat;
}


void myEdge::paramsOnGeo(const Geo& geom, const double* geo_params, 
		 double* eval_params){

  if(geom.getDim() != 0){
    assert("Trying to re-paramaterize something"
	   "other than a vertex on a line!");
  }

  
  const double u = BRep_Tool::
    Parameter(static_cast<const myVertex&>(geom).getTopoDS_Vertex(),edge);
  eval_params[0] = u;
  
}

bool myEdge::areParamsValid(const double* params){
  if(params[0] >= umin && params[0] <= umax) return true;
  else return false;
}
