#include "Mnode.h"
#include "Mel.h"
#include "myVertex.h"
#include "myEdge.h"
#include "myFace.h"


using namespace arma;
using namespace std;


//************** MNode Definitions ********************//
double MNode::computeSubmeshMerit(double &maxmerit){
  /*
  map<int,MEl*>::iterator emap;
  double merit=0.0;
  maxmerit = 0.0;
  for(emap=element_map.begin(); emap!=element_map.end(); ++emap){
    double elmerit = emap->second->computeMeritFunction();
    merit+= elmerit;
    maxmerit = max(maxmerit,elmerit);
  }
  merit/= element_map.size();
  
  return merit;
  */
}

double MNode::getSubmeshMerit(double &maxmerit){
  /*
  map<int,MEl*>::iterator emap;
  double merit=0.0;
  maxmerit = 0.0;
  for(emap=element_map.begin(); emap!=element_map.end(); ++emap){
    double elmerit = emap->second->meritFunction();
    merit+= elmerit;
    maxmerit = max(maxmerit,elmerit);
  }
  merit/= element_map.size();

 

  return merit;
  */
}


void MNode::saveParams(){
  //_uvOld[0] = _uv[0];
  //_uvOld[1] = _uv[1];

}

int MNode::writeTec(ofstream &out){
  for(int i=0; i<3; i++) out << xyz[i] << " ";
  out << endl;
}


void MNode::revertNode(){
  //_uv[0] = _uvOld[0];
  //_uv[1] = _uvOld[1];
  //updateParams(arma::vec(_uvOld,_dim));
}

void MNode::addElement(MEl *element){
  typedef std::pair<int,MEl*> value_type;
  //value_type val(element->getEl(), element);
  //element_map.insert(val);

}

void MNode::eraseElement(MEl *element){
  
  //element_map.erase(element->el());


}

void MNode::optimizeSubmesh(const double distortion_threshold){
  const int nit = 1;
  const int max_tries = 50;
  const double merit_threshold = 0.5*pow(distortion_threshold-1.0,2.0);

  /*

  double maxmerit0, maxmerit0D;
  double merit0 = getSubmeshMerit(maxmerit0);
  //double merit0 = computeSubmeshMerit(maxmerit0);

  //  cout << merit0 << " " << merit0D << " " << maxmerit0 << " " << 
  //  maxmerit0D << endl;

  //if(abs(merit0-merit0D) > 1.0e-8) cout << merit0 << " " << merit0D << endl;


  double merit = merit0;
  double alpha = 1.0;
  //double alpha = step_size;
  vec dparam;
  if(maxmerit0 > merit_threshold){
    int ngood = 0;
    int it = 0;
    while(ngood < nit && it < max_tries){
      //cout << "in here" << endl;
      saveParams();
      vec params_old = getParams();
      map<int,MEl*>::iterator emap;
      vec grad_param(_dim,fill::zeros);
      // mat hess_param(_dim,_dim,fill::zeros);
      //vec gradM(_dim);
      //mat HessM(_dim,_dim);
      for(emap=element_map.begin(); emap!=element_map.end(); ++emap){
	grad_param+= emap->second->getGradMerit(_nd);
	//emap->second->computeGradMeritNode(_nd,gradM,HessM);
	//gradM = 
	//grad_param+= gradM;
      }
      dparam = grad_param;
      if(norm(dparam,2.0) > 1.0) dparam/=norm(dparam,2.0);

      //vec dparam = inv(hess_param)*grad_param;
      vec paramnew = getParams() - alpha*dparam;
      
      updateParams(paramnew);
      double merit_old = merit;
      double maxmerit;
      merit = computeSubmeshMerit(maxmerit);
      //cout << merit << endl;
      
      if(merit > merit_old){
	alpha = 0.5*alpha;
	updateParams(params_old);
	merit = merit_old;
	//cout << "reducing step size to: " << alpha << endl;
      }
      else{
	alpha = min(2.0*alpha,1.0);
	ngood++;
      }
      it++;
    }
    //cout << dparam << endl;
    //cout << "alpha: " << alpha << endl;
    //if(it == max_tries) cout << "did not take " << nit << " good steps!" << endl;

  }
  //cout << "alpha: " << alpha << endl;



  if(merit > merit0){
    cout << "Merit is bigger than merit0!" << endl;
    //cout << merit/merit0 << " " << alpha << endl;
  }
  //if(alpha < 1.0e-2) cout << "alpha: " << alpha << endl;
  */

}

//**************** MVertNode Definitions ********************//

/*
arma::vec MVertNode::paramsOnEdge(const myEdge* edge) const{
  return edge->paramsOnEdge(_vertex);
}

arma::vec MVertNode::paramsOnFace(const myFace* face) const {
  return face->paramsOnFace(_vertex->getTopoDS_Vertex());
}
*/


//**************** MEdgeNode Definitions ********************//

/*
arma::vec MEdgeNode::paramsOnFace(const myFace* face) const {
  return face->paramsOnFace(_edge->getTopoDS_Edge(),_u);
}
*/

/*
vec3 MEdgeNode::param2xyz(double s){

  vec3 xyz;
  gp_Pnt pnt = _edge->getBRepAdaptor_Curve().Value(s);

  xyz[0] = pnt.X();
  xyz[1] = pnt.Y();
  xyz[2] = pnt.Z();

  return xyz;
}
*/

/*
void MEdgeNode::updateParams(arma::vec params){

  double u = params(0);
  _u = u;
  //gp_Pnt pnt = _curve->Value(u);
  gp_Pnt pnt = _edge->getBRepAdaptor_Curve().Value(u);
  
  _xyz[0] = pnt.X();
  _xyz[1] = pnt.Y();
  _xyz[2] = pnt.Z();

}


arma::mat MEdgeNode::D1() const{

  mat::fixed<3,1> d1;

  gp_Pnt P;
  gp_Vec V;
  // _curve->D1(_uv[0],P,V);
  _edge->getBRepAdaptor_Curve().D1(_u,P,V);
 
  d1(0,0) = V.X();
  d1(1,0) = V.Y();
  d1(2,0) = V.Z();

  return d1;
}
*/

//**************** MFaceNode Definitions ********************//

/*
vec3 MFaceNode::param2xyz(double u, double v){

  vec3 xyz;
 
  gp_Pnt pnt = _face->getBRepAdaptor_Surface().Value(u,v);
  xyz[0] = pnt.X();
  xyz[1] = pnt.Y();
  xyz[2] = pnt.Z();

  return xyz;
}
*/


/*
void MFaceNode::updateParams(arma::vec params){

  double u = params(0);
  double v = params(1);
  _uv[0] = u;
  _uv[1] = v;
  gp_Pnt pnt = _face->getBRepAdaptor_Surface().Value(u,v);
 
  _xyz[0] = pnt.X();
  _xyz[1] = pnt.Y();
  _xyz[2] = pnt.Z();
}


arma::mat MFaceNode::D1() const{

  mat::fixed<3,2> d1;

  gp_Pnt P;
  gp_Vec D1U, D1V;

  // _surface->D1(_uv[0],_uv[1],P,D1U,D1V);
  _face->getBRepAdaptor_Surface().D1(_uv[0],_uv[1],P,D1U,D1V);

  d1(0,0) = D1U.X();
  d1(1,0) = D1U.Y();
  d1(2,0) = D1U.Z();

  d1(0,1) = D1V.X();
  d1(1,1) = D1V.Y();
  d1(2,1) = D1V.Z();

  return d1;
}

*/
