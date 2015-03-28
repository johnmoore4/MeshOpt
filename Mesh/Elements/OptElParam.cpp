#include "OptElParam.h"
#include "OptEl.h"
#include "CompEl.h"
#include "GeometryContainer.h"

OptElParam::OptElParam(CompEl& compel, 
		       GeometryContainer *geometry, 
		       const arma::mat* ideal):
  compel(compel), geometry(geometry), ideal(ideal){
}

double OptElParam::computeGradMeritParam(arma::mat& gradMerit, 
					 double& distortion,
					 double factor, 
					 double minDetS) const{
  

  arma::mat grad_temp;
  double detS;
  double mer = optel->computeGradMerit(grad_temp,distortion,detS,factor,minDetS);

  ActiveMEl& active = compel.getActiveElement();
  const MEl* el = active.getMeshElement();
  const gind* nind = el->getNodes();


  node_map& nodes = active.getGlobalNodes();
  const int nn = active.Ndof();
  

  gradMerit.resize(3*nn,1);
  gradMerit.zeros();
  
  for(int i=0, dofcnt=0; i<nn; i++){
    const MNode* nd = nodes[nind[i]].get();
    const int dim = nd->getType();
    const short int entity = nd->getGeoEntity();
 
    if(nd->getType() == 1){
      std::vector<myEdge>& edges = geometry->getEdgesNC();
      myEdge& edge = edges[entity];
      
      const double* u = nd->getParametricCoords();
      arma::rowvec3 deriv;
      edge.D1(u,deriv);
      const double uderiv = arma::as_scalar(deriv*grad_temp.unsafe_col(i));
      gradMerit(dofcnt) = uderiv;
    }
    else if(nd->getType() == 2){
      if(1){
	std::vector<myFace>& faces = geometry->getFacesNC();
	myFace& face = faces[entity];
	const double* uv = nd->getParametricCoords();
	arma::mat::fixed<2,3> deriv;
	face.D1(uv,deriv);
	arma::vec2 uvderiv = deriv*grad_temp.unsafe_col(i);
	for(int j=0; j<2; j++){
	  gradMerit(dofcnt+j) = uvderiv[j];
	}
	//cout << "face: " << entity << endl;
      }
    }
    else if(nd->getType() == 3){
      for(int j=0; j<3; j++){
	gradMerit(dofcnt+j) = grad_temp.at(j,i);
      }
    }

    dofcnt+= dim;
    
    //if(i == ncn-1) gradMerit.resize(dofcnt);
  }
  
  //gradMerit.resize(dofcnt);
  return mer;

}

double OptElParam::computeDistortion() const{
  return optel->computeDistortion();
}
double OptElParam::computeGradMerit(arma::mat& gradMerit, double& distortion,
				    double& DetS, double factor,
				    double minDetS) const{
  return optel->computeGradMerit(gradMerit,distortion,DetS,factor,minDetS);
}
 

double OptElParam::computeMerit() const{
  return optel->computeMerit();
}

double OptElParam::computeMinDetS() const{
  return optel->computeMinDetS();
}
