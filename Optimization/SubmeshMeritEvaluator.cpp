#include "SubmeshMeritEvaluator.h"
#include "Mel.h"
#include "Mnode.h"
#include "GeometryContainer.h"

void SubmeshMeritEvaluator::GetCurrentState(double* state){
  //double* temp = nd->xyzptr();
  const double* temp = nd->getParametricCoords();

  for(int i=0; i<nd->getType(); i++) state[i] = temp[i];
  //for(int i=0; i<2; i++) temp[i] = state[i];
  
  //state = nd->xyzptr();

}

void SubmeshMeritEvaluator::SetCurrentState(const double* state){
  
  double temp[] = {0.0,0.0,0.0};
  
  for(int i=0; i<nd->getType(); i++) temp[i] = state[i];

  //for(int i=0; i<2; i++) temp[i] = state[i];

  //nd->setXYZ(temp);

  
  arma::vec3 newpos;
  const short int entity = nd->getGeoEntity();
  if(nd->getType() != 3){
    GeometryContainer& geometry = optel_manager.getGeometry();
    if(nd->getType() == 1){
      std::vector<myEdge>& edges = geometry.getEdgesNC();
      myEdge& edge = edges[entity];
      //const double* u = nd->getParametricCoords();
      if(edge.areParamsValid(state)){
	edge.param2xyz(state,newpos.memptr());
	nd->setParametricCoords(temp);
	nd->setXYZ(newpos.memptr());
      }
      else{
	//std::cout << state[0] << std::endl;
	//std::cout << "edge params are not valid" << std::endl;
      }
    }
    else if(nd->getType() == 2){
      std::vector<myFace>& faces = geometry.getFacesNC();
      myFace& face = faces[entity];
      //const double* uv = nd->getParametricCoords();
      if(face.areParamsValid(state)){
	face.param2xyz(state,newpos.memptr());
	nd->setParametricCoords(temp);
	nd->setXYZ(newpos.memptr());
      }
      else{
	//std::cout << "face params are not valid" << std::endl;
      }
    }
    //nd->setXYZ(newpos.memptr());
  }
  else{
    nd->setParametricCoords(temp);
  }
  

}

double SubmeshMeritEvaluator::EvaluateGradient(double* gradient){
  
  using std::cout;
  using std::endl;
  //std::cout << "in eval gradient" << std::endl;

  const node_map& nodes = optel_manager.getNodes();

  //cout << "begenning eval gradient" << endl;
  merits.resize(elements.size());
  double merit = 0.0;
  arma::mat gradMerit;
  //arma::vec3 gradNode = {0.0, 0.0, 0.0};
  arma::vec gradNode(nd->getType(),arma::fill::zeros);

  //arma::vec gradNode(2,arma::fill::zeros);

  /*
  double minDetS = 1.0/0.0;
  for(auto el = elements.begin(); el != elements.end(); ++el){
    double temp=1.0/0.0;
    double dist;
    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el->el,el->ideal);
    optel->computeGradMerit(gradMerit,dist, temp, factor);
    minDetS = std::min(minDetS,temp);
  }
  */
  if(elements.size() == 2){
    //std::cout << "element size is 2" << std::endl;
  }
  int cnt=0;
  for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
 
    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el->el,el->ideal);

    
    double temp, minDetS;
    //double mer = optel->computeGradMerit(gradMerit,merits[cnt], factor,minDetS);

    //double mer = optel->computeMerit();
    //gradMerit = optel->computeGradMeritFD();
    
    //std::cout << "mer: " << mer << std::endl;
    //std::cout << gradMerit << std::endl;

    double mer = optel->computeGradMeritParam(gradMerit,merits[cnt], factor,
    					      minDetS);
    
    //std::cout << mer << std::endl;
    //std::cout << gradMerit << std::endl;
    
    //double mer = optel->computeGradMeritParam(gradMerit,merits[cnt], factor);
  
    //merits[cnt] = (sqrt(mer) + 1.0)/10.0;
    //const gind* cn = el->el->getCornerNodes();
    //const int ncn = el->el->numCornerNodes();

    const gind* elnodes = el->el->getNodes();
    const int nn = el->el->NumNodes();

    bool found = false;
    for(int i=0, dofcnt=0; i<nn; i++){
 
      if(elnodes[i] == nd->getND()){
	//std::cout << nd->xyzvec3().t() << std::endl;
	//for(int j = 0; j < 2; j++){
	for(int j=0; j<nd->getType(); j++){
	  gradNode[j]+= gradMerit[dofcnt+j];
	  //std::cout << gradNode[j] << " ";
	}
	//std::cout << std::endl;
	found = true;
      }
      //dofcnt+= 3;
      dofcnt+= nodes.at(elnodes[i])->getType();
    }
    
    assert(found);

    merit+= mer;
  }

  //cout <<  "submesh merit: " << merit << endl;

  //cout << "endign eval gradient" << endl;
  for(int i=0; i<nd->getType(); i++) gradient[i] = gradNode[i];
  //for(int i = 0; i < 2; i++) gradient[i] = gradNode[i];

  
  return merit;
  
}
