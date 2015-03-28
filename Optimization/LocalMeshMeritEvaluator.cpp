#include "LocalMeshMeritEvaluator.h"
#include "Mel.h"
#include "Mnode.h"

//#include <Teuchos_Array.hpp>
//#include <Teuchos_ArrayView.hpp>

using namespace std;

int LocalMeshMeritEvaluator::NumDOFs(){
  return activenodes.size()*3;
}

void LocalMeshMeritEvaluator::GetCurrentState(double* state){
  node_map& nodes = optel_manager.getNodes();
  int cnt=0;
  for(auto nd = activenodes.begin(); nd != activenodes.end(); ++nd, ++cnt){
    double* temp = nodes.at(nd->first)->xyzptr();
    for(int i=0; i<3; i++) state[cnt*3+i] = temp[i];
  }
}

void LocalMeshMeritEvaluator::SetCurrentState(const double* state){
  node_map& nodes = optel_manager.getNodes();
  int cnt=0;
  for(auto nd = activenodes.begin(); nd != activenodes.end(); ++nd, ++cnt){
    nodes[nd->first]->setXYZ(state+cnt*3);
  }

}

double LocalMeshMeritEvaluator::EvaluateGradient(double* gradient){
  
  using std::cout;
  using std::endl;
  merits.resize(elements.size());

  for(int i=0; i<3*activenodes.size(); i++) gradient[i] = 0.0;

  double merit = 0.0;
  arma::mat gradMerit;
  arma::wall_clock timer;
  timer.tic();
  int cnt=0;
  for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
    std::unique_ptr<OptEl> optel = 
      optel_manager.CreateOptEl(*el,idealElements.at(*el));
    double temp;
    double mer = optel->computeGradMerit(gradMerit,merits[cnt],temp,factor);
    const gind* cn = (*el)->getCornerNodes();
    const int ncn = (*el)->numCornerNodes();
    for(int i=0; i<ncn; i++){
      auto it = activenodes.find(cn[i]);
      if(it != activenodes.end()){
	for(int j=0; j<3; j++){
	  gradient[it->second*3+j]+= gradMerit(j,i);
	}
      }
    }
    merit+= mer;
  }

 
  return merit;
  
}
/*
double LocalMeshMeritEvaluator::
EvaluateHessian(Teuchos::RCP<sparse_mat_type>& HESS, 
			 Teuchos::RCP<multivector_type>& GRAD){
 
  using std::cout;
  using std::endl;

  cout << "in evalute Hessian" << endl;
  bool was_fill_complete = false;
  if(HESS->isFillComplete()){
    HESS->resumeFill();
    was_fill_complete = true;
  }


  HESS->setAllToScalar(0.0);
  GRAD->putScalar(0.0);

  merits.resize(elements.size());
  cout << "Element size: " << elements.size() << endl;

 const int ndof = NumDOFs();

 //for(int i=0; i<ndof; i++) gradient[i] = 0.0;
 //for(int i=0; i<pow(ndof,2); i++) Hessian[i] = 0.0;
 

  double merit = 0.0;
  arma::mat gradMerit, HessianMerit;
  arma::wall_clock timer;
  timer.tic();
  int cnt=0;
  for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
    std::unique_ptr<OptEl> optel = 
      optel_manager.CreateOptEl(*el,idealElements.at(*el));
    //double mer = optel->computeGradMerit(gradMerit,merits[cnt],factor);

    double mer = optel->computeHessianMeritFD(gradMerit,HessianMerit,
					      merits[cnt],factor);

  

    const gind* cn = (*el)->getCornerNodes();
    const int ncn = (*el)->numCornerNodes();

    bool is_active[ncn];
    gind globind[ncn];
    int nact = 0;
    for(int i=0; i<ncn; i++){
      is_active[i] = false;
      auto it = activenodes.find(cn[i]);
      if(it != activenodes.end()){
	is_active[i] = true;
	globind[i] = it->second;
	nact++;
      }
    }
    //cout << "H0" << endl;
    Teuchos::Array<GO>   colindex (3*nact);
    for(int i=0, acnt=0; i<ncn; i++){
      if(is_active[i]){
	for(int j=0; j<3; j++){
	  colindex[acnt*3+j] = globind[i]*3+j;
	}
	acnt++;
      }
     
    }
    Teuchos::Array<ST> colvals(3*nact);

    for(int i=0; i<ncn; i++){
      if(is_active[i]){
	for(int j=0; j<3; j++){
	  for(int k=0, acnt=0; k<ncn; k++){
	    if(is_active[k]){
	      for(int l=0; l<3; l++){
		colvals[acnt*3+l] = HessianMerit(i*3+j,k*3+l);
	      }
	      acnt++;
	    }
	  }
	  GRAD->sumIntoGlobalValue(globind[i]*3+j,0,gradMerit(j,i));
	  //Teuchos::ArrayView<ST> colvals(HessianMerit.colptr(i*3+j),3*ncn);
	  if(was_fill_complete){
	    HESS->sumIntoGlobalValues(globind[i]*3+j,colindex,colvals);
	  }
	  else{
	    HESS->insertGlobalValues(globind[i]*3+j,colindex,colvals);
	  }
	}
      }
    } 


    merit+= mer;
  }
  HESS->fillComplete();
 
  return merit;
  
}
*/

double LocalMeshMeritEvaluator::
EvaluateHessianDense(double* gradient, double* Hessian){

  using std::cout;
  using std::endl;



  merits.resize(elements.size());
  cout << "Element size: " << elements.size() << endl;

 const int ndof = NumDOFs();

 for(int i=0; i<ndof; i++) gradient[i] = 0.0;
 for(int i=0; i<pow(ndof,2); i++) Hessian[i] = 0.0;
 

  double merit = 0.0;
  arma::mat gradMerit, HessianMerit;
  arma::wall_clock timer;
  timer.tic();
  int cnt=0;
  for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
    std::unique_ptr<OptEl> optel = 
      optel_manager.CreateOptEl(*el,idealElements.at(*el));
  

    double mer = optel->computeHessianMeritFD(gradMerit,HessianMerit,
					      merits[cnt],factor);

  

    const gind* cn = (*el)->getCornerNodes();
    const int ncn = (*el)->numCornerNodes();

    
    for(int i=0; i<ncn; i++){
      auto it = activenodes.find(cn[i]);
      if(it != activenodes.end()){
	for(int j=0; j<3; j++){
	  gradient[it->second*3+j]+= gradMerit(j,i);
	  for(int k=0; k<ncn; k++){
	    auto it2 = activenodes.find(cn[k]);
	    if(it2 != activenodes.end()){
	      for(int l=0; l<3; l++){
		Hessian[(it->second*3+j) + (it2->second*3+l)*ndof]+=
		  HessianMerit(i*3+j,k*3+l);
		//Hessian[(it->second*3+j)*ndof + (it2->second*3+l)]+=
		//  HessianMerit(i*3+j,k*3+l);
	      }
	    }
	  }
	}
      }
    }
    

    merit+= mer;
  }
 
 
  return merit;


}
