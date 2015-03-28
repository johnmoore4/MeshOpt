#include "LinearElasticityEvaluator.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"
#include "MeshContainer.h"
#include "ActiveMEl.h"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

#include <MatrixMarket_Tpetra.hpp>
/*
int LinearElasticityEvaluator::Assemble(){

  using arma::diagmat;

  const double mu = 1.0;
  const double lambda = 1.0;
 

  NodeIndexFactory index_factory;
  element_set& elements = mesh.getElementsNC();
  node_map& nodes = mesh.getNodesNC();

  const int ndof = nodemap.size();
  //const int ndof = nodes.size();

  arma::mat K[2][2];

  cout << "before element loop" << endl;
  arma::mat nodesel;
  for(auto el = elements.begin(); el != elements.end(); ++el){
    ActiveMEl active_el(el->get(),index_factory,nodes,NULL);
    const ShapeFunctionMatrices* sfh = 
      sf_factory.getShapeFunction((*el)->getElementType(),
				  (*el)->getOrder(),0);

    const arma::mat& sf = sfh->getQuadratureSF();
    const arma::cube& sfderiv = sfh->getQuadratureSFDeriv();
    const arma::vec& gw = sfh->getQuadratureWeights();


    nodesel = active_el.getNodesMatrix();

    //arma::mat nodes2d = nodesel.cols(0,1);
    

    arma::vec xxi = sfderiv.slice(0).t()*nodesel.unsafe_col(0);
    arma::vec xet = sfderiv.slice(1).t()*nodesel.unsafe_col(0);
    arma::vec yxi = sfderiv.slice(0).t()*nodesel.unsafe_col(1);
    arma::vec yet = sfderiv.slice(1).t()*nodesel.unsafe_col(1);
    arma::vec jacGW = (xxi%yet - xet%yxi)%gw;
    

    // Derivatives are pre-multiplied with detJ*gw
    arma::mat shapx = sfderiv.slice(0)*diagmat(yet) - 
      sfderiv.slice(1)*diagmat(yxi);

    arma::mat shapy = -sfderiv.slice(0)*diagmat(xet) + 
      sfderiv.slice(1)*diagmat(xxi);




    const double mu2lam = 2.0*mu + lambda;
    K[0][0] = 
      shapx*diagmat(mu2lam*gw)*shapx.t() + 
      shapy*diagmat(mu*gw)*shapy.t();

 
    K[0][1] = 
      shapx*diagmat(lambda*gw)*shapy.t() +
      shapy*diagmat(mu*gw)*shapx.t();

    K[1][0] = 
      shapy*lambda*diagmat(lambda*gw)*shapx.t() +
      shapx*diagmat(mu*gw)*shapx.t();
    
    K[1][1] = 
      shapy*diagmat(mu2lam*gw)*shapy.t() +
      shapx*diagmat(mu*gw)*shapx.t();


    const int ncn = (*el)->numCornerNodes();
    const gind* cn = (*el)->getCornerNodes();
    
    cout << "h0" << endl;
    arma::uvec globind(ncn), locind(ncn);
    int cnt=0;
    for(int n=0; n<ncn; n++){
      if(nodemap.find(cn[n]) != nodemap.end()){
	//globind(cnt) = cn[n];
	globind(cnt) = nodemap[cn[n]];
	locind(cnt) = n;
	cnt++;
      }
      else{

	arma::vec2 du = bcdef[cn[n]];
	for(int i=0; i<cnt; i++){
	  for(int j=0; j<2; j++){
	    for(int k=0; k<2; k++){
	      
	      RHS->sumIntoGlobalValue(globind[i]+ndof*j,0,
	      			      -K[j][k](locind(i),n)*du(j));
	    }
	  }
	}
      }
    }
    globind.resize(cnt);

    Teuchos::Array<GO> colind(cnt);
    Teuchos::Array<ST> colvals(cnt);

    cout << "h1" << endl;

    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
	for(int k=0; k<cnt; k++){
	  for(int l=0; l<cnt; l++){
	    colvals[l] = K[i][j](locind(k),locind(l));
	    colind[l] = globind(l)+j*ndof;
	  }
	  KGLOB->insertGlobalValues(globind(k)+i*ndof,colind,colvals);
	}
      }
    }

  }
  KGLOB->fillComplete();

  Tpetra::MatrixMarket::Writer<sparse_mat_type>::writeSparseFile("testfile",KGLOB);

  cout << "At end of assembler" << endl;
  return 1;
}
*/


int LinearElasticityEvaluator::Assemble(const double lambda, const double mu){
  /*
  using namespace std;

  using arma::diagmat;

  bool was_fill_complete = false;
  if(KGLOB->isFillComplete()){
    KGLOB->resumeFill();
    was_fill_complete = true;
  }


  KGLOB->setAllToScalar(0.0);
  RHS->putScalar(0.0);


  //const int ndof = nodemap.size();

  NodeIndexFactory index_factory;
  element_set& elements = mesh.getElementsNC();
  node_map& nodes = mesh.getNodesNC();
  const int ndof = nodes.size();

  arma::mat K[2][2];
  arma::vec F[2];

  arma::mat nodesel;
  for(auto el = elements.begin(); el != elements.end(); ++el){
    ActiveMEl active_el(el->get(),index_factory,nodes,NULL);
    
    const int ncn = (*el)->numCornerNodes();
    const gind* cn = (*el)->getCornerNodes();
    arma::mat f(ncn,2);
    for(int i=0; i<ncn; i++){
      arma::vec2& temp = forces[cn[i]];
      f(i,0) = temp(0);
      f(i,1) = temp(1);
    }

    const ShapeFunctionMatrices* sfh = 
      sf_factory.getShapeFunction((*el)->getElementType(),
				  (*el)->getOrder(),0);

    const arma::mat& sf = sfh->getQuadratureSF();
    const arma::cube& sfderiv = sfh->getQuadratureSFDeriv();
    const arma::vec& gw = sfh->getQuadratureWeights();


    nodesel = active_el.getNodesMatrix();

    //arma::mat nodes2d = nodesel.cols(0,1);
    
    arma::mat nodeselT = nodesel.t();

    arma::vec xxi = sfderiv.slice(0).t()*nodeselT.unsafe_col(0);
    arma::vec xet = sfderiv.slice(1).t()*nodeselT.unsafe_col(0);
    arma::vec yxi = sfderiv.slice(0).t()*nodeselT.unsafe_col(1);
    arma::vec yet = sfderiv.slice(1).t()*nodeselT.unsafe_col(1);
    //arma::vec jacGW = (xxi%yet - xet%yxi)%gw;
    arma::vec detjac = (xxi%yet - xet%yxi);
    for(int i=0; i<detjac.size(); i++){
      if(std::abs(detjac(i)) < 1.0e-10){
	cout << detjac << endl;
	cout << nodesel << endl;
	cout << "det is 0!" << endl;
      }
    }

    // Derivatives are pre-multiplied with detJ*gw
    arma::mat shapx = (sfderiv.slice(0)*diagmat(yet) - 
		       sfderiv.slice(1)*diagmat(yxi));

    arma::mat shapy = (-sfderiv.slice(0)*diagmat(xet) + 
		       sfderiv.slice(1)*diagmat(xxi));

    arma::mat shapxx = shapx*diagmat(gw/detjac)*shapx.t();
    arma::mat shapxy = shapx*diagmat(gw/detjac)*shapy.t();
    arma::mat shapyy = shapy*diagmat(gw/detjac)*shapy.t();

    //arma::mat test = shapy*diagmat(gw)*shapx.t();

 
    //const double area = arma::sum(detjac%gw);
    

    const double mu2lam = 2.0*mu + lambda;

    K[0][0] = (mu2lam*shapxx + mu*shapyy);

    K[0][1] = (lambda*shapxy + mu*shapxy.t());
    //K[0][1] = (lambda*shapxy.t() + mu*shapxy);

    K[1][0] = (lambda*shapxy.t() + mu*shapxy);
    
    K[1][1] = mu2lam*shapyy + mu*shapxx;


    arma::mat fgp = sf.t()*f;

    F[0] = sf*(detjac%gw%fgp.col(0));
    F[1] = sf*(detjac%gw%fgp.col(1));
    //const int ncn = (*el)->numCornerNodes();
    //const gind* cn = (*el)->getCornerNodes();
    
    arma::uvec globind(ncn), locind(ncn);
    int cnt=0;
    for(int n=0; n<ncn; n++){

      if(nodemap.find(cn[n]) != nodemap.end()){
	for(int i=0; i<2; i++){
	  //RHS->sumIntoGlobalValue(cn[n]+i*ndof,0,F[i](n));
	}
	globind(cnt) = cn[n];
	//globind(cnt) = nodemap[cn[n]];
	locind(cnt) = n;
	cnt++;
      }
      else{
	Teuchos::Array<GO> ti(1);
	Teuchos::Array<ST> tv(1);
	for(int i=0; i<2; i++){
	  ti[0] = cn[n] + i*ndof;
	  tv[0] = 1.0;
	  arma::vec2 du = bcdef[cn[n]];
	  KGLOB->insertGlobalValues(cn[n]+i*ndof,ti,tv);
	  RHS->sumIntoGlobalValue(cn[n]+i*ndof,0,du(i));
	}
      }
    }
    globind.resize(cnt);

    Teuchos::Array<GO> colind(ncn);
    Teuchos::Array<ST> colvals(ncn);


    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
	for(int k=0; k<cnt; k++){
	  for(int l=0; l<ncn; l++){
	    colvals[l] = K[i][j](locind(k),l);
	    colind[l] = cn[l]+j*ndof;

	    //colvals[l] = K[i][j](locind(k),locind(l));
	    //colind[l] = globind(l)+j*ndof;
	  }
	  if(was_fill_complete){
	    KGLOB->sumIntoGlobalValues(globind(k)+i*ndof,colind,colvals);
	  }
	  else{
	    KGLOB->insertGlobalValues(globind(k)+i*ndof,colind,colvals);
	  }
	}
      }
    }

  }
  KGLOB->fillComplete();

  //Tpetra::MatrixMarket::Writer<sparse_mat_type>::writeSparseFile("testfile",KGLOB);

  //cout << "At end of assembler" << endl;
  return 1;
  */
}

