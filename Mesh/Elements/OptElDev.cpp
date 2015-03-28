#include "OptEl.h"
#include "ShapeFunctionMatrices.h"

void computeHessDetJacobian(arma::mat& HessJ, const arma::mat& J) {

  
  //arma::mat99 HessJ;
  //HessJ.zeros();
  const unsigned int ind[3][2] = {{1,2},{0,2},{0,1}};
  /*
  arma::umat ind;
  ind << 1 << 2 << arma::endr <<
    0 << 2 << arma::endr <<
    0 << 1;
  */

  const int sgn[3] = {1, -1, 1};
  //arma::vec3 sgn = {1.0, -1.0, 1.0};
  arma::mat22 J22;

  const int temp[4] = {3,1,2,0};
  const int s[4] = {1,-1,-1,1};

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<2; k++){
	for(int l=0; l<2; l++){
	  J22.at(k,l) = J.at(ind[i][k],ind[j][l]);
	}
      }
      for(int k=0,cnt=0; k<2; k++){
	for(int l=0; l<2; l++,cnt++){
	  HessJ.at(i+3*j,ind[i][k]+3*ind[j][l]) = sgn[i]*sgn[j]*s[cnt]*
	    J22.at(temp[cnt]);
	}
      }
    }
  }

  //return HessJ;
  
}
/*
const arma::mat computeHessSigma(double detS, double delta, 
				 const arma::mat&  gradDetS,
				 const arma::mat& HessDetS){
  arma::mat temp = gradDetS.t();

  double denom = sqrt(detS*detS + 4*delta*delta);
  arma::mat gradijDetS(HessDetS.n_rows,HessDetS.n_cols);
  for(int i=0; i<HessDetS.n_rows; i++){
    for(int j=0; j<HessDetS.n_cols; j++){
      gradijDetS(i,j) = temp(i)*temp(j);
    }
  }
  return 0.5*(HessDetS*(1.0+detS/denom) + 
	      gradijDetS*(1.0/denom - pow(detS,2)/pow(denom,3)));

}
*/

void computeHessFroNorm(arma::mat& Hess, 
			const arma::mat& DS, 
			const arma::mat& invJI,
			const double fro){

  
  const int dim = DS.n_rows;
  const int dim2 = dim*dim;
  const double fro3Inv = 1.0/(fro*fro*fro);
   
  for(int i=0; i<dim2; i++){
    for(int j=0; j<dim2; j++){
      Hess.at(i,j) =  -DS.at(i)*DS.at(j)*fro3Inv;
      if(i/3 == j/3){
	Hess.at(i,j)+= dot(invJI.unsafe_col(i%3),invJI.unsafe_col(j%3))/fro;
      }
    }
  }
 
 
}
arma::mat computeHessEta(const double fronorm, const double sigma,
			 const arma::mat& gradFroNorm,
			 const arma::mat& GradSigma,
			 const arma::mat& HessSigma,
			 const arma::mat& HessFroNorm){

  
  return arma::mat();
}


void assembleHessian(double HessArea[18][18],
		     arma::mat& gradMerit,
		     arma::mat& gradArea,
		     arma::mat& HessianMerit,
		     arma::mat& HessianArea,
		     const arma::cube& dsf, 
		     const arma::mat& gradJ,
		     const arma::mat& DS,
		     const arma::mat& HessJ,
		     const arma::mat& HessFro,
		     const double gradJ_coeff,
		     const double gradDS_coeff,
		     const double HessHJ_coeff,
		     const double HessF_coeff,
		     const double HessJ2_coeff,
		     const double HessJDS_coeff,
		     const double HessDS2_coeff,
		     const double w, 
		     const int gpi,
		     const double dj){
  //double Hess_DJS = eta_HessJDS;
  // double Hess_DS2 = eta_HessDS2;

  const int dim = dsf.n_slices;
  const int ndof = dsf.n_rows;

  
  for(int j=0, jk1=0; j<ndof; j++){
    for(int k=0, cnt1=0; k<3; k++, jk1++){
      for(int l=0; l<3; l++,cnt1++){
	const int ind1 = k + dim*l;
	gradMerit(k,j)+= (gradJ_coeff*gradJ.at(k,l) + gradDS_coeff*DS.at(k,l))*
	  dsf.at(j,gpi,l)*w;
	gradArea.at(k,j)+= dj*gradJ.at(k,l)*dsf.at(j,gpi,l)*w;
	
	for(int j2 = 0, jk2=0; j2<ndof; j2++){
	  for(int k2=0,cnt2=0; k2<3; k2++, jk2++){
	    for(int l2=0; l2<3; l2++,cnt2++){
	      const int ind2 = k2+dim*l2;
	      HessianMerit.at(jk1,jk2)+= 
		(		 
		 HessHJ_coeff*HessJ.at(ind1,ind2) + 
		 HessF_coeff*HessFro.at(ind1,ind2) +
		 HessJ2_coeff*gradJ.at(k,l)*gradJ.at(k2,l2) +
		 HessJDS_coeff*(DS.at(k,l)*gradJ.at(k2,l2) + 
				DS.at(k2,l2)*gradJ.at(k,l)) +
		 HessDS2_coeff*DS.at(k,l)*DS.at(k2,l2)
		 )*dsf.at(j,gpi,l)*dsf.at(j2,gpi,l2)*w;
	      //HessArea[jk1][jk2]+= 
	      // 	HessJ.at(ind1,ind2)*dsf.at(j,gpi,l)*dsf.at(j2,gpi,l2);

	      HessianArea.at(jk1,jk2)+= 
		HessJ.at(ind1,ind2)*dsf.at(j,gpi,l)*dsf.at(j2,gpi,l2)*w;
	    }
	  }
	}
      }
    }
  }
  
}


const arma::mat OptEl::computeHessianMerit() const{
  ActiveMEl& active = compel.getActiveElement();

  const arma::mat nodes = active.getNodesMatrix();

  const double eps = 2.0*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::vec& gw = 
      (active.getShapeFunctionMatrices())->getQuadratureWeights();

 
  const int dim = dsf.n_slices;
  const int dim2 = dim*dim;

  arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal*dsf.slice(i);
  }
  const int ngp = dsf.n_cols;
  const int ndof = dsf.n_rows;
  //double element_area = 0.0;
  double eta_shape = 0.0;
  double area = 0.0;
  arma::vec3 derivpt[3], derivptI[3];

  /*
  double J_mem[dim2], invJt_mem[dim2], JI_mem[dim2];
  double S_mem[dim2], invJI_mem[dim2], HessJ_mem[dim2*dim2];
  double DS_mem[dim2], HessFro_mem[dim2*dim2];
  */
  /*
  double J_mem[9], invJt_mem[9], JI_mem[9];
  double S_mem[9], invJI_mem[9], HessJ_mem[81];
  double DS_mem[9], HessFro_mem[81];

  double gradMerit_mem[dim*ndof];

  
  arma::mat J(J_mem,dim,dim,false);
  arma::mat JI(JI_mem,dim,dim,false);
  arma::mat invJI(invJI_mem,dim,dim,false);
  arma::mat S(S_mem,dim,dim,false);
  arma::mat DS(DS_mem,dim,dim,false);
  arma::mat invJt(invJt_mem,dim,dim,false);
  arma::mat HessJ(HessJ_mem,dim2,dim2,false);
  arma::mat gradMerit(gradMerit_mem,ndof,3,false);
  arma::mat HessFro(HessFro_mem,dim2,dim2,false);
  */

  
  arma::mat J, JI, invJI, S, invJt, DS;
  
 
  arma::mat HessJ(dim2,dim2,arma::fill::zeros);
  arma::mat HessFro(dim2,dim2);
 
  arma::mat gradMerit(3,ndof,arma::fill::zeros);
  arma::mat gradArea(3,ndof,arma::fill::zeros);
  arma::mat HessianMerit(3*ndof,3*ndof,arma::fill::zeros);
  arma::mat HessianArea(3*ndof,3*ndof,arma::fill::zeros);
  //arma::mat::fixed<18,18> HessianArea;
  double HessArea[18][18];
  HessianArea.zeros();

  gradMerit.zeros();
  HessJ.zeros();

  for(int gp=0; gp<ngp; gp++){

    
    for(int j=0; j<dim; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(gp);
      derivptI[j] = xyzderivI[j].unsafe_col(gp);
    }
    J = compel.computeJacobianPoint(derivpt);
    JI = compel.computeJacobianPoint(derivptI);
 

    invJt = inv(J).t();
    invJI = inv(JI);
    S = invJI*J;
    DS = invJI.t()*S;


    double detJ = det(J);
    double detS = det(S);
    double detInvJI = det(invJI);

    //arma::mat gradDetJ =computeGradDeterminant(J);
    //arma::mat gradDetS = det(invJI)*gradDetJ;

    const double dj = std::abs(det(J));
    // arma::mat gradDetJ = computeGradDeterminant(J);
 
    double delta;
    if(detS < eps) delta = sqrt(eps*(eps-detS));
    else delta = 0.0;
 

    const double stab = sqrt(detS*detS + 4.0*delta*delta);
    const double sigma = 0.5*(detS + stab);

    const double sigma_gradJ = 0.5*detJ*detInvJI*(1.0 + detS/stab);
    const double sigma_HessHJ = 0.5*detInvJI*(1.0 + detS/stab);
    const double sigma_HessJ2 = 0.5*detJ*(1.0/stab - detS*detS/pow(stab,3));
		
    const double fronorm = norm(S,"fro");

    const double d = double(dim);
    const double d3 = double(dim)/3.0;
    const double eta = fronorm*fronorm/(d*pow(sigma,d3));



    const double eta_gradJ = 
      -pow(fronorm,2)/(3.0*pow(sigma,d3+1.0))*sigma_gradJ;
    const double eta_gradDS = 2.0/(d*pow(sigma,d3));    
    const double eta_HessHJ = 
      -pow(fronorm,2)/(3.0*pow(sigma,d3+1.0))*sigma_HessHJ;

    const double eta_HessJ2 = 
      (d3+1.0)*pow(fronorm,2)/(3.0*pow(sigma,d3+2.0))*
      pow(sigma_gradJ,2) - pow(fronorm,2)/(3.0*pow(sigma,d3+1.0))*sigma_HessJ2;

    const double eta_HessJDS = -2.0/3.0*sigma_gradJ/pow(sigma,d3+1.0);
    const double eta_HessF = 2.0*fronorm/(d*pow(sigma,d3));
    const double eta_HessDS2 = 2.0/(d*pow(sigma,d3)*pow(fronorm,2.0));
 
  
    computeHessDetJacobian(HessJ,J);
    computeHessFroNorm(HessFro,DS,invJI,fronorm);

    const double w = gw.at(gp);

    const double grad_J = 2.0*eta*dj*eta_gradJ + dj*eta*eta;
    const double grad_DS = 2.0*eta*dj*eta_gradDS;
    const double Hess_HJ = (2.0*dj*eta*eta_HessHJ + pow(eta,2.0));
    const double Hess_F = 2.0*dj*eta*eta_HessF;
    const double Hess_J2 = (2.0*dj*(eta*eta_HessJ2 + pow(eta_gradJ,2)) + 
			      4.0*eta*detJ*eta_gradJ);
    const double Hess_DJS = 
      (2.0*dj*(eta*eta_HessJDS + eta_gradJ*eta_gradDS) + 
	 2.0*eta*detJ*eta_gradDS);
    const double Hess_DS2 = 2.0*dj*(eta*eta_HessDS2 + pow(eta_gradDS,2));
    
    
    //arma::mat HessDetS = det(invJI)*computeHessDetJacobian(J);
    assembleHessian(HessArea,gradMerit,gradArea,HessianMerit,HessianArea,
		    dsf,invJt,DS,HessJ,HessFro,grad_J,grad_DS,
     		    Hess_HJ,Hess_F,Hess_J2,Hess_DJS,Hess_DS2,w,gp,dj);
    
    
    area+= dj*w;
    eta_shape+= pow(eta,2)*dj*w;
    
    /*
    for(int j=0; j<ndof; j++){
      for(int k=0; k<3; k++){
	for(int l=0; l<3; l++){
	  for(int j2 = 0; j2<ndof; j2++){
	    for(int k2=0; k2<3; k2++){
	      for(int l2=0; l2<3; l2++){
		HessianMerit.at(k+j*3,k2+j2*3)+= 
		  (
		   Hess_HJ*HessJ.at(k+dim*l,k2+dim*l2) + 
		   Hess_F*HessFro.at(k+dim*l,k2+dim*l2) +
		   Hess_J2*invJt.at(k,l)*invJt.at(k2,l2) +
		   Hess_DJS*(DS.at(k,l)*invJt.at(k2,l2) + 
			     DS.at(k2,l2)*invJt.at(k,l)) +
		   Hess_DS2*DS.at(k,l)*DS.at(k2,l2)
		   )*dsf.at(j,gp,l)*dsf.at(j2,gp,l2);

	      }
	    }
	  }
	}
      }
    }
    */

    /*

 
    
    double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));
    arma::mat gradSigma = 0.5*(gradDetS + detS*gradDetS/
     			       sqrt(detS*detS + 4.0*delta*delta)); 

    arma::mat HessSigma = computeHessSigma(detS,delta,gradDetS,HessDetS);
    //accumulator+= computeHessSigma(detS,delta,gradDetS,HessDetS);
    double fronorm = norm(S,"fro");
    arma::mat gradFroNorm = S.t()*invJI/fronorm;
    arma::mat HessFroNorm = computeHessFroNorm(S,invJI);

    double denom = dim*pow(sigma,double(dim)/double(3));
    //arma::mat gradDenom = 
    //  dim*(dim/3)*pow(sigma,double(dim)/double(3)-1)*gradSigma;
    

    double eta = fronorm*fronorm/denom;

    arma::mat HessEta = computeHessEta(fronorm,sigma,gradFroNorm,gradSigma,
				       HessSigma,HessFroNorm);

    //accumulator+= HessFroNorm;

    //arma::mat HessDetJ = computeHessDetJacobian(J);
    */
    /*
    arma::mat HessDetJ(9,9);
    {

      double step = 1.0e-8;
      for(int i=0; i<9; i++){
	arma::mat33 JFD = J;
	JFD(i)+= step;
	arma::mat gradDetJFD = computeGradDeterminant(JFD);
	HessDetJ.col(i) = vectorise(gradDetJFD.t()-gradDetJ.t())/step;
      }
    }
    */
  
  }
 
  arma::mat temp = 
    HessianMerit*area + vectorise(gradMerit)*vectorise(gradArea).t() +
    (vectorise(gradMerit)*vectorise(gradArea).t()).t() + eta_shape*HessianArea;

  


  return temp;
  return HessianMerit;
  
}

const double OptEl::computeHessianMeritFD(arma::mat& gradMerit, 
					  arma::mat& HessianMerit,
					  double& dist, double factor) const{
  
  /*
  double DetS = 0.0;
  ActiveMEl& active = compel.getActiveElement();
  node_map& nodes = active.getGlobalNodes();



  //const arma::mat gradMerit = debugGradMerit();
  
  const MEl* el = active.getMeshElement();
  const int ncn = el->numCornerNodes();
  const gind* cn = el->getCornerNodes();

  double step = 1.0e-8;
  const int ndof = ncn*3;
  HessianMerit.resize(ndof,ndof);
  gradMerit.resize(3,ncn);

  const double mer = computeGradMerit(gradMerit,dist,DetS,factor);
  for(int i=0; i<ncn; i++){
    arma::vec3 xyz = nodes[cn[i]]->getXYZ();

    for(int j=0; j<3; j++){
 
      arma::vec3 xyzfd = xyz;
      xyzfd(j)+= step;
      nodes[cn[i]]->setXYZ(xyzfd.memptr());
      arma::mat gradMeritFDplus;
      double temp;

      computeGradMerit(gradMeritFDplus,temp,DetS,factor);

      xyzfd(j)-= 2.0*step;
      nodes[cn[i]]->setXYZ(xyzfd.memptr());
      arma::mat gradMeritFDminus;
      computeGradMerit(gradMeritFDminus,temp,DetS,factor);
      for(int k=0; k<ncn; k++){
	for(int l=0; l<3; l++){
	  HessianMerit(j+i*3,l+k*3) = 
	    (gradMeritFDplus(l,k)-gradMeritFDminus(l,k))/step/2.0;
	}
      }
      nodes[cn[i]]->setXYZ(xyz.memptr());
    }
  }

  return mer;
  */
}

/// ----------------------------------------------

const arma::mat OptEl::debugGradMerit() const{
  /*
  ActiveMEl& active = compel.getActiveElement();
  const arma::mat nodes = active.getNodesMatrix();

  const double eps = 2.0*std::numeric_limits<double>::epsilon();

  const arma::cube& dsf = 
    (active.getShapeFunctionMatrices())->getQuadratureSFDeriv();

  const arma::vec& gw = 
      (active.getShapeFunctionMatrices())->getQuadratureWeights();


 
  const int dim = dsf.n_slices;


  arma::mat xyzderiv[3];
  for(int i=0; i<dim; i++){
    xyzderiv[i] = nodes*dsf.slice(i);
  }
 
  arma::mat xyzderivI[3];
  for(int i=0; i<dim; i++){
    xyzderivI[i] = ideal*dsf.slice(i);
  }
  const int ngp = dsf.n_cols;
  const int ndof = dsf.n_rows;
  double element_area = 0.0;
  double eta_shap = 0.0;
  arma::mat J, JI, invJI, S;
  arma::mat iden = arma::eye(dim,dim);
  arma::mat gradsf(3,ndof);
  arma::mat gradDebug(3,ndof,arma::fill::zeros);
  arma::mat gradEtaShapeAccum(3,ndof,arma::fill::zeros);
  arma::mat gradAreaAccum(3,ndof,arma::fill::zeros);

 
  for(int i=0; i<ngp; i++){
    arma::vec3 derivpt[3], derivptI[3];
    for(int j=0; j<dim; j++){
      derivpt[j] = xyzderiv[j].unsafe_col(i);
      derivptI[j] = xyzderivI[j].unsafe_col(i);
    }
    J = compel.computeJacobianPoint(derivpt);
    JI = compel.computeJacobianPoint(derivptI);
 
    
    double dj = det(J);
    //double dj = std::abs(det(J));
    arma::mat gradDetJ =computeGradDeterminant(J);
  
    

 

    for(int j=0; j<dim; j++){
      for(int k=0; k<ndof; k++) gradsf(j,k) = dsf(k,i,j);
    }


    
    invJI = inv(JI);
    S = invJI*J;
    double detS = det(S);
    arma::mat gradDetS = det(invJI)*gradDetJ;

    //gradDebug+= gradDetS*gradsf;
    
    double delta;
    if(detS < eps) delta = sqrt(eps*(eps-detS));
    else delta = 0.0;
    
    double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));
    arma::mat gradSigma = 0.5*(gradDetS + detS*gradDetS/
			       sqrt(detS*detS + 4.0*delta*delta)); 
    
    // gradDebug+= gradSigma*gradsf;
    double fronorm = norm(S,"fro");
    arma::mat gradFroNorm = invJI.t()*S/fronorm;

    //gradDebug+= gradFroNorm*gradsf;
    if(i == -2){
      arma::mat99 Hess;
      Hess.zeros();
      double step = 1.0e-8;
      for(int i=0; i<9; i++){
	arma::mat33 JFD = J;
	JFD(i)+= step;
	arma::mat SFD = invJI*JFD;
	double froFD = norm(SFD,"fro");
	arma::mat gradFroFD = SFD.t()*invJI/froFD;
	Hess.col(i) = vectorise(gradFroFD.t()-gradFroNorm.t())/step;
      }
      cout << "Hess FD: " << endl;
      cout << Hess << endl;
    }
 


    double denom = dim*pow(sigma,double(dim)/double(3));
    arma::mat gradDenom = dim*(dim/3)*pow(sigma,double(dim)/double(3)-1)*gradSigma;
    

    double eta = fronorm*fronorm/denom;

    arma::mat gradEta = 2.0*fronorm*gradFroNorm/denom -
    fronorm*fronorm/pow(denom,2)*gradDenom;
    //arma::mat gradEta = 2.0*fronorm*gradFroNorm;
    //arma::mat gradEta = 2.0*fronorm*gradFroNorm/sigma - 
    //  pow(fronorm,2)/pow(sigma,2)*gradSigma;
    //arma::mat gradEta = 2.0*fronorm*gradFroNorm/sigma;
    //arma::mat gradEta = -pow(fronorm,2)/pow(sigma,2)*gradSigma;
    //arma::mat gradEta = gradFroNorm*sigma;
    //arma::mat gradEta = gradFroNorm*sigma;
    //arma::mat gradEta = arma::ones<arma::mat>(3,3)*sigma;
    //

    //arma::mat gradEtaShape = 2.0*eta*gradEta*dj + 0.0*pow(eta,2)*gradDetJ;
    //arma::mat gradEtaShape = 2.0*eta*gradEta*dj;
    //arma::mat gradEtaShape = gradDetJ;
    gradEtaShapeAccum+= 
      (2.0*eta*gradEta*dj + pow(eta,2)*gradDetJ)*gradsf*gw(i);
    gradAreaAccum+= gradDetJ*gradsf*gw(i);
    element_area+= dj*gw(i);
    eta_shap+= eta*eta*dj*gw(i);
  

    //gradDebug+= gradDetJ*gradsf*gw(i);

 
  }
  gradDebug = gradEtaShapeAccum*element_area + eta_shap*gradAreaAccum;



  return gradDebug;

*/
}
