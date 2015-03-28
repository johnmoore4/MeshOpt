#pragma once
#include "Mel.h"
#include "ActiveMEl.h"
#include "CompEl.h"

class myFace;
class myEdge;
class GeometryContainer;

template<class T>
inline T mySign(const T s) { return ((s >= 0.0) ? 1.0 : -1.0); }

template<class T>
inline T Kronecker(const T i, const T j){ return (i == j) ? 1 : 0; } 

class OptEl: public Element{
 private:

protected:
  const arma::mat& ideal;
  const ShapeFunctionMatrices* sf_ideal;
  
  //ActiveMEl& active;
  //const CompEl& compel;
  CompEl& compel;
  GeometryContainer* geometry;

  //std::unique_ptr<CompEl>& compel;
  virtual void computeJacobianPoint(const double *deriv[], 
				    arma::mat& J) const = 0;
  virtual void computeGradDetS(const arma::mat& J, 
			       arma::mat& grad,
			       const double coeff) const = 0;

  virtual void computeGradFrobNorm(const arma::mat& J,
				   const arma::mat& invJI,
				   arma::mat& grad) const = 0;

public:
 OptEl(CompEl& compT,
       const arma::mat& idealT, 
       const ShapeFunctionMatrices* sf_ideal,
       GeometryContainer *geometry_t): 
  compel(compT), ideal(idealT), sf_ideal(sf_ideal),geometry(geometry_t){}
  double computeDistortion() const;
  const double computeGradMerit(arma::mat&, double& distortion,
				double& DetS, double factor=1.0,
				double minDetS=0.0) const;
 
  const double computeGradMeritParam(arma::mat&, double& distortion,
				     double factor= 1.0, 
				     double minDetS=0.0) const;

  const arma::mat computeGradMeritFD() const;
  const double computeMerit() const;
  const double computeMinDetS() const;

  const double computeDebug() const;
  const arma::mat debugGradMerit() const;
  const arma::mat debugGradMeritFD() const;



  const arma::mat computeHessianMerit() const;
  const double computeHessianMeritFD(arma::mat& gradMerit,
					arma::mat& HessianMerit,
					double& dist, double factor) const;


};


class OptEl2D: public OptEl{
 private:
  //arma::vec3 normal;
  std::unique_ptr<CompEl> compel2D;
 
 protected:
  void computeJacobianPoint(const double* deriv[], arma::mat& J) const;

  void computeGradDetS(const arma::mat& J, 
		       arma::mat& grad,
		       const double coeff = 1.0) const;

  void computeGradFrobNorm(const arma::mat& J,
			   const arma::mat& invJI,
			   arma::mat& grad) const;

  //const arma::mat computeGradDeterminant(const arma::mat& J) const;
  //const arma::mat computeGradFrobNorm(const arma::mat& J) const;
  CompEl2D compel2d;
 public:
  OptEl2D();
  OptEl2D(const MEl* el, 
	  NodeIndexFactory& index_factory, 
	  node_map& nodes, 
	  const ShapeFunctionMatrices* sf, 
	  const ShapeFunctionMatrices* sf_ideal,
	  const arma::mat& ideal,
	  GeometryContainer *geometry_t=NULL):

  OptEl(compel2d,ideal,sf_ideal,geometry_t), 
    compel2d(ActiveMEl(el,index_factory,nodes,sf)){}


  /*
 OptEl2D(ActiveMEl& elT, const arma::mat& idealT): 
  OptEl(elT,compel2D,idealT){
    
    
    if(!idealT.is_empty()){
      normal = arma::cross(idealT.col(1)-idealT.col(0),
			   idealT.col(2)-idealT.col(0));
      normal/=arma::norm(normal,2.0);
    }
    compel2D =
      std::unique_ptr<CompEl2D>(new CompEl2D(elT,normal,!idealT.is_empty()));
  }
  */


};


class OptEl3D: public OptEl{
  // const CompEl3D& compel;
  //std::unique_ptr<CompEl> compel3D;
 private:

 protected:
  void computeJacobianPoint(const double *deriv[], arma::mat& J) const;

  void computeGradDetS(const arma::mat& J, 
		       arma::mat& grad,
		       const double coeff = 1.0) const;
 
  void computeGradFrobNorm(const arma::mat& J,
			   const arma::mat& invJI,
			   arma::mat& grad) const;

  //const arma::mat computeGradDeterminant(const arma::mat& J) const;
  //const arma::mat computeGradFrobNorm(const arma::mat& J) const;
  //ActiveMEl active;
  CompEl3D compel3d;
 public:
  OptEl3D();
  OptEl3D(const MEl* el, 
	  NodeIndexFactory& index_factory, 
	  node_map& nodes, 
	  const ShapeFunctionMatrices* sf, 
	  const ShapeFunctionMatrices* sf_ideal,
	  const arma::mat& ideal,
	  GeometryContainer *geometry_t=NULL):
  OptEl(compel3d,ideal,sf_ideal,geometry_t), 
    compel3d(ActiveMEl(el,index_factory,nodes,sf)){}

  // OptEl3D(ActiveMEl& elT, const arma::mat& idealT): 
  //OptEl(elT,compel3D,idealT){
  //compel = std::unique_ptr<CompEl3D>(new CompEl3D(el));
  

};

