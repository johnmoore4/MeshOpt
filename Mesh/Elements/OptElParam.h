#pragma once
#include <memory>

#include <armadillo>

class OptEl;
class CompEl;
class GeometryContainer;

class OptElParam{
 public:
  OptElParam(CompEl& compel, 
	     GeometryContainer *geometry, 
	     const arma::mat* ideal=NULL);

  double computeGradMeritParam(arma::mat& gradMerit, 
			       double& distortion,
			       double factor= 1.0, 
			       double minDetS=0.0) const;


  double computeDistortion() const;
  double computeGradMerit(arma::mat&, double& distortion,
			  double& DetS, double factor=1.0,
			  double minDetS=0.0) const;
 
  double computeMerit() const;
  double computeMinDetS() const;

  //const double computeDebug() const;
  //const arma::mat debugGradMerit() const;
  //const arma::mat debugGradMeritFD() const;

  
 private:
  std::shared_ptr<OptEl> optel;

  const arma::mat* ideal;
  CompEl& compel;
  GeometryContainer* geometry;
  
};
