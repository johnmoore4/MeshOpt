#pragma once
#include "Mel.h"
#include "ActiveMEl.h"

class ActiveMEl;

class CompEl: public Element{
protected:
  ActiveMEl el;
 

 CompEl(ActiveMEl elT): el(elT){}
public:
 

    arma::vec computeGaussPoints(){};
    virtual const arma::vec computeGaussNormals() const = 0;
    virtual const arma::mat computeJacobianPoint(const arma::vec3[]) const = 0;
    const arma::vec computeDetJacobianGP() const;
    double computeVolume() const;
    ActiveMEl& getActiveElement() { return el; }

};

class CompEl1D: public CompEl{

 public:
  CompEl1D(ActiveMEl el);
  const arma::vec computeGaussNormals() const;
  const arma::vec computeDetJacobianGP() const;
  const arma::mat computeJacobianPoint(const arma::vec3 grad[]) const;

 private:
  arma::rowvec3 pseudoInverse;
};


class CompEl2D: public CompEl{
 private:
  const arma::vec3 lin_normal;
 arma::mat::fixed<3,2> pseudoInverse;

 const arma::mat::fixed<3,2> computePseudoInverse() const;
 public:
 CompEl2D(ActiveMEl elT): 
 CompEl(elT){
   pseudoInverse = computePseudoInverse();
 }

  const arma::vec computeGaussNormals() const {}
  //const arma::vec computeDetJacobianGP() const;
  
  // Maybe put this in a more basic class
  const arma::mat computeJacobianPoint(const arma::vec3[]) const;

  const arma::mat::fixed<3,2>& getPseudoInverse() const{ return pseudoInverse; }
  

};

class CompEl3D: public CompEl{

 public:
 CompEl3D(ActiveMEl elT): CompEl(elT) {}
  const arma::vec computeGaussNormals() const {}
  const arma::vec computeDetJacobianGP() const;
  const arma::mat computeJacobianPoint(const arma::vec3[]) const;
};
