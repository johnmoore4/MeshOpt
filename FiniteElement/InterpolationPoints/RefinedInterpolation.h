#pragma once

#include <map>
//#include <vector>
#include <array>
#include <armadillo>
#include <memory>

class InterpolationPointFactory;
class PolynomialBasisFactory;
class PolynomialBasis;
class InterpolationPoints;

class RefinedInterpolation{
 public:
  arma::mat& getInterpolationMatrix(){ return interpolationMatrix; }
  arma::umat& getConnectivity(){ return connectivity; }
  int Refine(int Nref);

 protected:
 RefinedInterpolation(int Nref, int Nnode, int Nadd, int Nsub, int p,
		      const PolynomialBasis& basis, 
		      const InterpolationPoints& points): 
  nref(Nref), nnode(Nnode), nadd(Nadd), nsub(Nsub), p(p), basis(basis), 
    points(points) {}
  //virtual int Refine() = 0;

  std::map<std::array<double,3>, int> refinedPoints;
  std::vector<std::array<double,3> > refinedPointsVec;
  std::vector<std::vector<int> > recursionConnectivity;

  arma::mat interpolationMatrix;
  arma::umat connectivity;
  arma::umat edgesList;
  const int nref;
  const int nnode;
  const int nadd;
  const int nsub;
  const int p;
  const PolynomialBasis& basis;
  const InterpolationPoints& points;
  
  int point_counter = 0;

 private:

};

class TriRefinedInterpolation: public RefinedInterpolation{
 public:
  TriRefinedInterpolation(int Nref, int p, 
			  const PolynomialBasis& basis,
			  const InterpolationPoints& points);

 protected:
 
 private:

};

class QuadRefinedInterpolation: public RefinedInterpolation{
 public:
  QuadRefinedInterpolation(int Nref, int p, 
			  const PolynomialBasis& basis,
			  const InterpolationPoints& points);

 protected:
 
 private:

};

class TetRefinedInterpolation: public RefinedInterpolation{
 public:
  TetRefinedInterpolation(int Nref, int p, 
			  const PolynomialBasis& basis,
			  const InterpolationPoints& points);

 protected:
 
 private:

};

class PrismRefinedInterpolation: public RefinedInterpolation{
 public:
  PrismRefinedInterpolation(int Nref, int p, 
			  const PolynomialBasis& basis,
			  const InterpolationPoints& points);

 protected:
 
 private:

};


class RefinedInterpolationFactory{
public:
  RefinedInterpolationFactory();
  RefinedInterpolation& getRefinedInterpolation(int Nref, int p,int eltype);
protected:

private:
  std::vector<std::shared_ptr<RefinedInterpolation> > interpolations;

  std::shared_ptr<InterpolationPointFactory> interpolation_factory;
  std::shared_ptr<PolynomialBasisFactory> basis_factory;

};
