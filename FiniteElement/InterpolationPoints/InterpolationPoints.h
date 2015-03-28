#pragma once
#include <armadillo>

class InterpolationPoints{


 public:
  virtual arma::mat ComputePoints(const int p, const int type) const = 0;

};

class LinePoints: public InterpolationPoints{
  
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class TriPoints: public InterpolationPoints{
  
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class QuadPoints: public InterpolationPoints{
  
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class TetPoints: public InterpolationPoints{
  
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class HexPoints: public InterpolationPoints{
  
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class PrismPoints: public InterpolationPoints{
  
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class PyramidPoints: public InterpolationPoints{
 public:
  arma::mat ComputePoints(const int p, const int type) const;
};

class InterpolationPointFactory{

 private:
  LinePoints Line;
  TriPoints Tri;
  QuadPoints Quad;
  TetPoints Tet;
  HexPoints Hex;
  PrismPoints Prism;
  PyramidPoints Pyramid;

 public:
  const InterpolationPoints& getInterpolationPoints(const int type){
    switch(type){
    case 1:
      return Line;
    case 2:
      return Tri;
    case 3:
      return Quad;
    case 4:
      return Tet;
    case 5:
      return Hex;
    case 6:
      return Prism;
    case 7:
      return Pyramid;
    default:
      throw std::runtime_error("Unknown Element type of "+std::to_string(type));
    }
  }
};
