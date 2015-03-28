#pragma once

/*
namespace arma{
  template <class T> class Mat;
  template <class T> class Vec;
}
*/
#include <armadillo>

class Quadrature{
 private:
  
 protected:

 public:
  virtual void ComputeQuadrature
    (const int p, arma::mat& x, arma::vec& w) const = 0;

};

class LineQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;

};

class TriQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;
};

class QuadQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;
};

class TetQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;
};

class HexQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;
};

class PrismQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;
};

class PyramidQuadrature: public Quadrature{

 public:
  void ComputeQuadrature(const int p, arma::mat& x, arma::vec& w) const;
};


class QuadratureFactory{
 private:
  LineQuadrature Line;
  TriQuadrature Tri;
  QuadQuadrature Quad;
  TetQuadrature Tet;
  HexQuadrature Hex;
  PrismQuadrature Prism;
  PyramidQuadrature Pyramid;

 public:
  const Quadrature& getQuadrature(const int type){

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
    }
  }
  
};
