#pragma once
#include <armadillo>

class PolynomialBasis{


 public:
  virtual void EvalBasis(const int p, const arma::mat& evalpts, arma::mat& f, 
			 arma::cube& df, arma::cube& df2) const = 0;

};

class LinePolynomialBasis: public PolynomialBasis{
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class TriPolynomialBasis: public PolynomialBasis{
 public:
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class QuadPolynomialBasis: public PolynomialBasis{
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class TetPolynomialBasis: public PolynomialBasis{
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class HexPolynomialBasis: public PolynomialBasis{
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class PrismPolynomialBasis: public PolynomialBasis{
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class PyramidPolynomialBasis: public PolynomialBasis{
 public:
  void EvalBasis(const int p, const arma::mat& evalpts, 
		 arma::mat& f, arma::cube& df, arma::cube& df2) const;
};

class PolynomialBasisFactory{
 private:
  LinePolynomialBasis Line;
  TriPolynomialBasis Tri;
  QuadPolynomialBasis Quad;
  TetPolynomialBasis Tet;
  HexPolynomialBasis Hex;
  PrismPolynomialBasis Prism;
  PyramidPolynomialBasis Pyramid;

 public:

  const PolynomialBasis& getPolynomialBasis(const int type){
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
