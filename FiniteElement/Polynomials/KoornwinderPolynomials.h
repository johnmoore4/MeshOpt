#pragma once
#include <armadillo>

class KoornwinderPolynomial{

 public:
  virtual void EvaluatePolynomial(const int p, const arma::mat& x, 
			  arma::mat &f, arma::cube &df) = 0;

};

class Koornwinder2D: public KoornwinderPolynomial{
 public:
  void EvaluatePolynomial(const int p, const arma::mat& x, 
			  arma::mat &f, arma::cube &df);
};

class Koornwinder3D: public KoornwinderPolynomial{
 public:
  void EvaluatePolynomial(const int p, const arma::mat& x, 
			  arma::mat &f, arma::cube &df);
};

