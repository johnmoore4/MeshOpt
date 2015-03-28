#pragma once
#include "armadillo"


void polyval(const arma::vec &c, const arma::vec &x, arma::vec& val);
void polyder(const arma::vec &c, arma::vec& dr);
arma::vec jacobi(const int n, const double a, const double b);
void Jacobi1D(const int p, const arma::mat& x, arma::mat &f, 
	      arma::cube& df);

