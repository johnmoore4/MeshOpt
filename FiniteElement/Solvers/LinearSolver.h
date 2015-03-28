#pragma once

#include "armadillo"

class SolverParams;

class LinearSolver{
 private:
  SolverParams& params;

 public:
 LinearSolver(SolverParams& params_t): params(params_t) {}
  virtual arma::vec Solve() = 0;
  virtual void AddToSystem(arma::mat& A, arma::vec& b, arma::uvec& ind) = 0;
  virtual void Reset() = 0;
  virtual const SolverParams& getSolverParams() const{ return params; } 


};

class IterativeSolver: public LinearSolver {

 public:
 IterativeSolver(SolverParams& params_t): LinearSolver(params_t) {}
  virtual void AddToSystem(arma::mat& A, arma::vec& b, arma::uvec& ind) = 0;
};


