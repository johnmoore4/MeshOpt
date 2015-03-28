#pragma once
#include "LinearSolver.h"
#include "SolverParams.h"

class Ifpack2Solver: public IterativeSolver{
  

 public:
 Ifpack2Solver(SolverParams& params_t): IterativeSolver(params_t){}
  void AddToSystem(arma::mat& A, arma::vec& b, arma::uvec& ind);
};
