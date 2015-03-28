#pragma once

class SolverParams{


 public:
  

};

class IterativeParams: public SolverParams{

 public:
  int MaxIters;
  double tol;
  

};

class Ifpack2Params: public IterativeParams{

 public:
  std::string method;
  std::string preconditioner;
  bool AdditaveSchwarz;
  int AS_overlap;
  int RILU_k;
  

};
