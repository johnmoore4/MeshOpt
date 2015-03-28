#include "LinearElasticityEvaluator.h"
#include "NewtonSolver.h"

//#include "computeSolutionBelos.h"


int LinearElasticityEvaluator::Solve(){

  //ElasticityAssembler assembler(mesh,sf_factory,nodemap,bcdef);

  NewtonSolver Newton(*this);

  //NewtonSolver Newton(assembler);
  //Newton.Solve(1.0e-14);

  /*
  arma::wall_clock timer;
  std::cout << "In solve" << std::endl;

  timer.tic();
  // Solve the system
  computeSolutionBelos(KGLOB,RHS,DU,1,0,1000,1000,1.0e-14,false,"RILUK");
  std::cout << "Solve time: " << timer.toc() << std::endl;

  */ 

  

  return 1;
}

