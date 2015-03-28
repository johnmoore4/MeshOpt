#include "ElasticitySystem.h"
#include "MeshContainer.h"
#include "ShapeFunctionMatrices.h"
#include "NewtonSolver.h"


int ElasticitySystem::Solve(){
  NewtonSolver Newton(*this);
  Newton.Solve(1.0e-14);

}
