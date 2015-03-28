#pragma once
#include <map>
#include <armadillo>
#include "TrilinosTypedefs.h"
#include "NodeIndexer.h"

class MeshContainer;
class ShapeFunctionMatricesFactory;

class NonlinearElasticitySolver{
 private:
  MeshContainer& mesh;
  ShapeFunctionMatricesFactory& sf_factory;
  double lambda, mu;

  std::map<int,int> nodemap;
  std::map<int,arma::vec2> bcdef;
  std::map<int,arma::vec2> x0;


  Teuchos::RCP<sparse_mat_type> KGLOB;
  Teuchos::RCP<multivector_type> F, DU; 

 protected:

 public:
  NonlinearElasticitySolver(MeshContainer& mesh_container, 
			    ShapeFunctionMatricesFactory& sf_factory):
  mesh(mesh_container), sf_factory(sf_factory){}
  
  int SetUpProblem(const double h, const double lambda_l, const double mu_l);
  int Solve();
  int Assemble();

};
