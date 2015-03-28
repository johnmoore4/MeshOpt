#pragma once
#include <map>
#include <armadillo>
#include "TrilinosTypedefs.h"
#include "NodeIndexer.h"
#include "System.h"
#include "GlobalDefines.h"
#include "TrilinosLinearSystem.h"

class MeshContainer;
class ShapeFunctionMatricesFactory;
class MEl;
//class LinearElasticityAssembler;

// Define
// Grad(u), Div(u)
// Div(q)
// u, q


// Linear Elasticity
// UResidual = (lambda+mu)*GradDiv(u) + mu*Laplacian(u);
// Or UResidual = (lambda+mu)*Grad(Div(u)) + mu*Div(Grad(u));

// Div u...
// Utype divU = shapx*u(:,1) + shapy*u(:,2) + shapz*u(:,3);
// Grad(divU) : [divU*shapx.t(), divU*shapy.t(), divU*shapz.t()];


/// Grad u....
// Grad(u) : ;


/// Navier-Stokes
// QResidual = q - Grad(u);
// UResidual = Div(F(q,u)) + DotNormal(F(q,u)) + tau*(u-uhat);
// UhatResidual = DotNormal(F(q,u)) + tau*(u-uhat) + BC;

// Assembler.assemble(UResidual)
// Assembler.assemble(QResidual)


class CurrentState{
 private:
  std::vector<double> Q, U, UHAT, SOURCE;
 protected:

 public:

};

class ElementBlock{
 private:

 public:

 protected:

};

class LocalState{
 private:
  arma::mat q, u, uhat, source;
 protected:

 public:
  
};

class Assembler{
 private:
  // DistributedMesh& mesh;
  //node_map& nodes;
  //element_set& elements;
  LinearSystem& system;
  CurrentState& state;

  std::string Galerkin_type;
  int N_unknowns;

  //LocalSatate local_state;
 protected:

 public:
  int InitializeLinearSystem(int N_unknowns=1, std::string Galerkin_type="CG"){}
  int AddToSystem(ElementBlock& block){}
  const LocalState& GetLocalState(const MEl* el){};
};

class ElasticitySystem: public System{
 private:
  // DistributedMesh& mesh
  MeshContainer& mesh;
  ShapeFunctionMatricesFactory& sf_factory;

 
  std::map<gind,arma::vec>& Dirichlet, Forces;

  TrilinosLinearSystem system;
  CurrentState state;

 protected:

 public:
   ElasticitySystem(MeshContainer& mesh,
			     ShapeFunctionMatricesFactory& sf_factory,
			     std::map<gind,arma::vec>& Dirichlet,
			     std::map<gind,arma::vec>& Forces):
  mesh(mesh), sf_factory(sf_factory), Dirichlet(Dirichlet), Forces(Forces){}

  LinearSystem& getLinearSystem(){ return system; }

  void Assemble(bool eval_Jacobian){};
  void updateR(){};
  void updateU(double step){};




  //int SetUpProblem(const double h);
  //int Assemble(const double lambda, const double mu);
  int Solve();
  //int UpdateMesh();

};
