#pragma once
#include "Assembler.h"
#include "TrilinosLinearSystem.h"

#include <map>
#include <armadillo>

class MeshContainer;
class ShapeFunctionMatricesFactory;

class ElasticityAssembler: public Assembler{
 private:
  MeshContainer& mesh;
  ShapeFunctionMatricesFactory& sf_factory;

  TrilinosLinearSystem system;

  std::map<int,int>& nodemap;
  std::map<int,arma::vec2>& bcdef;
  std::string nonlin_type;

  //std::string nonlin_type = "linear";
  double mu = 1.0;
  double lambda = 1.0;
  std::vector<double> U;

 protected:
 
  

 public:
 ElasticityAssembler(MeshContainer& mesh,
		     ShapeFunctionMatricesFactory& sf_factory,
		     std::map<int,int>& nodemap,
		     std::map<int,arma::vec2>& bcdef,
		     std::string nonlin_type="linear"):
  mesh(mesh), sf_factory(sf_factory), nodemap(nodemap), bcdef(bcdef),
    nonlin_type(nonlin_type){}

  void Initialize();
  void Assemble(bool eval_Jacobian);
  void updateR();
  void updateU(double step);

  LinearSystem& getLinearSystem(){ return system; }

};

/*
class LinearElasticityAssembler: public ElasticityAssembler{
 private:

 protected:

 public:
  void Assemble(){}


};
*/
