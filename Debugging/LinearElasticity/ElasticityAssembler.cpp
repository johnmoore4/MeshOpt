#include "LinearElasticityAssembler.h"
#include "MeshContainer.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"
#include "ActiveMEl.h"

#include "FEFieldOperator.h"


void ElasticityAssembler::Initialize(){
  node_map& nodes = mesh.getNodesNC();
  const int ndof = nodes.size()*mesh.MeshDimension();
  U.resize(ndof);
  for(int i=0; i<ndof; i++) U[i] = 0.0;
  
}

void ElasticityAssembler::updateR(){

  Assemble(false);


}

void ElasticityAssembler::updateU(double step){
  double* du = system.X();
  
  node_map& nodes = mesh.getNodesNC();

  for(int i=0; i<U.size(); i++){
    U[i]+= du[i];
  }

  const int N = nodemap.size();
  int cnt=0;
  for(auto nd = nodemap.begin(); nd != nodemap.end(); ++nd){
    MNode* node = nodes.at(nd->first).get();
    arma::vec3 temp = node->getXYZ();
    for(int i=0; i<2; i++, cnt++){
      temp[i]+= step*U[cnt];
    }

    node->setXYZ(temp.memptr());
  }
  
}

void ElasticityAssembler::Assemble(bool eval_Jacobian){

  NodeIndexFactory index_factory;
  element_set& elements = mesh.getElementsNC();
  node_map& nodes = mesh.getNodesNC();

  const int ndof = nodemap.size();

  for(auto el = elements.begin(); el != elements.end(); ++el){
    ActiveMEl active_el(el->get(),index_factory,nodes,NULL);
    const ShapeFunctionMatrices* sfh = 
      sf_factory.getShapeFunction((*el)->getElementType(),
				  (*el)->getOrder(),0);

    
    arma::mat u;

    using namespace FEFieldOperator;
    arma::mat** Klin = Laplacian(u);

    //arma::mat** Klin = FE.Div(u) + FE.Laplacian(u);

    //arma::mat** Klin = FE.GradDiv(u) + FE.Laplacian(u);
    
    //arma::mat** Knonlin = FE.Div(F(u));

    //std::pair<arma::vec*,arma::mat**> KFlin = FE.GradDiv(u) + FE.Laplacian(u);

    if(nonlin_type == "linear"){
      
    }
    else if(nonlin_type == "UL"){

    }
    else{
      throw std::runtime_error("Unsupported Elasticity model: " + nonlin_type);
    }

  }



}
