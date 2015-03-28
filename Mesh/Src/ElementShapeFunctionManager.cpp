#include "ElementShapeFunctionManager.h"

#include "ShapeFunctionMatrices.h"
#include "Mel.h"
#include "CEl.h"

int ElementShapeFunctionManager::
getJacobianDegree(int dim, int order, bool is_curved){
  int total_degree = 0;
  if(is_curved){
    total_degree+= dim*(order-1);
  }
  return total_degree;
}

ElementShapeFunctionManager::ElementShapeFunctionManager(){
  sf_factory = std::make_shared<ShapeFunctionMatricesFactory>();
}

const ShapeFunctionMatrices* ElementShapeFunctionManager::
getSFMatrices(int type, int order, int interp_type, int degree){
  return sf_factory->getShapeFunction(type,order,interp_type,degree);
}

const ShapeFunctionMatrices* ElementShapeFunctionManager::
getSFofIntegrationDegree(const MEl& el, int degree_int){
  int total_degree = 
    getJacobianDegree(el.getDim(),el.getOrder(),el.isCurved()) + degree_int;

  return getSFMatrices(el.getElementType(),el.getOrder(),0,total_degree);
}
  
const ShapeFunctionMatrices* ElementShapeFunctionManager::
getSFofIntegrationDegree(const CEl& el, int degree_int){
  //std::cout << el.is_curved << " ";
  int total_degree = 
    getJacobianDegree(el.dim,el.porder,el.is_curved) + degree_int;
  
  return getSFMatrices(el.el_type,el.porder,0,total_degree);

}
