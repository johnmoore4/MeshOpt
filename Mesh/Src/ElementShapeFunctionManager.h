#pragma once

#include <memory>

class ShapeFunctionMatricesFactory;
class ShapeFunctionMatrices;
class MEl;
class CEl;

class ElementShapeFunctionManager{
 public:
  ElementShapeFunctionManager();

  int getJacobianDegree(int dim, int order, bool is_curved);

  const ShapeFunctionMatrices* getSFMatrices(int type, int order, 
					     int interp_type, int degree);
  
  const ShapeFunctionMatrices* getSFofIntegrationDegree(const MEl& el, 
							int degree_int);
  
  const ShapeFunctionMatrices* getSFofIntegrationDegree(const CEl& el, 
							int degree_int);
 private:
  std::shared_ptr<ShapeFunctionMatricesFactory> sf_factory;

  
};
