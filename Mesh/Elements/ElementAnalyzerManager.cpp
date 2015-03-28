#include "ElementAnalyzerManager.h"
#include "ElementAnalyzer.h"
#include "MeshContainer.h"
#include "ShapeFunctionMatrices.h"
#include "NodeIndexer.h"
#include "Mel.h"
#include "ActiveMEl.h"
#include "CompEl.h"

ElementAnalyzerManager::
ElementAnalyzerManager(MeshContainer& mesh, 
		       ShapeFunctionMatricesFactory& sf_factory,
		       NodeIndexFactory& index_factory):
  mesh(mesh), sf_factory(sf_factory), index_factory(index_factory){
  
}

ElementAnalyzer ElementAnalyzerManager::getElementAnalyzer(const MEl* el){
  int el_type = el->getElementType();
  int order = el->getOrder();
  int quad_degree = 3*order;


  const ShapeFunctionMatrices* sf = 
    sf_factory.getShapeFunction(el_type,order,0,quad_degree);
  
  const ShapeFunctionMatrices* sf_lin = 
    sf_factory.getShapeFunction(el_type,1,0,quad_degree);

  node_map& nodes = mesh.getNodesNC();
  
  std::shared_ptr<CompEl> compel;
  if(el->getDim() == 1){
    compel = std::make_shared<CompEl1D>(ActiveMEl(el,index_factory,nodes,sf));
  }
  if(el->getDim() == 2){
    compel = std::make_shared<CompEl2D>(ActiveMEl(el,index_factory,nodes,sf));
  }
  else if(el->getDim() == 3){
    compel = std::make_shared<CompEl3D>(ActiveMEl(el,index_factory,nodes,sf));
  }
  
  return ElementAnalyzer(compel,sf_lin);

}
