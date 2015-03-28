#include "MeshAnalyzer.h"

#include "MeshContainer.h"
//#include "OptElManager.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"
#include "ElementAnalyzerManager.h"
#include "ElementAnalyzer.h"

MeshAnalyzer::MeshAnalyzer(MeshContainer& mesh): mesh(mesh){

}

const element_set& getElements(MeshContainer& mesh, int dim){
  int mesh_dim = mesh.MeshDimension();
  if(dim == mesh_dim) return mesh.getElements();
  else if (dim == mesh_dim-1) return mesh.getSubElements();
  else return mesh.getSubSubElements();
}

int MeshAnalyzer::ComputeMeshQualities(int dim){
  std::cout << "h0" << std::endl;
  NodeIndexFactory index_factory;
  ShapeFunctionMatricesFactory sf_factory;
  
  ElementAnalyzerManager 
    element_analyzer_manager(mesh,sf_factory,index_factory);


  int elcnt = 0;
  const element_set& elements = getElements(mesh,dim);

  std::cout << "h1" << std::endl;
  //const element_set& elements = mesh.getElements();
  qualities[dim-1].resize(elements.size());
  for(auto it = elements.begin(); it != elements.end(); ++it, ++elcnt){
    const MEl* el = it->get();
    ElementAnalyzer analyzer = element_analyzer_manager.getElementAnalyzer(el);
    qualities[dim-1][elcnt] = 1.0/analyzer.computeDistortion();
  }


  if(elements.size() > 0){
    auto minmax = std::minmax_element(qualities[dim-1].begin(),
				      qualities[dim-1].end());

    std::cout << dim << "D min/max quality: " << *minmax.first << " " << 
      *minmax.second << std::endl;
  }
}

int MeshAnalyzer::Analyze(){
  int mesh_dim = mesh.MeshDimension();
  
  for(int dim = mesh_dim; dim > 0; dim--){
    ComputeMeshQualities(dim);
  }
  
  return 0;
}

int MeshAnalyzer::SetMeshCurvedElements(){
  for(int dim = 1; dim <= mesh.MeshDimension(); dim++){
    const element_set& elements = getElements(mesh,dim);
    int elcnt = 0;
    for(auto el = elements.begin(); el != elements.end(); ++el, ++elcnt){
      if(1-qualities[dim-1][elcnt] > 1.0e-14){
	(*el)->setIsCurved(false);
      }
    }
  }

}
