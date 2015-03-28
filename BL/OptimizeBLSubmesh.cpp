#include "BoundaryLayer.h"
#include "MeshContainer.h"
#include "ActiveMEl.h"
#include "OptEl.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"
#include "MeritEvaluator.h"

#include <lbfgs.h>
#include <map>
#include <vector>
#include "Mel.h"

typedef std::map<gind,std::vector<MEl*> > node2element_type;

class SubmeshOptimizer{
private:
  node2element_type Node2ElementMap;
  //MeshContainer& mesh;
  ShapeFunctionMatricesFactory sf_factory;
  NodeIndexFactory index_factory;
public:
  //SubmeshOptimizer(MeshContainer& mesh_t): mesh(mesh_t){}
  void ComputeNode2ElementMap(int Nlayers, double threshold);
  double OptimizeSubmesh(gind nd, std::vector<MEl*> elements);
};

class SubmeshMeritEvaluator: public MeritEvaluator{
private:
  std::vector<MEl*>& elements;
protected:

public:
  //SubmeshMeritEvaluator(gind nd, std::vector<MEl*>& elements_t, MeshContainer&
  //			mesh): elements(elements_t), MeritEvaluator(mesh){}
};

void SubmeshOptimizer::ComputeNode2ElementMap(int Nlayers, double threshold){
  /*
  ideal_map& idealElements = mesh.getIdealElementsNC();
  const element_set& elements = mesh.getElements();
  node_map& nodes = mesh.getNodesNC();
  std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();

  std::unordered_set<MEl*> active_elements;
  std::map<MNode*,gind> active_nodes;
  
  MEl* debel = idealElements.begin()->first;
  const ShapeFunctionMatrices* sf = 
    sf_factory.getShapeFunction(debel->getElementType(),
				debel->getOrder(),0);

  ActiveMEl activeEl(debel,index_factory,nodes,sf);
  OptEl3D optel(activeEl,idealElements.at(debel));


 
  for(auto el = idealElements.begin(); el != idealElements.end(); ++el){
    active_elements.insert(el->first);
    //if(all_merits[el->first] > threshold) active_elements.insert(el->first);
  }

  for(int lay=0; lay<=Nlayers; ++lay){
    
    for(auto el = active_elements.begin(); el != active_elements.end(); ++el){
      const int ncn = (*el)->numCornerNodes();
      const gind* cn = (*el)->getCornerNodes();
      for(int i=0; i<ncn; i++){
	active_nodes[nodes.at(cn[i]).get()] = 0;
      }
    }
 
    for(auto el = elements.begin(); el != elements.end(); ++el){
   
      const int ncn = (*el)->numCornerNodes();
      const gind* cn = (*el)->getCornerNodes();
      bool active = false;
      for(int i=0; i<ncn; i++){
	auto it = active_nodes.find(nodes.at(cn[i]).get());
	if(it != active_nodes.end()){
	  active = true;
	  break;
	}
      }
      if(active){
	active_elements.insert(el->get());
	auto it = idealElements.find(el->get());
	if(it == idealElements.end()){
	  ActiveMEl activeEl(el->get(),index_factory,nodes,NULL);
	  arma::mat ideal = activeEl.getNodesMatrix();
	  idealElements.insert(it,std::make_pair(el->get(),ideal));
	}
      }
    }
  }

  for(auto el = active_elements.begin(); el != active_elements.end(); ++el){
    const int ncn = (*el)->numCornerNodes();
    const gind* cn = (*el)->getCornerNodes();
    for(int i=0; i<ncn; i++){
      std::vector<MEl*>& temp = Node2ElementMap[cn[i]];
      temp.push_back(*el);
    }
  }
  */
}

double SubmeshOptimizer::OptimizeSubmesh(gind nd, std::vector<MEl*> elements){

  //SubmeshMeritEvaluator evaluator(nd, elements,mesh);
  //ObjectiveFunction objective_function(&merit_evaluator);
  //ObjectiveFunction.run();

  /*
  arma::mat gradMerit;
  node_map& nodes = mesh.getNodesNC();
  ideal_map& idealElements = mesh.getIdealElementsNC();

  
  arma::mat Grad(3,1,arma::fill::zeros);
  for(auto el = elements.begin(); el != elements.end(); ++el){
      const int ncn = (*el)->numCornerNodes();
      const gind* cn = (*el)->getCornerNodes();
      int grad_node = -1;
      for(int i=0; i<ncn; i++){
	if(cn[i] == nd) grad_node = i;
      }
      assert(grad_node == -1);

    const ShapeFunctionMatrices* sf = 
      sf_factory.getShapeFunction((*el)->getElementType(),
				  (*el)->getOrder(),0);
	
    ActiveMEl activeEl(*el,index_factory,nodes,sf);
    OptEl3D optel(activeEl,idealElements.at(*el));

    double mer = optel.computeGradMerit(gradMerit,0.0);
    Grad.col(0)+= gradMerit.col(grad_node);

  }
  */
}
int BoundaryLayerGenerator::OptimizeBLSubmesh(){

  //SubmeshOptimizer submesh_optimizer(mesh);
  //submesh_optimizer.ComputeNode2ElementMap(1,1.0);


  return 1;
}
