#pragma once

#include "Mel.h"
#include "NodeIndexer.h"

#include <fstream>

class ShapeFunctionMatrices;
class NodeIndexFactory;


class ActiveMEl: public Element{
protected:
  const MEl* el;
  NodeIndexFactory& index_factory;
  const node_map& globalnodes;
  const NodeIndexer* ndind;
  const ShapeFunctionMatrices* _sf;
  
  const std::vector<indtype>& getChildNodes(const int child) const{ 
    return ndind->getChildNodes(child);
  }
  const std::vector<indtype>& getReversedNodes() const{
    return ndind->getReversedNodes();
  }
  const std::vector<indtype>& getInteriorNodes() const{
    return ndind->getInteriorNodes();
  }
  const std::vector<indtype> getOrientedNodes(const int orient) const{
    return ndind->getOrientedNodes(orient);
  }

  inline const int NumChildren() const{ return ndind->NumChildren(); }
  /*
  const gind* getCornerNodes() const{ 
    return el->getCornerNodes();
  }
  */


public:
 ActiveMEl(const MEl* elT, NodeIndexFactory& index, const node_map& nodes, 
	const ShapeFunctionMatrices* sfT):
  el(elT), index_factory(index), globalnodes(nodes), 
    ndind(index_factory.getNodeIndexer(el->getElementType(),el->getOrder())), _sf(sfT){}

  const arma::mat getNodesMatrix() const;

  //const std::vector<gind> getAllNodes() const;

  const std::vector<indtype>& getCornerNodeIndex() const{
    return ndind->getCornerNodes();
  }
  const ShapeFunctionMatrices* getShapeFunctionMatrices() const {
    return _sf;
  }

  inline const int Ndof() const {
    return ndind->Ndof(el->getOrder());
  }

  const MEl* getMeshElement() const{ return el; }
  
  const node_map& getGlobalNodes() const { return globalnodes; }

  void writeGMSH(std::ofstream &out, int num);

  arma::mat getLinearNodes() const;

};
