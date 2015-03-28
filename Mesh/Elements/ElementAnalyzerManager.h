#pragma once

class MeshContainer;
class ShapeFunctionMatricesFactory;
class NodeIndexFactory;
class ElementAnalyzer;
class MEl;

class ElementAnalyzerManager{
 public:
  ElementAnalyzerManager(MeshContainer& mesh, 
			 ShapeFunctionMatricesFactory& sf_factory,
			 NodeIndexFactory& index_factory);
  
  ElementAnalyzer getElementAnalyzer(const MEl* el);
  
 private:
  MeshContainer& mesh;
  ShapeFunctionMatricesFactory& sf_factory;
  NodeIndexFactory& index_factory;

};
