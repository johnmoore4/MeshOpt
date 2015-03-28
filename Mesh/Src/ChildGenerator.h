#pragma once
#include "NodeIndexer.h"
#include <memory>
#include <vector>
#include "MeshTypedefs.h"

class MEl;
//typedef std::unique_ptr<MEl> element_ptr;
typedef unique_element_ptr element_ptr;

class ChildGenerator{
 private:
  std::vector<element_ptr> children;
  NodeIndexFactory index_factory;

 public:
  
std::vector<element_ptr>& GenerateChildren(const element_ptr& el);
  
};
