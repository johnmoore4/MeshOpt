#pragma once
#include "GlobalDefines.h"
#include "NodeIndexer.h"

#include <map>
#include <set>
#include <armadillo>

class MNode;
class OptEl;
class MEl;

//typedef std::unique_ptr<MNode> unique_node_ptr;
//typedef std::unique_ptr<MEl> unique_element_ptr;
//typedef std::unique_ptr<OptEl> unique_optel_ptr;

typedef std::shared_ptr<MNode> unique_node_ptr;
typedef std::shared_ptr<MEl> unique_element_ptr;
typedef std::shared_ptr<OptEl> unique_optel_ptr;

class ElementEqual{
 private:
  NodeIndexFactory factory;
  
 public:

  bool operator()(const unique_element_ptr& el1, 
		  const unique_element_ptr& el2) const ;

};


class ElementLess{
 public:
  bool operator()(const unique_element_ptr& el1, 
		  const unique_element_ptr& el2) const;

};

class ElementHasher{
private:
  NodeIndexFactory factory;

public:
  std::size_t operator()(const unique_element_ptr& el) const;

};



//typedef std::set<unique_element_ptr,ElementLess> element_set;

//typedef std::unordered_set<unique_element_ptr,ElementHasher,ElementEqual> 
//  element_set;

//typedef std::unordered_map<gind,unique_node_ptr> node_map;

typedef std::set<unique_element_ptr,ElementLess> element_set;
typedef std::map<gind,unique_node_ptr> node_map;

//typedef std::unordered_map<MEl*, unique_optel_ptr> ideal_map;
typedef std::map<MEl*, arma::mat> ideal_map;

typedef std::map<MEl*,double> ElDouble_map;
