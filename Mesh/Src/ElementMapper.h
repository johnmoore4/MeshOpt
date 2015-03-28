#pragma once


class ElementMapper{
public:
  
  unique_element_ptr& getMappedElement(const unique_element_ptr& from);

  int insertElementPair(const unique_element_ptr& from,
			const unique_element_ptr& to);

  int Erase();

private:
  std::unordered_map<unique_element_ptr,unique_element_ptr> elmap[3];
  
};
