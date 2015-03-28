#include <vector>
#include <memory>
#include <iostream>

class A{

public: 
  A(){};
  double dbl[20];
};

typedef std::shared_ptr<A> A_ptr;

class B{
  public:
  const std::vector<A_ptr> createAVector(){
    std::vector<A_ptr> vec;
    for(int i=0; i<4; i++){
      vec.push_back(A_ptr( new A() ));
    }
    return vec;
  }
};

int myfunc(){

  // Do Stuff...

  std::vector<A_ptr> globvec;

  B b;
  for(int i=0; i<1e6; i++){
    const std::vector<A_ptr> locvec = b.createAVector();

    for(int i=0; i<locvec.size(); i++) globvec.push_back(locvec[i]);

  }

  globvec.clear();
  globvec.shrink_to_fit();

  // Do more stuff...

  return 1;
}


int main(){

  //std::vector<int> deb;
  //for(int i=0; i<10; i++) deb.push_back(i);

  std::vector<int> deb(10);
  for(int i=0; i<10; i++) deb[i] = i;
  for(auto it=deb.begin(); it < deb.end(); it++){
    std::cout << *it << std::endl;
  }

  myfunc();

  for(auto i=0; i<3; i++){
    myfunc();
  }
  return 1;
}

