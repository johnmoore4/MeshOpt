#pragma once
#include <utility>
#include <algorithm>
#include <vector>
#include "armadillo"

template <class T>
class set_oper{
public:
  T x;
  arma::uvec ia, ib;
 set_oper(T xT): x(xT){};
  T unique(void);
  T sort(void);
  static bool pair_equal(std::pair<T,int> pr1, std::pair<T,int> pr2);
  static bool pair_greater(std::pair<T,int> pr1, std::pair<T,int> pr2);

};



using namespace std;
using namespace arma;

template <class T>
bool set_oper<T>::pair_equal(pair<T,int> pr1, pair<T,int> pr2){ 
  int nel = pr1.first.n_elem;
  int cnt = 0;
  for(int i=0; i<nel; i++){
    if( abs(pr1.first[i] - pr2.first[i]) < max(1.0e-10*(min(abs(pr1.first[i]), abs(pr2.first[i]))), 1.0e-16) ) cnt++;
    else return false;
  }
  // If all elements are equal, return true
  return true;

}

template <class T>
bool set_oper<T>::pair_greater(pair<T,int> pr1, pair<T,int> pr2){ 
  int nel = pr1.first.n_elem;
  int nel2 = pr2.first.n_elem;
  if(nel2 != nel) cout << "ERROR! nel2 != nel" << endl;

  for(int i=0; i<nel; i++){
    if(pr1.first[i] < pr2.first[i]) return true;
    else if(pr1.first[i] > pr2.first[i]) return false;
  }
  // If they are equal, return false
  return false;
}


template <class T>
T set_oper<T>::sort(void){
  T srt;
  
 
  int n = x.n_rows;
  int m = x.n_cols;
  srt.resize(n,m);

  ia.zeros(n);

  vector< pair<T,int> > vals(n);
  for(int i=0; i<n; i++){
    vals[i].first = x.row(i);
    vals[i].second = i;
  }

  std::sort(vals.begin(),vals.end(),pair_greater);

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++) srt(i,j) = vals[i].first(j);
    ia(i) = vals[i].second;
  }

  return srt;
  
}


template <class T>
T set_oper<T>::unique(void){
  int n = x.n_rows;
  int m = x.n_cols;

  T unq;

  vector< pair<T,int> > vals(n);
  for(int i=0; i<n; i++){
    vals[i].first = x.row(i);
    vals[i].second = i;
  }
 

  std::sort(vals.begin(),vals.end(),pair_greater);

  ib.resize(n);
  int cnt = 0;
  ib(vals[0].second) = cnt;
  for(int i=1; i<n; i++){
    if(!pair_equal(vals[i-1],vals[i])) cnt++;
    //if(vals[i-1].first != vals[i].first) cnt++;
    ib[vals[i].second] = cnt;
  }

  typename std::vector< std::pair<T,int> >::iterator it;
  //std::vector< int >::iterator it;
  it = std::unique(vals.begin(),vals.end(),pair_equal);

  vals.resize(distance(vals.begin(),it));
  int nun = vals.size();
  unq.resize(nun,m);
  ia.resize(nun);
 
  for(int i=0; i<nun; i++){
    for(int j=0; j<m; j++) unq(i,j) = vals[i].first(j);
    //unq.row(i) = vals[i].first;
    ia(i) = vals[i].second;
  }

  return unq;
}


typedef std::pair<int,int> intpair;
typedef std::vector< unsigned int > stduvec;
arma::uvec ismember(arma::uvec a, arma::uvec b);

/*
class set_oper{
public:
  arma::uvec x, ia, ib;
  set_oper(arma::uvec xx){ x = xx;}
  arma::uvec unique(void);
  arma::uvec sort(void);
};


arma::uvec ismember(arma::uvec a, arma::uvec b);

typedef std::pair<int,int> intpair;
typedef std::vector< unsigned int > stduvec;
*/

arma::mat unique(const arma::mat &x,arma::uvec &ia, arma::uvec &ib); 
