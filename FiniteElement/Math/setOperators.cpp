#include "setOperators.h"
using namespace std;
using namespace arma;


bool pair_equal(pair<int,int> pr1, pair<int,int> pr2){ 
  return (pr1.first == pr2.first); 
}

bool pair_greater(pair<int,int> pr1, pair<int,int> pr2){ 
  return (pr1.first < pr2.first); 
}


/*
uvec set_oper::unique(void){
  int n = x.n_rows;

  uvec unq;

  vector< pair<int,int> > vals(n);
  for(int i=0; i<n; i++){
    vals[i].first = x(i);
    vals[i].second = i;
  }
 

  std::sort(vals.begin(),vals.end(),pair_greater);

  ib.resize(n);
  int cnt = 0;
  ib(vals[0].second) = cnt;
  for(int i=1; i<n; i++){
    if(vals[i-1].first != vals[i].first) cnt++;
    ib[vals[i].second] = cnt;
  }

  std::vector< pair<int,int> >::iterator it;
  it = std::unique(vals.begin(),vals.end(),pair_equal);

  vals.resize(distance(vals.begin(),it));
  int nun = vals.size();
  unq.resize(nun);
  ia.resize(nun);
 
  for(int i=0; i<nun; i++){
    unq(i) = vals[i].first;
    ia(i) = vals[i].second;
  }

  return unq;
}
*/

uvec ismember(uvec a, uvec b){

  
  uint Na = a.n_elem;
  uint Nb = b.n_elem;

  vector< intpair > c(Na);
  for(int i=0; i<Na; i++){
    c[i].first = a[i];
    c[i].second = i;
  }

  // Finds indices of a which are members of b
  stduvec d = conv_to< stduvec >::from(b);

  sort(c.begin(),c.end(),pair_greater);

  sort(d.begin(),d.end());

  uvec mem(Na);

  int memcnt = 0, j=0;
  for(int i=0; i<Nb; i++){
    while(c[j].first <= d[i]){
      if(c[j].first == d[i]){
	mem(memcnt) = c[j].second;
	memcnt++;
      }
      j++;
    }
  }
  mem.resize(memcnt);
  return mem;
  
}

typedef std::pair<int,arma::vec> pair_type;


bool arma_pair_equal(const pair_type &p1, const pair_type &p2){
  double av_norm = (norm(p1.second,2.0) + norm(p2.second,2.0))/2.0;
  if(norm(p2.second-p1.second) < 1.0e-12*av_norm) return true;
  else return false;
}

bool arma_pair_greater(const pair_type &p1, const pair_type &p2){
  
  const int sv = p1.second.n_rows;
  if(p1.second.n_elem!= p2.second.n_elem){
    cout << "elements must be equal!" << endl;
    cout << "num elemes: " << p1.second.n_rows << " " << 
      p2.second.n_rows << endl;
    cout << p1.second << endl;
    cout << p2.second << endl;
  }
  for(int i=0; i<sv; i++){
    if(p2.second(i) > p1.second(i)) return true;
    else if(p2.second(i) == p1.second(i)){}
    else return false;
  }
  return false;
  
}


arma::mat unique(const arma::mat &x,arma::uvec &ia, arma::uvec &ib){
  const int n = x.n_cols;
  const int sv = x.n_rows;

  std::vector<pair_type> values(n);
  for(int i=0; i<n; i++){
    values[i].first = i;
    values[i].second = x.col(i);
  }

  std::sort(values.begin(),values.end(),arma_pair_greater);

  arma::mat uniquex(sv,n);

  ib.resize(n);
  ib(values[0].first) = 0;
  uniquex.col(0) = values[0].second;
  int cnt=0;
  for(int i=1; i<n; i++){
    if(!arma_pair_equal(values[i-1],values[i])){
      cnt++;
      uniquex.col(cnt) = values[i].second;
    }
    ib(values[i].first) = cnt;
  }
  uniquex.resize(sv,cnt+1);


  return uniquex;


}
