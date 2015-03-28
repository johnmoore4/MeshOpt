#pragma once

template<class T1, class T2=int>
std::vector<std::pair<T1,T2> > indexed_sort(T1* begin, T1* end){
  std::vector<std::pair<T1,T2> > sorted(std::distance(begin,end));
  int cnt = 0;
  for(T1* it = begin; it != end; ++it, ++cnt){
    sorted[cnt] = std::make_pair(*it,cnt);
  }
  typedef std::pair<T1,T2> mypair;
  
  std::sort(sorted.begin(),sorted.end(),[](mypair const &a, mypair const &b) { return a.first < b.first; });

  return sorted;
}
