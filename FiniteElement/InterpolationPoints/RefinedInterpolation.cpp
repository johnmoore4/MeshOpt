#include "RefinedInterpolation.h"
#include "InterpolationPoints.h"
#include "PolynomialBasis.h"

#include <set>

template <typename T1>
std::array<T1,3> operator + (std::array<T1,3>& a, 
			     std::array<T1,3>& b){

  std::array<T1,3> arr;
  for(int i = 0; i < 3; i++) arr[i] = a[i]+b[i];
  return arr;
}

template <typename T1>
std::array<T1,3> operator / (const std::array<T1,3>& a, T1 b){

  std::array<T1,3> arr;
  for(int i = 0; i < 3; i++) arr[i]  = a[i]/b;
  return arr;
}

RefinedInterpolationFactory::RefinedInterpolationFactory(){
  interpolations.resize(8);
  interpolation_factory = std::make_shared<InterpolationPointFactory>();
  basis_factory = std::make_shared<PolynomialBasisFactory>();
}

RefinedInterpolation& RefinedInterpolationFactory::
getRefinedInterpolation(int Nref, int p,int eltype){
  if(interpolations[eltype]){
    return *interpolations[eltype];
  }
  else{
    //std::cout << "creating refinement for type: " << eltype << std::endl;
    const InterpolationPoints& points = 
      interpolation_factory->getInterpolationPoints(eltype);
    const PolynomialBasis& basis = basis_factory->getPolynomialBasis(eltype);

    switch(eltype){
    case 2:
      interpolations[eltype] = 
	std::make_shared<TriRefinedInterpolation>(Nref,p,basis,points);
      break;
    case 3:
      interpolations[eltype] = 
	std::make_shared<QuadRefinedInterpolation>(Nref,p,basis,points);
      break;
    case 4:
      interpolations[eltype] =
	std::make_shared<TetRefinedInterpolation>(Nref,p,basis,points);
      break;
    case 5:
      throw std::runtime_error("Hex refinement not implemented!");
      break;
    case 6:
      interpolations[eltype] =
	std::make_shared<PrismRefinedInterpolation>(Nref,p,basis,points);
      break;
    }
    //std::cout << "before computing refinement for : " << eltype << std::endl;
    interpolations[eltype]->Refine(Nref);
    //std::cout << "after refinement for " << eltype << std::endl;
    return *interpolations[eltype];
  }
}

int RefinedInterpolation::Refine(int Nref){

  refinedPointsVec.resize(refinedPoints.size());
  for(auto it = refinedPoints.begin(); it != refinedPoints.end(); ++it){
    refinedPointsVec[it->second] = it->first;
  }

  connectivity.resize(1,nnode);

  for(int i = 0; i < nnode; i++) connectivity[i] = i;

  for(int ref = 0; ref < Nref; ref++){
    int Nel_curr = connectivity.n_rows;
    arma::umat connectivity_new(Nel_curr*nsub,nnode);
    for(int el = 0; el < connectivity.n_rows; el++){
      int newpt_indices[20];
      for(int i = 0; i < nnode; i++) newpt_indices[i] = connectivity(el,i);
      for(int eg = 0; eg < edgesList.n_cols; eg++){
	int e1 = connectivity(el,edgesList(0,eg));
	int e2 = connectivity(el,edgesList(1,eg));

	std::array<double,3> newpt = 
	  (refinedPointsVec[e1] + refinedPointsVec[e2])/2.0;


	auto it = refinedPoints.find(newpt);
	if(it == refinedPoints.end()){
	  it = refinedPoints.insert(it,std::make_pair(newpt,point_counter++));
	  refinedPointsVec.push_back(newpt);
	}
	newpt_indices[nnode+eg] = it->second;
	
      }


      for(int se = 0; se < nsub; se++){
	for(int i = 0; i < nnode; i++){
	  connectivity_new(el*nsub+se,i) = 
	    newpt_indices[recursionConnectivity[se][i]];
	}
      }
    }
    connectivity = connectivity_new;

  }
  // end of refinement loop
  const arma::mat& nodal_pts = points.ComputePoints(p,0);

  arma::mat refined_pts(refinedPointsVec.size(),3);
  for(int pt = 0; pt < refinedPointsVec.size(); pt++){
    for(int i = 0; i < 3; i++) refined_pts(pt,i) = refinedPointsVec[pt][i];
  }
  /*
  if(connectivity.n_cols == 4){
    std::cout << connectivity << std::endl;
    std::cout << refined_pts << std::endl;
  }

  if(connectivity.n_cols == 6){
    std::cout << connectivity << std::endl;
    std::cout << refined_pts << std::endl;
  }
  */
  //connectivity.save("connectivity",arma::raw_ascii);
  //refined_pts.save("refined_pts",arma::raw_ascii);

  
  arma::mat A, f;
  arma::cube dA, df;
  arma::cube dA2, df2;
  basis.EvalBasis(p,nodal_pts,A,dA,dA2);
  basis.EvalBasis(p,refined_pts,f,df,df2);


  interpolationMatrix = arma::solve(A,f);


  
}

TriRefinedInterpolation::
TriRefinedInterpolation(int Nref,
			int p,
			const PolynomialBasis& basis,
			const InterpolationPoints& points):
  RefinedInterpolation(Nref,3,3,4,p,basis,points){
  
  refinedPoints[{0,0,0}] = point_counter++;
  refinedPoints[{1,0,0}] = point_counter++;
  refinedPoints[{0,1,0}] = point_counter++;

  
  edgesList = {0,1,1,2,2,0};
  edgesList.reshape(2,3);


  recursionConnectivity.push_back({0,3,5});
  recursionConnectivity.push_back({3,1,4});
  recursionConnectivity.push_back({4,2,5});
  recursionConnectivity.push_back({3,4,5});

}

QuadRefinedInterpolation::
QuadRefinedInterpolation(int Nref,
			int p,
			const PolynomialBasis& basis,
			const InterpolationPoints& points):
  RefinedInterpolation(Nref,4,5,4,p,basis,points){
  
  refinedPoints[{0,0,0}] = point_counter++;
  refinedPoints[{1,0,0}] = point_counter++;
  refinedPoints[{0,1,0}] = point_counter++;
  refinedPoints[{1,1,0}] = point_counter++; 

  edgesList = {0,1,1,3,3,2,2,0,1,2};
  edgesList.reshape(2,nadd);

  recursionConnectivity.push_back({0,4,7,8});
  recursionConnectivity.push_back({4,1,8,5});
  recursionConnectivity.push_back({8,5,6,3});
  recursionConnectivity.push_back({7,8,2,6});


}

TetRefinedInterpolation::
TetRefinedInterpolation(int Nref,
			int p,
			const PolynomialBasis& basis,
			const InterpolationPoints& points):
  RefinedInterpolation(Nref,4,6,8,p,basis,points){

  refinedPoints[{0,0,0}] = point_counter++;
  refinedPoints[{1,0,0}] = point_counter++;
  refinedPoints[{0,1,0}] = point_counter++;
  refinedPoints[{0,0,1}] = point_counter++;

  edgesList = {0,1,1,2,2,0,0,3,1,3,2,3};
  edgesList.reshape(2,nadd);

  recursionConnectivity.push_back({0,4,6,7});
  recursionConnectivity.push_back({4,1,5,8});
  recursionConnectivity.push_back({5,2,6,9});
  recursionConnectivity.push_back({7,8,9,3});
  recursionConnectivity.push_back({8,9,4,7});
  recursionConnectivity.push_back({9,6,4,7});
  recursionConnectivity.push_back({6,9,4,5});
  recursionConnectivity.push_back({9,8,4,5});

  /*
  recursionConnectivity.push_back({7,9,8,0});
  recursionConnectivity.push_back({4,1,5,8});
  recursionConnectivity.push_back({5,2,6,9});
  recursionConnectivity.push_back({7,8,9,3});
  recursionConnectivity.push_back({8,9,4,0});
  recursionConnectivity.push_back({9,6,4,0});
  recursionConnectivity.push_back({6,9,4,5});
  recursionConnectivity.push_back({9,8,4,5});
  */

}

PrismRefinedInterpolation::
PrismRefinedInterpolation(int Nref,
			  int p,
			  const PolynomialBasis& basis,
			  const InterpolationPoints& points):
  RefinedInterpolation(Nref,6,12,8,p,basis,points){

  refinedPoints[{0,0,0}] = point_counter++;
  refinedPoints[{1,0,0}] = point_counter++;
  refinedPoints[{0,1,0}] = point_counter++;
  refinedPoints[{0,0,1}] = point_counter++;
  refinedPoints[{1,0,1}] = point_counter++;
  refinedPoints[{0,1,1}] = point_counter++;

  edgesList = {0,1,1,2,2,0,3,4,4,5,5,3,0,3,1,4,2,5,0,4,1,5,2,3};
  edgesList.reshape(2,nadd);

  recursionConnectivity.push_back({0,6,8,12,15,17});
  recursionConnectivity.push_back({6,1,7,15,13,16});
  recursionConnectivity.push_back({8,7,2,17,16,14});
  recursionConnectivity.push_back({6,7,8,15,16,17});
  recursionConnectivity.push_back({12,15,17,3,9,11});
  recursionConnectivity.push_back({15,13,16,9,4,10});
  recursionConnectivity.push_back({17,16,14,11,10,5});
  recursionConnectivity.push_back({15,16,17,9,10,11});

}
