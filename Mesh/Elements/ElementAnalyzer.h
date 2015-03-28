#pragma once
#include <armadillo>
#include <memory>

class CompEl;
class ShapeFunctionMatrices;

class ElementAnalyzer{
 public:
  ElementAnalyzer(std::shared_ptr<CompEl>& compel,
		  const ShapeFunctionMatrices* sf_lin,
		  const arma::mat* ideal=NULL);
  double computeDistortion();
  double Distortion(){ return distortion; };
  double Quality(){ return 1.0/distortion; }

 private:
  std::shared_ptr<CompEl> compel;
  double distortion;
  const arma::mat* ideal;
  const ShapeFunctionMatrices* sf_lin;

};
