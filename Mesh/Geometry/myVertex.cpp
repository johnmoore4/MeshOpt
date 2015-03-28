#include "myVertex.h"
#include "armadillo"
#include <assert.h>

using namespace std;
using namespace arma;

void myVertex::paramsOnGeo(const Geo& geom, const double* geo_params, 
			   double* eval_params){
  assert("Cannot parameterize an entity of dimension <=0");

}
