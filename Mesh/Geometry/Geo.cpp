#include "Geo.h"

void Geo::param2xyz(const double *params, double* xyz) const{
  const gp_Pnt point = Value(params);
  xyz[0] = point.X();
  xyz[1] = point.Y();
  xyz[2] = point.Z();
}
