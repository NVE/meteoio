#include "libinterpol1D.h"

using namespace std;

double Interpol1D::linearInterpolation(const double& d1, const double& d2, const double& weight){
  double tmp = abs(d1 - d2);
  
  if (d1 < d2){
    return (d1 + tmp*weight);
  } else {
    return (d1 - tmp*weight); 
  }
  return 0;
}
