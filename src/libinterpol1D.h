#ifndef __LIBINTERPOL1D_H__
#define __LIBINTERPOL1D_H__

#include <cmath>

class Interpol1D {
 	public:
  		static double linearInterpolation(const double& d1, const double& d2, const double& weight=1.0);
};


#endif
