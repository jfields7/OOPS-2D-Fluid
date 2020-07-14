#ifndef MINMOD_H
#define MINMOD_H
#include <cmath>
#include <recon/recon.h>

/*************************************************************************
 *
 * Class: Minmod
 * Author: Jacob Fields
 * Date Modified: 16-6-2020
 *
 * Description: Reconstruction via the minmod slope-limiting algorithm.
 *
 ************************************************************************/

class Minmod : public Recon{
  private:
    inline double minmod(double x, double y){
      return 0.5*(copysign(1.0,x) + copysign(1.0,y)) * fmin(fabs(x),fabs(y));
    }
  public:
    Minmod();
    virtual ~Minmod();

    virtual void reconstruct(const unsigned int n, const double* const RESTRICT u,
                             double* const RESTRICT ul, double* const RESTRICT ur);
};

#endif
