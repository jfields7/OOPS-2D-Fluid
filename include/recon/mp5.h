#ifndef MP5_H
#define MP5_H
#include <cmath>
#include <recon/recon.h>

/*************************************************************************
 *
 * Class: MP5
 * Author: Jacob Fields
 * Date Modified: 1-7-2020
 *
 * Description: Reconstruction via MP5.
 *
 ************************************************************************/

class MP5 : public Recon{
  private:
    inline double minmod(double x, double y){
      return 0.5*(copysign(1.0,x) + copysign(1.0,y)) * fmin(fabs(x),fabs(y));
    }
    inline double minmod4(double w, double x, double y, double z){
      return 0.125*(copysign(1.0,w) + copysign(1.0,x))*fabs((copysign(1.0,w)+copysign(1.0,y)) *
                   (copysign(1.0,w) + copysign(1.0,z)))*fmin(fabs(w),fmin(fabs(x),fmin(fabs(y),fabs(z))));
    }

    inline double doMP5(const double am2, const double am1, const double a,
                      const double ap1, const double ap2,
                      const double anorm, const double mp5_eps, const double mp5_alpha){
      const double vl = (2.0*am2 - 13.0*am1 + 47.0*a + 27.0*ap1 - 3.0*ap2) * 1.0/60.0;
      const double vmp = a + minmod(ap1-a, mp5_alpha*(a-am1));
      if((vl-a)*(vl-vmp) < mp5_eps*anorm){
        return vl;
      }
      else{
        const double djm1 = am2 - 2.0*am1 + a;
        const double dj   = am1 - 2.0*a + ap1;
        const double djp1 = a - 2.0*ap1 + ap2;
        const double dm4jph = minmod4(4.0*dj-djm1, 4.0*djm1-dj, dj, djm1);
        const double dm4jmh = minmod4(4.0*dj-djm1, 4.0*djm1-dj, dj, djm1);
        const double vul = a + mp5_alpha*(a-am1);
        const double vav = 0.5*(a+ap1);
        const double vmd = vav - 0.5*dm4jph;
        const double vlc = a + 0.5*(a-am1) + 4.0/3.0*dm4jmh;
        const double vmin = fmax(fmin(a,fmin(ap1,vmd)),fmin(a,fmin(vul,vlc)));
        const double vmax = fmin(fmax(a,fmax(ap1,vmd)),fmax(a,fmax(vul,vlc)));
        return vl + minmod(vmin-vl, vmax-vl);
      }
      return 0;
    }
    const bool adaptive = true;
  public:
    MP5();
    virtual ~MP5();

    virtual void reconstruct(const unsigned int n, const double* const RESTRICT u,
                             double* const RESTRICT ul, double* const RESTRICT ur);

    virtual void reconstructPt(const unsigned int i, const unsigned int n, 
                               const double* const RESTRICT u,
                               double* const RESTRICT ul, double* const RESTRICT ur);
};

#endif
