#ifndef RECON_H
#define RECON_H
#include <types.h>

/****************************************************************************
 *
 * Class: Recon
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 * 
 * Description: An abstract class representing a reconstruction method for
 *              an HRSC scheme.
 *
 ***************************************************************************/

class Recon{
  protected:
    const unsigned int stencil;
  public:
    Recon(unsigned int n) : stencil(n){}
    virtual ~Recon(){}

    virtual void reconstruct(const unsigned int n, const double* const RESTRICT u,
                        double* const RESTRICT ul, double* const RESTRICT ur)=0;

    virtual void reconstructPt(const unsigned int i, const unsigned int n, 
                               const double* const RESTRICT u,
                               double* const RESTRICT ul, double* const RESTRICT ur)=0;

    inline unsigned int getStencilSize() const{
      return stencil;
    }
};

#endif
