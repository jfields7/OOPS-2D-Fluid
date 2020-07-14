#ifndef NO_RECON_H
#define NO_RECON_H
#include <recon/recon.h>

/********************************************************************************
 *
 *
 *
 *******************************************************************************/

class NoRecon : public Recon{
  public:
    NoRecon();
    virtual ~NoRecon();

    virtual void reconstruct(const unsigned int n, const double* const RESTRICT u,
                        double* const RESTRICT ul, double* const RESTRICT ur);
};

#endif
