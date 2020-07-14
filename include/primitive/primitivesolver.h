#ifndef PRIMITIVESOLVER_H
#define PRIMITIVESOLVER_H

#include <fluidconfig.h>
#include <types.h>

/********************************************************************************
 *
 * Name: PrimitiveSolver
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 *
 * Description: An abstract class representing a primitive solver for a specific
 *              EOS for a relativistic fluid.
 *
 *******************************************************************************/

class Metric;

class PrimitiveSolver{
  protected:
    // Variable labels
    const unsigned int NU    = 5;
    const unsigned int U_D   = 0;
    const unsigned int U_SX  = 1;
    const unsigned int U_SY  = 2;
    const unsigned int U_SZ  = 3;
    const unsigned int U_TAU = 4;

    const unsigned int NV    = 5;
    const unsigned int V_RHO = 0;
    const unsigned int V_VX  = 1;
    const unsigned int V_VY  = 2;
    const unsigned int V_VZ  = 3;
    const unsigned int V_P   = 4;

    Metric* metric;

    double vacuum;

    #ifdef DEBUG_PRIMITIVE_SOLVER
    unsigned int totalIterations;
    unsigned int calls;
    #endif
  public:
    PrimitiveSolver(Metric* m);
    virtual ~PrimitiveSolver();
    
    double calcW(double vu[3]);
    
    /**
     * Calculate the conserved variables from the primitives at a point.
     */
    void primToConPt(double *u, double *v);

    // Abstract methods.
    /**
     * Calculate the specific enthalpy. This is EOS-dependent and therefore virtual.
     */
    virtual double calcH(double rho, double P)=0;
    /**
     * Convert the primitive variables to the conserved variables at a point.
     * This transformation is EOS-dependent, so it is virtual.
     */
    virtual bool conToPrimPt(double *u, double *v)=0;
    /**
     * The sound speed is also dependent on the EOS, so it must be implemented
     * by the child class.
     */
    virtual double calcSoundSpeed(double rho, double p)=0;

    inline void setMetric(Metric* m){
      metric = m;
    }
    inline Metric* getMetric(){
      return metric;
    }

    inline double getVacuum() const{
      return vacuum;
    }
    inline void setVacuum(double v){
      vacuum = v;
    }

    #ifdef DEBUG_PRIMITIVE_SOLVER
    inline void resetNewtonAverage(){
      totalIterations = 0;
      calls = 0;
    }
    inline double getNewtonAverage(){
      return (double)totalIterations/(double)calls;
    }
    #endif
};

#endif
