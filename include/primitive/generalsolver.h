#ifndef GENERAL_SOLVER_H
#define GENERAL_SOLVER_H

#include <primitive/primitivesolver.h>

/***************************************************************
 * 
 * Class: GeneralSolver
 * Author: Jacob Fields
 * Date Modified: 14-6-2021
 *
 * Description: An EOS solver for a general equation of state.
 *              A specific EOS should inherit from this class
 *              and use its primitive solver.
 *
 * Reference:   Galeazzi et al., Phys. Rev D 88, 064009 (2013).
 *
 **************************************************************/

class GeneralSolver : public PrimitiveSolver {
  protected: 
    unsigned int maxIterations;
    double tol;

    double maxRho;
    double maxEps;
  public:
    GeneralSolver(Metric* m);
    virtual ~GeneralSolver();

    virtual double calcH(double rho, double P);
    virtual bool conToPrimPt();
    virtual double calcSoundSpeed(double rho, double p);

    inline getTolerance() const{
      return tol;
    }
    inline void setTolerance(double t){
      tol = t;
    }

    inline unsigned int getMaxIterations() const{
      return maxIterations;
    }
    inline void setMaxIterations(unsigned int i){
      maxIterations = i;
    }

    inline double getMaxRho() const{
      return maxRho;
    }
    inline void setMaxRho(double rho){
      maxRho = rho;
    }

    virtual double calcP(double rho, double eps) = 0;
    virtual double calcEps(double rho, double P) = 0;

    double calcA(double rho, double eps);
};

#endif
