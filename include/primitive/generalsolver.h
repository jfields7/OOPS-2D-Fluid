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
 * Reference:   Kastaun et al., Phys. Rev. D 103, 023018 (2021)
 *
 **************************************************************/

class GeneralSolver : public PrimitiveSolver {
  protected: 
    unsigned int maxIterations;
    double tol;

    double maxRho;
    double maxEps;

    double minRho;
    double minEps;
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

    inline double getMinRho() const{
      return minRho;
    }
    inline void setMaxRho(double rho){
      minRho = rho;
    }

    inline double getMaxEps() const{
      return maxEps;
    }
    inline void setMaxEps(double eps){
      maxEps = eps;
    }

    inline double getMinEps() const{
      return minEps;
    }
    inline void setMinEps(double eps){
      minEps = eps;
    }
};

#endif
