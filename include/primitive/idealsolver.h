#ifndef IDEAL_SOLVER_H
#define IDEAL_SOLVER_H

#include <primitive/primitivesolver.h>

/****************************************************************************
 *
 * Class: IdealSolver
 * Author: Jacob Fields
 * Date Modified: 22-6-2020
 *
 * Description: A solver for an ideal gas with fixed gamma. This is based
 *              on the solver used in DeBuhr et al. (2018), 
 *              https://doi.org/10.3847/1538-4357/aae5f9
 *              As a general rule, this is a very robust solver, with the
 *              final Newton solve typically taking less than 5 iterations
 *              to achieve an accuracy of 10^-14.
 *
 ***************************************************************************/
class IdealSolver : public PrimitiveSolver{
  protected:
    double gamma;

    double f(double x, double D, double tau, double Ssq);
    void fd(double& fx, double& dfx, double x, double D, double tau, double Ssq);
    void fdd(double& fx, double& dfx, double& ddfx, double x, double D, double tau, double Ssq);

    double tol;
    unsigned int maxIterations;
  public:
    IdealSolver(Metric* m);
    virtual ~IdealSolver();

    virtual double calcH(double rho, double P);

    virtual bool conToPrimPt(double *u, double *v);

    virtual double calcSoundSpeed(double rho, double p);

    inline void setGamma(double g){
      gamma = g;
    }
    inline double getGamma(){
      return gamma;
    }

    inline double getTolerance() const{
      return tol;
    }
    inline void setTolerance(double t){
      tol = t;
    }

    inline unsigned int getMaxIterations(){
      return maxIterations;
    }
    inline void setMaxIterations(unsigned int i){
      maxIterations = i;
    }
};

#endif
