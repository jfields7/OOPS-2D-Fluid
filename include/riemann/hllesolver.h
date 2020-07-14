#ifndef HLLESOLVER_H
#define HLLESOLVER_H

#include <riemann/riemannsolver.h>

/****************************************************************************
 *
 * Class: HLLESolver
 * Author: Jacob Fields
 * Date Modified: 17-6-2020
 *
 * Description: An implementation of the HLLE solver for the Riemann 
 *              problem.
 *
 ***************************************************************************/

class HLLESolver : public RiemannSolver{
  private:
    bool calcMaxSpeeds(double *bp, double *bm, double *vl, double *vr, const int dir);
    void calcDefaultSpeeds(double *bp, double *bm, const int dir); 
  public:
    HLLESolver(PrimitiveSolver* primitive, Metric* m);
    virtual ~HLLESolver();
    virtual void solveRiemannProblem(double *F, double *fl, double *fr,
                                     double *ul, double *ur, double *vl, double *vr,
                                     const int dir);
};

#endif
