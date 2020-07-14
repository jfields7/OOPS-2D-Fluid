#ifndef ENTHALPY_SOLVER_H
#define ENTHALPY_SOLVER_H

#include <primitive/idealsolver.h>

/****************************************************************************
 *
 * Class: Enthalpy Solver 
 * Author: Jacob Fields
 * Date Modified: 25-6-2020
 *
 * Description: A solver for an ideal gas with fixed gamma. This is based
 *              on the solver used in
 *              Kastaun, W. 2007, Ph.D. Thesis, University of Tübingen
 *              It is based around rewriting the conserved and primitive
 *              variables in terms of specific enthalpy alone, then solving
 *              performing a root solve until the predicted enthalpy matches
 *              the calculated enthalpy.
 *
 * Known Issues: This solver is not nearly as fast as IdealSolver. It has
 *               the theoretical advantage that W is calculated from a
 *               quantity that can never be complex, but the lower bound
 *               is too conservative. One possible improvement might be
 *               scaling h by 1/(tau + D), which hypothetically could make 
 *               the quantity closer to order 1.
 *
 ***************************************************************************/

class EnthalpySolver : public IdealSolver{
  protected:
    double calcRho(double D, double W);
    double calcP(double D, double tau, double h, double W);

  public:
    EnthalpySolver(Metric* m);
    virtual ~EnthalpySolver();

    virtual bool conToPrimPt(double *u, double *v);
};

#endif
