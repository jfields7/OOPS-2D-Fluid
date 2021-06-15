#include <primitive/generalsolver.h>
#include <cmath>
#include <geom/metric.h>
#include <iostream>
#include <mpicommunicator.h>

// Constructor {{{
GeneralSolver::GeneralSolver(Metric* m) : PrimitiveSolver(m){
  maxIterations = 50;
  tol = 1e-14;
#ifdef DEBUG_PRIMITIVE_SOLVER
  calls = 0;
  totalIterations = 0;
#endif
}
// }}}

// Destructor {{{
GeneralSolver::~GeneralSolver(){

}
// }}}

// {{{
double GeneralSolver::calcH(double rho, double P){
  return rho*(1.0 + calcEps(rho, P)) + P;
}
// }}}

// {{{
double GeneralSolver::calcA(double rho, double eps){
  return calcP(rho, eps)/(rho
}
// }}}

// {{{
/**
 * The method used for the primitive solver here is based on
 * that present in
 * Galeazzi et al., Phys. Rev. D 88, 064009 (2013).
 * It should work for a generic SRHD equation of state.
 */
bool GeneralSolver::conToPrimPt(double *u, double *v){
  const double VELOCITY_FACTOR = 0.9999;
  double upt[NU];
  double sdetg = metric->getDeterminantRoot();
  double isdetg = 1.0/(sdetg + 1.0e-15);
  unsigned int errors = 0;

  // First undensitize the variables
  upt[U_D  ] = isdetg * u[U_D  ];
  upt[U_SX ] = isdetg * u[U_SX ];
  upt[U_SY ] = isdetg * u[U_SY ];
  upt[U_SZ ] = isdetg * u[U_SZ ];
  upt[U_TAU] = isdetg * u[U_TAU];

  // Some utility variables for calculation.
  double Sd[3];
  Sd[0] = upt[U_SX];
  Sd[1] = upt[U_SY];
  Sd[2] = upt[U_SZ];
  double Ssq = metric->squareForm(Sd);

  // Apply the floor to the undensitized variables
  if (upt[U_D] < vacuum){
    upt[U_D  ] = vacuum;
    upt[U_SX ] = 0.0;
    upt[U_SY ] = 0.0;
    upt[U_SZ ] = 0.0;
    Ssq = 0.0;
    errors++;

    v[V_RHO] = D;
    v[V_VX ] = 0.0;
    v[V_VY ] = 0.0;
    v[V_VZ ] = 0.0;
  }

  double q = upt[U_TAU]/upt[U_D];
  double r = sqrt(Ssq)/upt[U_D];
  double k = r/(1.0 + q);

  // Check some physical bounds.
  // If q < 0, we have a negative energy density.
  if(q < 0){
    // Set q = 0, but leave k and D constant. This effectively rescales
    // r and thus Ssq.
    q = 0;
    r = k;
    Ssq = r*upt[U_D];
  }
  double k_max = 2.0*VELOCITY_FACTOR/(1.0 + VELOCITY_FACTOR*VELOCITY_FACTOR);
  // If we exceed the maximum allowed k, rescale k and the momentum density.
  if(k >= k_max){
    k = k_max;
    r = k_max*(1.0 + q);
    Ssq = r*upt[U_D];
  }
}
// }}}
