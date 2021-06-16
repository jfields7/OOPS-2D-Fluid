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
/**
 * The method used for the primitive solver here is based on
 * that present in
 * Kastaun et al., Phys. Rev. D 103, 023018 (2021)
 * It's intended for MHD, but it should work for the special
 * relativistic limit, too.
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
  double D = upt[U_D];
  double tau = upt[U_TAU];
  double Sd[3];
  Sd[0] = upt[U_SX];
  Sd[1] = upt[U_SY];
  Sd[2] = upt[U_SZ];
  double Ssq = metric->squareForm(Sd);
  // Some variables that show up in neutron star MHD calculations
  // but aren't important for us.
  double Ye = 0;
  double Bu[3] = {0.0, 0.0, 0.0};

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

  // Define some auxiliary quantities
  // Note that we don't have a magnetic field, so B = 0.
  // Similarly, we don't have any additional chemical or
  // nuclear species, so Y_i = 0.
  double sqrtD = sqrt(D);
  double bu[3] = {Bu[0]/sqrtD, Bu[1]/sqrtD, Bu[2]/sqrtD};
  double rl[3] = {Sd[0]/D, Sd[1]/D, Sd[2]/D};

  double ru[3] = metric->raiseForm(rl);
  double rsqr  = metric->squareForm(ru);
  double rb    = rl[0]*bu[0] + rl[1]*bu[1] + rl[2]*bu[2];
  double rbsqr = rb*rb;
  double bsqr  = metric->squareVector(bu);
  double q     = tau/D;
  double Ye0  = Ye/D;

  /** This is stuff that doesn't have a purpose in this
      code, but will be needed in Athena++.
  // Kastaun's code complains if we have NaNs. We won't
  // do that here because there shouldn't be NaNs coming in.
  // But it's not a bad idea to check.

  // Check if the magnetic field is too large. Complain if it is.
  // (Obviously that won't be a problem here.)

  // Restrict the range of Ye or other chemical species.
  **/


  // Construct the bracket.
  

  // Check for bizarre cases, like rho being too small or too large.
  // This is something that could be handled by a special policy object.

  // Do the root solve.

  // Complain if we couldn't find a root.

  // Apply the floor?

  // Retrieve the primitive variables.

  // Correct the conserved variables as needed.
}
// }}}
