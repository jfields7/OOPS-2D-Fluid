#include <primitive/generalsolver.h>
#include <cmath>
#include <geom/metric.h>
#include <iostream>
#include <mpicommunicator.h>
#include <numtoolsroot.h>

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

// conToPrimPt {{{
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
  double rsq  = metric->squareForm(ru);
  double rb    = rl[0]*bu[0] + rl[1]*bu[1] + rl[2]*bu[2];
  double rbsq = rb*rb;
  double bsq  = metric->squareVector(bu);
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
  // This solver looks for a solution in terms of mu = 1/Wh. The lower
  // bound for mu is trivially 0, corresponding to either a high-energy,
  // high-pressure solution or an ultrarelativistic solution.
  // If r < minH, where minH is the smallest enthalpy in the system, the 
  // upper bound is just 1/minH. This corresponds to the Newtonian limit.
  // If r > minH, we can get a tighter bound by finding the root to the
  // auxiliary function, fa(mu) = mu*sqrt(minH^2 + rbar(mu)^2) - 1.
  // It's not important that we know the upper-bound super precisely,
  // so we're satisfied with a fixed tolerance of 1e-3 and 10 iterations.
  // However, this upper root is important because there is a kink beyond
  // it that makes the main root finding more difficult.
  double mul = 0.0;
  double muh = 1.0/minH;
  if(rsq > minH*minH){
    double btol = 1e-3;
    double max_iter = 10;
    double count = 0;
    double f;
    do {
      double x = 1.0/(1.0 + muh*bsq);
      double xsq = x*x;
      double rbarsq = rsq*xsq + muh*x*(1.0 + x)*rbsq;
      // The function fa is well behaved, so we can solve for it using a
      // derivative-based solution. For computational efficiency, we
      // need to calculate some extra quantities for both f and the
      // derivative df.
      double dis = sqrt(h0*h0 + rbarsq);
      double dx = -bsqr*xsq;
      double drbarsq = rbsq*x*(1.0 + x) + (muh*rbsq + 2.0*(muh*rbsq + rsq)*x)*dx;
      f = muh*dis - 1.0;
      double df = dis + mu*drbarsq/(2.0*dis);
      muh = muh - f/df;
      count++;
    }
    while (fabs(f) > btol && count < max_iter);
    // Throw an error if we haven't converged.
    if (fabs(f) > btol){
      errors++;
      printf("C2P: Having difficulty calculating the upper root.");
    }
  }
  

  // Kastaun checks for some nonphysical cases using the bracketed root afterward,
  // though I'm not sure why. I see no reason why we should expect the primitive
  // variables to have a physical state after bracketing, so it just seems like
  // wasted computational effort to me.

  // Do the root solve.
  // Again, we're solving for mu = 1/Wh. The actual master function used for the root
  // solve is fairly complicated. For testing purposes, it is written out explicitly
  // here, but it would be wise to wrap it up in a few functions for the sake of
  // appearances.
  // The root solve here is performed with the Illinois variant of false position,
  // but it would be good to experiment with different solvers.


  // Complain if we couldn't find a root.

  // Apply the floor?

  // Retrieve the primitive variables.

  // Correct the conserved variables as needed.
}
// }}}
