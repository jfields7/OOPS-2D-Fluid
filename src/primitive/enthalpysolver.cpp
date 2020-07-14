#include <primitive/enthalpysolver.h>
#include <cmath>
#include <geom/metric.h>
#include <iostream>
#include <mpicommunicator.h>

// Constructor {{{
EnthalpySolver::EnthalpySolver(Metric *m) : IdealSolver(m){
  maxIterations = 100;
}
// }}}

// Destructor {{{
EnthalpySolver::~EnthalpySolver(){

}
// }}}

// calcRho {{{
double EnthalpySolver::calcRho(double D, double W){
  return D/W;
}
// }}}

// calcP {{{
double EnthalpySolver::calcP(double D, double tau, double h, double W){
  return D*(h*W - 1.0) - tau;
}
// }}}

// conToPrimPt {{{
bool EnthalpySolver::conToPrimPt(double *u, double *v){
  const double VELOCITY_FACTOR = 0.9999;
  double upt[NU];
  double sdetg = metric->getDeterminantRoot();
  double isdetg = 1.0/(sdetg + 1.0e-15);
  unsigned int errors = 0;

  // Undensitize the variables
  upt[U_D  ] = isdetg * u[U_D  ];
  upt[U_SX ] = isdetg * u[U_SX ];
  upt[U_SY ] = isdetg * u[U_SY ];
  upt[U_SZ ] = isdetg * u[U_SZ ];
  upt[U_TAU] = isdetg * u[U_TAU];

  // Some variables we'll need.
  double Sd[3];
  Sd[0] = upt[U_SX];
  Sd[1] = upt[U_SY];
  Sd[2] = upt[U_SZ];
  double Ssq = metric->squareForm(Sd);

  // Apply the floor to the undensitized variables.
  if(upt[U_D] < vacuum){
    upt[U_D  ] = vacuum;
    upt[U_SX ] = 0.0;
    upt[U_SY ] = 0.0;
    upt[U_SZ ] = 0.0;
    upt[U_TAU] = vacuum;
  }
  else{
    if(upt[U_TAU] < vacuum){
      upt[U_TAU] = vacuum;
    }
    double Ssq_max = (2.0*upt[U_D] + upt[U_TAU])*upt[U_TAU];
    if(Ssq_max < Ssq){
      double t_min = VELOCITY_FACTOR * sqrt(Ssq_max/Ssq);
      upt[U_SX] *= t_min;
      upt[U_SY] *= t_min;
      upt[U_SZ] *= t_min;
      Ssq *= t_min*t_min;
    }
  }

  double D = upt[U_D];
  double tau = upt[U_TAU];

  // We first need to bound h.
  double lb, ub, b, asq;
  asq = Ssq/(D*D);
  b = (tau + D)/D;

  // However, if asq is effectively zero, the transformation is actually analytic.
  if(asq < 1e-16){
    v[V_RHO] = D;
    v[V_P  ] = (gamma - 1.0)*tau;
    v[V_VX ] = 0.0;
    v[V_VY ] = 0.0;
    v[V_VZ ] = 0.0;
  }
  else{
    // For the lower bound to be valid, we need Ssq < (tau + D)^2. This
    // inequality is weaker than the one enforced during the floor, so
    // we are already guaranteed to have a physical value.
    lb = fmax(1.0,sqrt(b*b - asq));
    double g = (gamma - 1.0)/gamma;
    /*double adjb = b - g;
    double k = 1.0 - g*g;
    double dis = adjb*adjb - asq*k;
    if (dis < 0){
      std::cout << "EnthalpySolver: The discriminant for the upper bound is negative.\n";
      return false;
    }
    ub = (-(g*(1.0 - b*g)) + sqrt(dis))/k;*/
    ub = 1.0 + gamma*(b - 1.0);

    // Calculate the bounds.
    double Wub = sqrt(1.0 + asq/(ub*ub));
    double fub = 1.0 + Wub*(ub*Wub - b)/g - ub;
    double Wlb = sqrt(1.0 + asq/(lb*lb));
    double flb = 1.0 + Wlb*(lb*Wlb - b)/g - lb;
    
    // Check that the bounds are good. We initialize h to the midpoint, which
    // we can adjust later on.
    double h = 0.5*(ub + lb);
    if(fub < 0){
      if(fub > -tol){
        // This is close enough to give us a good answer, so we'll just
        // push ub up by a hair and let the solve complete.
        h = ub;
        ub += tol;
        Wub = sqrt(1.0 + asq/(ub*ub));
        fub = 1.0 + Wub*(ub*Wub - b)/g - ub;
      }
      else{
        std::cout << "EnthalpySolver: There's a problem with the upper bound.\n";
        return false;
      }
    }
    else if(fub < tol){
      h = ub;
      ub += tol;
      Wub = sqrt(1.0 + asq/(ub*ub));
      fub = 1.0 + Wub*(ub*Wub - b)/g - ub;
    }
    if(flb > 0){
      if(flb < tol){
        // This is close enough to give us a good answer, so we'll just
        // push lb down by a hair and let the solve complete.
        h = lb;
        lb -= tol; 
        Wlb = sqrt(1.0 + asq/(lb*lb));
        flb = 1.0 + Wlb*(lb*Wlb - b)/g - lb;
      }
      else{
        std::cout << "EnthalpySolver: There's a problem with the lower bound.\n";
      }
    }
    else if(flb > -tol){
      h = lb;
      lb -= tol;
      Wlb = sqrt(1.0 + asq/(lb*lb));
      flb = 1.0 + Wlb*(lb*Wlb - b)/g - lb;
    }
    // By now, our bounds should be good, so we can go ahead with the solve.
    double h0 = h;
    double ub0 = ub;
    double lb0 = lb;
    double fub0 = fub;
    double flb0 = flb;
    double fh = 1.0;
    double dfh = 1.0;
    unsigned int count = 0;
    while(fabs(fh) > tol && count < maxIterations){
      double oldh = h;
      double W = sqrt(1.0 + asq/(h*h));
      double hW = h*W;
      fh = 1.0 + W*(hW - b)/g - h;
      dfh = (1.0 + asq*(b-hW)/(h*h*hW))/g - 1.0;
      h = h - fh/dfh;
      // Check for convergence by ensuring that the new q is
      // still inside the bounds. If it's not, apply a variant
      // of the false position method.
      if((h > ub || h < lb) && fabs(fh) > tol){
        if(fh*fub > 0){
          double m = 1.0;
          h = (fub*lb - m*flb*ub)/(fub - m*flb);
        }
        else{
          double m = 1.0;
          h = (m*fub*lb - flb*ub)/(m*fub - flb);
        }
      }
      // Move the bounds according to the sign of fh.
      if(fh > 0){
        ub = oldh;
        fub = fh;
      }
      else{
        lb = oldh;
        flb = fh;
      }
      count++;
    }
    // In the event of failure, specify what went on.
    if(fabs(fh) > tol){
      printf("EnthalpySolver: The Newton solve failed.\n");
      printf("EntahlpySolver: Current values include \n");
      printf("  lb = %25.20e \n",lb);
      printf("  ub = %25.20e \n",ub);
      printf("  h  = %25.20e \n",h);
      printf("  lb0 = %25.20e \n",lb0);
      printf("  ub0 = %25.20e \n",ub0);
      printf("  h0  = %25.20e \n",h0);
      printf("  flb0 = %25.20e \n",flb0);
      printf("  fub0 = %25.20e \n",fub0);
      printf("EnthalpySolver: Some conserved quantities and \n");
      printf("EnthalpySolver: inequalities, \n");
      printf("  D                     = %25.20e\n", D);
      printf("  (tau+D)/D             = %25.20e\n", b);
      printf("  S^2/D^2               = %25.20e\n", asq);
      printf("  (Ssq+D^2)/(tau + D)^2 = %25.20e\n", (Ssq+D*D)/(tau + D)*(tau + D));

      MPI_Abort(MPI_COMM_WORLD, 0);
    }
    #ifdef DEBUG_PRIMITIVE_SOLVER
    totalIterations += count;
    calls++;
    #endif

    // If we get to this point, we probably have a solution!
    double W = sqrt(1.0 + asq/(h*h));
    double hW = h*W;
    double HWsq = D*hW;

    // We need the "up" components of S.
    double Su[3];
    metric->raiseForm(Su, Sd);

    // Now we can calculate the primitive variables.
    v[V_RHO] = D/W;
    v[V_P  ] = D*(hW - 1.0) - tau;
    v[V_VX ] = Su[0] / HWsq;
    v[V_VY ] = Su[1] / HWsq;
    v[V_VZ ] = Su[2] / HWsq;

    // In some primitive solvers, the pressure ends up identically zero.
    // Let's see if that's an issue here.
    #ifdef DEBUG_PRIMITIVE_SOLVER
    if(v[V_P] == 0){
      std::cout << "Pressure is identically zero.\n";
      std::cout << "Add an analytical solution for this case.\n";
    }
    #endif
  }
  // Apply the floor
  if(v[V_RHO] < vacuum){
    v[V_RHO] = vacuum;
  }
  if(v[V_P] < vacuum){
    v[V_P] = vacuum;
  }

  // Redensitize the variables.
  upt[U_D  ] = sdetg*u[U_D  ];
  upt[U_SX ] = sdetg*u[U_SX ];
  upt[U_SY ] = sdetg*u[U_SY ];
  upt[U_SZ ] = sdetg*u[U_SZ ];
  upt[U_TAU] = sdetg*u[U_TAU];

  return errors == 0;
}
// }}}
