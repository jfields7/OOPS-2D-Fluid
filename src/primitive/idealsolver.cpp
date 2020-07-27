#include <primitive/idealsolver.h>
#include <cmath>
#include <geom/metric.h>
#include <iostream>
#include <mpicommunicator.h>

// Constructor {{{
IdealSolver::IdealSolver(Metric *m) : PrimitiveSolver(m){
  tol = 1e-14;
  maxIterations = 50;
  gamma = 5.0/3.0;
  #ifdef DEBUG_PRIMITIVE_SOLVER
  calls = 0;
  totalIterations = 0;
  #endif
}
// }}}

// Destructor {{{
IdealSolver::~IdealSolver(){}
// }}}

// calcH {{{
double IdealSolver::calcH(double rho, double P){
  return 1.0 + gamma/(gamma - 1.0)*P/rho;
}
// }}}

// calcSoundSpeed {{{
double IdealSolver::calcSoundSpeed(double rho, double p){
  const double factor = 1.0e3;
  double y, cs;

  if(rho*factor < p){
    y = (gamma - 1.0*rho/(gamma*p));
    cs = sqrt(gamma - 1.0)*(1.0 + (-0.5+(0.375+(-0.3125+0.2734375*y)*y)*y)*y);
  }
  else{
    cs = sqrt(gamma*(gamma - 1.0)*p / ((gamma - 1.0)*rho + gamma*p));
  }

  return cs;
}
// }}}

// conToPrimPt {{{
bool IdealSolver::conToPrimPt(double *u, double *v){
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
    upt[U_TAU] = fmax(vacuum_tau,upt[U_TAU]);
    Ssq = 0.0;
    errors++;
  }
  else{
    if(upt[U_TAU] < vacuum_tau){
      upt[U_TAU] = vacuum_tau;
      errors++;
    }
    double Ssq_max = (2.0*upt[U_D] + upt[U_TAU])*upt[U_TAU];
    if(Ssq_max < Ssq){
      double t_min = VELOCITY_FACTOR * sqrt(Ssq_max/Ssq);
      upt[U_SX] *= t_min;
      upt[U_SY] *= t_min;
      upt[U_SZ] *= t_min;
      Sd[0] = upt[U_SX];
      Sd[1] = upt[U_SY];
      Sd[2] = upt[U_SZ];
      Ssq *= t_min*t_min;
      errors++;
    }
  }

  double D = upt[U_D];
  double tau = upt[U_TAU];

  unsigned int j;
  double qf_1, qf_2, qg_2, qf_max;
  double ql, qh, f_at_qh, f_at_qh_v2;
  double f_at_ql = 0.0;
  double rescale_factor1;
  double a = tau + D;
  double t1 = 0.5*gamma;
  double beta_sq = Ssq/(a*a);
  double t2 = t1*t1 - (gamma-1.0)*beta_sq;
  // t2 shows up under a square root, so we check to see if it is < 0.
  if (t2 < 0.0){
    std::cout << "primToConPt: Problem with upper bound. t2 = " << t2 << "\n";
    MPI_Abort(MPI_COMM_WORLD,0);
    errors++;
  }
  double sqrt_t2 = sqrt(t2);
  qf_1 = t1 - sqrt_t2;  // The smaller root of f
  qf_2 = t1 + sqrt_t2;  // The larger root of f
  qf_max = t1;          // The location of the max of f
  //qg_2 = sqrt(beta_sq); // The larger root of g
  qg_2 = 1.0;
  
  if(std::isnan(qg_2)){
    std::cout << "qg_2 is a nan.\n";
    errors++;
  }

  // If the velocity vector is effectively zero, then we have an analytic
  // transformation.
  if(beta_sq < 1.e-16){
    v[V_RHO] = D;
    v[V_P  ] = (gamma - 1.0)*tau;
    v[V_VX ] = 0.0;
    v[V_VY ] = 0.0;
    v[V_VZ ] = 0.0;
  }
  else{
    if(fabs((qg_2 - qf_1)/(qg_2 + qf_1 + 1.0e-15)) < 1.e-6){
      qg_2 += 1.e-5;
    }
    // Set the lower bound of qg_2 unless qg_2 is smaller than the location
    // of the maximum of f.
    ql = fmax(qg_2, qf_max);

    // This temporarily initializes the upper bound to a possible value.
    // The following code checks this to make sure it works.
    qh = qf_2;

    // If the inequality is obeyed, then we can move on.
    if(qg_2 < qf_2 && beta_sq < 1.0){
      qh = qf_2;
    }
    // An inequality is violated, so we rescale S.
    else if (beta_sq + (D/a)*(D/a) - 1.0 >= 0.0){
      errors++;
      //FIXME: This scaling factor is arbitrary. It would be well
      //       to experiment with the value added in.
      rescale_factor1 = 1.0/(sqrt(beta_sq + (D/a)*(D/a)) + 1.e-14);

      // Multiply Sd by the scaling factor.
      for(j = 0; j < 3; j++){
        Sd[j] *= rescale_factor1;
      }

      // Recalculate Ssq.
      Ssq = metric->squareForm(Sd);

      // Recalculate the roots of f and g with this new Ssq.
      beta_sq = Ssq/(a*a);
      t2 = t1*t1 - (gamma-1.0)*beta_sq;
      if(t2 < 0.0){
        std::cout << "primToConPt: t2 < 0, i.e. t2=" << t2 << std::endl;
      }
      qf_2 = t1 + sqrt(t2); // new root of f; upper bound on q
      qg_2 = sqrt(beta_sq); // new root of g; lower bound on q

      // Check that these new roots are indeed good upper and lower bounds.
      if (qg_2 < qf_2){
        ql = qg_2;
        qh = qf_2;
      }

      // Because we made a rescaling of the components of S and of Ssq, we need
      // to save these.
      upt[U_SX] = Sd[0];
      upt[U_SY] = Sd[1];
      upt[U_SZ] = Sd[2];
    }
    // If qg_2 and qf_2 are exactly the same, push the roots apart.
    else if(qg_2 == qf_2){
      qg_2 -= 1.e-16;
      qf_2 += 1.e-16;

      ql = qg_2;
      qh = qf_2;
    }
    else{
      // We don't know what happened. Scream and quit.
      printf("primToConPt:  Something unexpected (1) \n");
      printf("  qg_2 = %25.20e \n", qg_2);
      printf("  qf_2 = %25.20e \n", qf_2);
      printf("  beta_sq = %25.20e \n", beta_sq);
      printf("  D = %25.20e \n", D);
      printf("  tau = %25.20e \n", tau);
      printf("  Ssq = %25.20e \n", Ssq);
      printf("  beta_sq + (D/a)^2 = %25.20e \n", beta_sq+(D/a)*(D/a));
      MPI_Abort(MPI_COMM_WORLD, 0);
      errors++;
    }

    // At this point, we should have the two values, namely qg_2 and qf_2,
    // that should form a legitimate bracket for finding the root of F(q).
    // Let's check.
    //
    // We want to check that at ql, the function F(q) := f(q)-(const)*sqrt(g(q))
    // is positive. First, just calculate F(ql).
    f_at_ql = f(ql, D, tau, Ssq);

    // Now check that at qh, the function F(q) := f(q)-(const)*sqrt(g(q))
    // is negative. For now, just calculate F(qh).
    f_at_qh = f(qh, D, tau, Ssq);

    // If, indeed, F(ql) > 0 and F(qh) < 0, then we should have a bracket
    // for the root. We will check for the failure of this and do
    // what we can.
    // FIXME: Actually, if there's a problem, we don't do anything here.
    if(f_at_ql < 0.0 || f_at_qh > 0.0){
    
    }

    if(f_at_ql > 0.0 && f_at_qh > 0.0){
      qh += 6.e-16;

      f_at_qh = f(qh, D, tau, Ssq);
    }

    // At this point, we should have a legitimate bracket for our root, namely
    // [ql, qh] := [qg_2, qf_2]. We should now be able to proceed with the
    // Newton solve.

    // As our initial guess for the root, take the midpoint of the bracket.
    //double qi = 0.5*(ql + qh);
    // We have the function values already, so we estimate the root using
    // false position instead.
    double qi = (ql*f_at_qh - qh*f_at_ql)/(f_at_qh - f_at_ql);
    double q = qi;

    // Incorporate a little redundancy in the event things go bad.
    double ql0 = ql;
    double qh0 = qh;
    double q0 = q;

    // Perform the Newton solve.
    double fq = 1.0;
    double dfq = 1.0;
    double ddfq = 1.0;
    unsigned int count = 0;
    int held = 0;
    while(fabs(fq) > tol && count < maxIterations){
      double oldq = q;
      //fd(fq, dfq, q, D, tau, Ssq);
      //q = q - fq/dfq;
      // Halley's method
      fdd(fq, dfq, ddfq, q, D, tau, Ssq);
      q = q - 2.0*fq*dfq/(2.0*dfq*dfq - fq*ddfq);
      //q = q - 2.0*fq/dfq
      // Check for convergence by ensuring that the new q is
      // still inside the bounds. If not, apply the Illinois
      // variant of the false position method.
      if(q > qh || q < ql){
        // Calculate the new root.
        if(fq*f_at_qh > 0){
          double m;
          if(held == -1){
            m = 0.5;
          }
          else{
            m = 1.0;
          }
          //double m = 1.0;
          q = (f_at_qh*ql - m*f_at_ql*qh)/(f_at_qh - m*f_at_ql);
        }
        else{
          double m;
          if(held == 1){
            m = 0.5;
          }
          else{
            m = 1.0;
          }
          //double m = 1.0;
          q = (m*f_at_qh*ql - f_at_ql*qh)/(m*f_at_qh - f_at_ql);
        }
        //q = 0.5*(qh + ql);
      }
      // Move the bounds according to the sign of fq.
      if(fq < 0){
        held = -1;
        qh = oldq;
        f_at_qh = fq;
      }
      else{
        held = 1;
        ql = oldq;
        f_at_ql = fq;
      }
      count++;
    }
    // In the event of failure, specify what went on.
    if(fabs(fq) > tol){
      printf("primToConPt: The Newton solve failed.\n");
      printf("primToConPt: Current values include \n");
      printf("  ql = %25.20e \n",ql);
      printf("  qh = %25.20e \n",qh);
      printf("  qi = %25.20e \n",qi);
      printf("  q  = %25.20e \n",q);
      printf("  ql0 = %25.20e \n", ql0);
      printf("  qh0 = %25.20e \n", qh0);
      printf("  q0  = %25.20e \n", q0);
      printf("primToConPt: Some conserved quantities and \n");
      printf("primToConPt: inequalities, \n");
      printf("  D               = %25.20e\n", D);
      printf("  tau+D           = %25.20e \n",a);
      printf("  D^2/a^2a        = %25.20e \n", (D*D)/(a*a));
      printf("  Ssq/(a*a)       = %25.20e \n", Ssq/(a*a));
      printf("  (Ssq+D^2)/(a*a) = %25.20e \n", (Ssq+D*D)/(a*a));

      MPI_Abort(MPI_COMM_WORLD,0);
      errors++;
    }
    #ifdef DEBUG_PRIMITIVE_SOLVER
    totalIterations += count;
    calls++;
    #endif

    // If we get to this point, we may actually have a solution.
    // We have to rescale back by a = tau + D.
    q *= a;

    // Knowing the value for q (:=hW^2, i.e., related to the enthalpy), we
    // can now get the remaining primitive variables.  Start with the velocity
    // squared: v^2.
    double vsq = Ssq / (q*q);
    #ifdef DEBUG_PRIMITIVE_SOLVER
    if(vsq >= 1.0){
      printf("primToConPt: Unphysical velocity.\n");
      printf("  vsq = %g\n",vsq);
    }
    #endif
    // We have the "down" components of S; now get the "up" components.
    double Su[3];
    metric->raiseForm(Su, Sd);

    // Now calculate rho, pressure, and the "up" components of v^i.
    v[V_RHO] = D*sqrt(1.0 - vsq);
    v[V_P  ] = (gamma-1.0)/gamma*(q*(1.0-vsq) - v[V_RHO]);
    v[V_VX ] = Su[0] / q;
    v[V_VY ] = Su[1] / q;
    v[V_VZ ] = Su[2] / q;

    #ifdef DEBUG_PRIMITIVE_SOLVER
    double vu[3] = {v[V_VX], v[V_VY], v[V_VZ]};
    if(metric->squareVector(vu) >= 1.0){
      printf("primToConPt: Unphysical velocity.\n");
      printf("  vx = %g\n", vu[0]);
      printf("  vy = %g\n", vu[1]);
      printf("  vz = %g\n", vu[2]);
    }
    #endif

    // In some cases, if the velocity is very small, P can be
    // identically zero. Therefore, we'll define it analytically.
    if(v[V_P] == 0){
      v[V_P] = (gamma - 1.0)*tau;
    }
  }
  // Apply the floor and move on.
  if(v[V_RHO] < vacuum){
    v[V_RHO] = vacuum;
    v[V_VX ] = 0.0;
    v[V_VY ] = 0.0;
    v[V_VZ ] = 0.0;
    upt[U_D  ] = vacuum;
    upt[U_SX ] = 0.0;
    upt[U_SY ] = 0.0;
    upt[U_SZ ] = 0.0;
    errors++;
  }
  if(v[V_P] < vacuum_tau){
    v[V_P] = vacuum;
    errors++;
  }

  // Redensitize the variables.
  u[U_D  ] = sdetg*upt[U_D  ];
  u[U_SX ] = sdetg*upt[U_SX ];
  u[U_SY ] = sdetg*upt[U_SY ];
  u[U_SZ ] = sdetg*upt[U_SZ ];
  u[U_TAU] = sdetg*upt[U_TAU];

  // To keep things consistent, apply primToConPt and move on.
  //primToConPt(u, v);
  return errors == 0;
}
// }}}

// f {{{
double IdealSolver::f(double x, double D, double tau, double Ssq){
  //double betainv = (gamma - 1.0)/gamma;
  //return x*(x*(1.0 - betainv) + betainv) + asq - b*sqrt(x*x + asq);
  double a = tau + D;
  double betasq = Ssq/(a*a);
  double dis = x*x - betasq;
  if( dis < 1.e-18){
    dis = 1.e-18;
  }
  double xoff = x - 0.5*gamma;

  return -xoff*xoff + gamma*gamma/4.0
     - (gamma-1.0)*betasq - (gamma-1.0)*D/a*sqrt(dis);
}
// }}}

// fd {{{
void IdealSolver::fd(double& fx, double& dfx, double x, double D, double tau, double Ssq){
  //double betainv = (gamma - 1.0)/gamma;
  //return 2.0*x*(1.0 - betainv) + betainv - b*x/sqrt(x*x+asq);
  double a = tau + D;
  double betasq = Ssq/(a*a);
  double dis = x*x - betasq;
  if( dis < 1.e-18){
    dis = 1.e-18;
  }
  double sdis = sqrt(dis);
  double xoff = x - 0.5*gamma;
  double coeff = (gamma-1.0)*D/a;

  fx = -xoff*xoff + gamma*gamma/4.0 - (gamma-1.0)*betasq - coeff*sdis;
  dfx = -2.0*xoff - coeff*x/sdis;
}
// }}}

// fdd {{{
void IdealSolver::fdd(double& fx, double& dfx, double& ddfx, double x, double D, double tau, double Ssq){
  double a = tau + D;
  double betasq = Ssq/(a*a);
  double dis = x*x - betasq;
  if(dis < 1.e-18){
    dis = 1.e-18;
  }
  double sdis = sqrt(dis);
  double xoff = x - 0.5*gamma;
  double coeff = (gamma-1.0)*D/a;

  fx = -xoff*xoff + gamma*gamma/4.0 - (gamma-1.0)*betasq - coeff*sdis;
  dfx = -2.0*xoff - coeff*x/sdis;
  ddfx = -2.0 + coeff*betasq/(dis*sdis);
}
// }}}
