#include <riemann/riemannsolver.h>
#include <geom/metric.h>

// Constructor {{{
RiemannSolver::RiemannSolver(PrimitiveSolver* ps, Metric* m){
  primitive = ps;
  metric = m;
}
// }}}

// Destructor {{{
RiemannSolver::~RiemannSolver(){
  
}
// }}}

// calcFluxX {{{
void RiemannSolver::calcFluxX(double *f, double *u, double *v){
  const double Alpha = 1.0;
  const double Beta[3] = {0.0};

  double sdetg = metric->getDeterminantRoot();

  if(sdetg < 1.0e-7){
    f[U_D  ] = 0.0;
    f[U_SX ] = 0.0;
    f[U_SY ] = 0.0;
    f[U_SZ ] = 0.0;
    f[U_TAU] = 0.0;
  }

  double D, Tau, P;
  double Sd[3], Su[3], vd[3], vu[3];
  double isdetg = 1.0/(sdetg + 1.0e-15);

  D     = isdetg * u[U_D  ];
  Sd[0] = isdetg * u[U_SX ];
  Sd[1] = isdetg * u[U_SY ];
  Sd[2] = isdetg * u[U_SZ ];
  Tau   = isdetg * u[U_TAU];

  vu[0] = v[V_VX];
  vu[1] = v[V_VY];
  vu[2] = v[V_VZ];
  P     = v[V_P ];

  metric->lowerVector(vd, vu);
  metric->raiseForm(Su, Sd);

  f[U_D  ] = Alpha*sdetg*D*(vu[0] - Beta[0]/Alpha);
  f[U_SX ] = Alpha*sdetg*(Sd[0]*vu[0] + P - Beta[0]/Alpha*Sd[0]);
  f[U_SY ] = Alpha*sdetg*(Sd[1]*vu[0] - Beta[0]/Alpha*Sd[1]);
  f[U_SZ ] = Alpha*sdetg*(Sd[2]*vu[0] - Beta[0]/Alpha*Sd[2]);
  f[U_TAU] = Alpha*sdetg*(Su[0] - D*vu[0]) - sdetg*Beta[0]*Tau;

}
// }}}

// calcFluxY {{{
void RiemannSolver::calcFluxY(double *f, double *u, double *v){
  const double Alpha = 1.0;
  const double Beta[3] = {0.0};

  double sdetg = metric->getDeterminantRoot();

  if(sdetg < 1.0e-7){
    f[U_D  ] = 0.0;
    f[U_SX ] = 0.0;
    f[U_SY ] = 0.0;
    f[U_SZ ] = 0.0;
    f[U_TAU] = 0.0;
  }

  double D, Tau, P;
  double Sd[3], Su[3], vd[3], vu[3];
  double isdetg = 1.0/(sdetg + 1.0e-15);

  D     = isdetg * u[U_D  ];
  Sd[0] = isdetg * u[U_SX ];
  Sd[1] = isdetg * u[U_SY ];
  Sd[2] = isdetg * u[U_SZ ];
  Tau   = isdetg * u[U_TAU];

  vu[0] = v[V_VX];
  vu[1] = v[V_VY];
  vu[2] = v[V_VZ];
  P     = v[V_P ];

  metric->lowerVector(vd, vu);
  metric->raiseForm(Su, Sd);

  f[U_D  ] = Alpha*sdetg*D*(vu[1] - Beta[1]/Alpha);
  f[U_SX ] = Alpha*sdetg*(Sd[0]*vu[1] - Beta[1]/Alpha*Sd[0]);
  f[U_SY ] = Alpha*sdetg*(Sd[1]*vu[1] + P - Beta[1]/Alpha*Sd[1]);
  f[U_SZ ] = Alpha*sdetg*(Sd[2]*vu[1] - Beta[1]/Alpha*Sd[2]);
  f[U_TAU] = Alpha*sdetg*(Su[1] - D*vu[1]) - sdetg*Beta[1]*Tau;

}
// }}}
