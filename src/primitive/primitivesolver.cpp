#include <primitive/primitivesolver.h>
#include <iostream>
#include <geom/metric.h>
#include <cmath>

// Constructor {{{
PrimitiveSolver::PrimitiveSolver(Metric* m){
  metric = m;
  vacuum = 1e-10;
  vacuum_tau = 1e-15;
}
// }}}

// Destructor {{{
PrimitiveSolver::~PrimitiveSolver(){}
// }}}

// calcW {{{
double PrimitiveSolver::calcW(double vu[3]){
  double vsq = metric->squareVector(vu);
  return 1.0/sqrt(1.0 - vsq);
}
// }}}

// primToConPt {{{
void PrimitiveSolver::primToConPt(double *u, double *v){
  double rho = v[V_RHO];
  double P   = v[V_P];

  if(rho < vacuum){
    std::cout << "primToConPt: rho < vacuum. rho=" << rho << std::endl;
    rho = vacuum;
  }
  if(P < vacuum_tau){
    std::cout << "primToConPt: P < vacuum. P=" << rho << std::endl;
    P = vacuum_tau;
  }

  double vu[3], vd[3];
  vu[0] = v[V_VX];
  vu[1] = v[V_VY];
  vu[2] = v[V_VZ];

  double vsq = metric->squareVector(vu);
  metric->lowerVector(vd, vu);

  double Wsq = 1.0/(1.0 - vsq);
  double h = rho*calcH(rho, P);
  double hWsq = h*Wsq;

  double sdetg = metric->getDeterminantRoot();

  u[U_D  ] = sdetg*(rho*sqrt(Wsq));
  u[U_TAU] = sdetg*(hWsq - P) - u[U_D];
  u[U_SX ] = sdetg*(hWsq*vd[0]);
  u[U_SY ] = sdetg*(hWsq*vd[1]);
  u[U_SZ ] = sdetg*(hWsq*vd[2]);
}
// }}}
