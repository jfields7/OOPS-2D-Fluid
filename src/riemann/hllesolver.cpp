#include <riemann/hllesolver.h>
#include <cmath>
#include <geom/metric.h>
#include <primitive/primitivesolver.h>

// Constructor {{{
HLLESolver::HLLESolver(PrimitiveSolver* primitive, Metric* m) : RiemannSolver(primitive, m){}
// }}}

// Destructor {{{
HLLESolver::~HLLESolver(){

}
// }}}

// calcDefaultSpeeds {{{
void HLLESolver::calcDefaultSpeeds(double *bp, double *bm, const int dir){
  double **gu = metric->getInverseMetric();

  const double Alpha = 1.0;
  const double Beta[3] = {0.0};

  *bp = Alpha * sqrt(gu[dir][dir] - Beta[dir]);
  *bm = -Alpha * sqrt(gu[dir][dir] - Beta[dir]);
}
// }}}

// calcMaxSpeeds {{{
bool HLLESolver::calcMaxSpeeds(double *bp, double *bm, double *vl, double *vr, const int dir){
  const double Alpha = 1.0;
  const double Beta[3] = {0.0};

  double vlu[3], vru[3];

  double rhol = vl[V_RHO];
  double pl   = vl[V_P  ];
  vlu[0] = vl[V_VX];
  vlu[1] = vl[V_VY];
  vlu[2] = vl[V_VZ];

  double rhor = vr[V_RHO];
  double pr   = vr[V_P  ];
  vru[0] = vr[V_VX];
  vru[1] = vr[V_VY];
  vru[2] = vr[V_VZ];

  double csl, csr, cslsq, csrsq;
  csl = primitive->calcSoundSpeed(rhol, pl);
  csr = primitive->calcSoundSpeed(rhor, pr);
  cslsq = csl*csl;
  csrsq = csr*csr;

  double **gu = metric->getInverseMetric();
  double hkk = gu[dir][dir];
  double betak = Beta[dir];
  double vlk = vlu[dir];
  double vrk = vru[dir];

  double vlsq = metric->squareVector(vlu);
  double vrsq = metric->squareVector(vru);

  double dis;
  dis = cslsq*(1.0-vlsq)*(hkk*(1.0-vlsq*cslsq) - vlk*vlk*(1.0-cslsq));
  if (dis < 0.0){
    return false;
  }

  double Lpl = Alpha*(((1.0 - cslsq)*vlk+sqrt(dis))/(1.0 - vlsq*cslsq)) - betak;
  double Lml = Alpha*(((1.0 - cslsq)*vlk-sqrt(dis))/(1.0 - vlsq*cslsq)) - betak;

  dis = csrsq*(1.0-vrsq)*(hkk*(1.0-vrsq*csrsq) - vrk*vrk*(1.0-csrsq));
  if (dis < 0.0){
    return false;
  }
  
  double Lpr = Alpha*(((1.0 - csrsq)*vrk+sqrt(dis))/(1.0 - vrsq*csrsq)) - betak;
  double Lmr = Alpha*(((1.0 - csrsq)*vrk-sqrt(dis))/(1.0 - vrsq*csrsq)) - betak;

  double t1 = fmax(0.0, Lpl);
  *bp = fmax(t1, Lpr);

  t1 = fmin(0.0, Lml);
  *bm = fmin(t1, Lmr);

  return true;
}
// }}}

// solveRiemannProblem {{{
void HLLESolver::solveRiemannProblem(double *F, double *fl, double *fr,
                                     double *ul, double *ur, double *vl, double *vr,
                                     const int dir){
  double bp, bm;

  bool result = calcMaxSpeeds(&bp, &bm, vl, vr, dir);

  if(!result){
    calcDefaultSpeeds(&bp, &bm, dir);
  }

  //double rec = 1.0/(bp - bm);
  for(unsigned int m = 0; m < NU; m++){
    F[m] = ((bp*fl[m] - bm*fr[m]) + (ur[m] - ul[m])*bp*bm)/(bp - bm);
  }
}
// }}}
