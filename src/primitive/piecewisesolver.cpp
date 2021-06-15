#include <primitive/piecewisesolver.h>
#include <primitive/eospiece.h>

// Constructor {{{
PiecewiseSolver::PiecewiseSolver(Metric* m, unsigned int n) : PrimitiveSolver(m){
  nPieces = n;
  tol = 1e-14;
  maxIterations = 50;
#ifdef DEBUG_PRIMITIVE_SOLVER
  calls = 0;
  totalIterations = 0;
#endif

  pieces = new EOSPiece[nPieces];
}
// }}}

// Destructor {{{
PiecewiseSolver::~PiecewiseSolver(){
  delete[] pieces;
}
// }}}

// calcH {{{
double PiecewiseSolver::calcH(double rho, double p){
  // First, we need to identify which EOSPiece this agrees with.
  // We traverse the pieces in reverse order so we can find the
  // first piece that we match with.
  for(int i = nPieces-1; i >= 0; i--){
    if(pieces[i].getRho() < rho){
      // Once we've identified which piece we're using, just
      // calculate the enthalpy directly using that piece.
      return pieces[i].calcH(rho, P);
    }
  }
}
// }}}

double PiecewiseSolver::calcSoundSpeed(double rho, double p){
  // We'll basically do the same thing that we did for enthalpy:
  // find the piece that actually matters, then calculate it
  // directly from that piece.
  for(int i = nPieces-1; i >= 0; i--){
    if(pieces[i].getRho() < rho){
      return pieces[i].calcSoundSpeed(dens, p);
    }
  }
}

void PiecewiseSolver::setPiece(unsigned int i, double rho, double kappa, double gamma, double a){
  pieces[i] = EOSPiece(rho, kappa, gamma, a);
}

bool PiecewiseSolver::conToPrimPt(double *u, double *v){
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
    upt[U_D] = vacuum;
    upt[U_SX] = 0.0;
    upt[U_SY] = 0.0;
    upt[U_SZ] = 0.0;
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
  }
}
