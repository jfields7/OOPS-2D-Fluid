#include <eospiece.h> 
#include <cmath>


EOSPiece::EOSPiece(double r, double k, double g, double c): rho(r), kappa(k), gamma(g), a(c){
  pressure = kappa*std::pow(rho,gamma);
  enthalpy = rho*(1.0 + a) + gamma/(gamma - 1.0)*pressure;
}

double EOSPiece::getTrialSsq(double D, double tau){
  double x = tau + D + pressure;
  return x*(x + enthalpy);
}

double EOSPiece::calcH(double dens, double P){
  return dens*(1.0 + a) + gamma/(gamma - 1.0)*P;
}

double EOSPiece::calcSoundSpeed(double dens, double P){
  // This factor is an arbitrary threshold that we use
  // to define what "much less than" is.
  const double factor = 1.0e3;
  double y, cs;

  if(rho*factor < p){
    y = (gamma - 1.0)*(a + 1.0)*dens/(gamma*P);
    cs = sqrt(gamma - 1.0)*(1.0 + (-0.5+(0.375+(-0.3125+0.2734375*y)*y)*y)*y);
  }
  else{
    cs = sqrt(gamma*(gamma - 1.0)*P/((gamma - 1.0)*(a + 1.0)*dens + gamma*P));
  }

  return cs;
}
