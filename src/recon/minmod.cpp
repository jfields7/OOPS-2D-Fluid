#include <recon/minmod.h>

Minmod::Minmod() : Recon(3){

}

Minmod::~Minmod(){
  
}

void Minmod::reconstruct(const unsigned int n, const double* const RESTRICT u,
                         double* const RESTRICT ul, double* const RESTRICT ur){
  ur[0] = u[0];
  ul[0] = u[0];

  for(unsigned int i = 1; i < n-1; i++){
    double df1 = u[i] - u[i-1];
    double df2 = u[i+1] - u[i];

    double slp = 0.5*minmod(df1, df2);
    ul[i+1] = u[i] + slp;
    ur[i] = u[i] - slp;
  }
  ul[n-1] = u[n-1];
  ur[n-1] = u[n-1];
}
