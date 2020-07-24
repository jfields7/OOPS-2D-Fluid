#include <recon/norecon.h>

NoRecon::NoRecon() : Recon(2){
  
}

NoRecon::~NoRecon(){
  
}

void NoRecon::reconstruct(const unsigned int n, const double* const RESTRICT u,
                          double* const RESTRICT ul, double* const RESTRICT ur){
  ul[0] = u[0];
  for(unsigned int i = 0; i < n-1; i++){
    ul[i+1] = u[i];
    ur[i] = u[i];
  }
  ur[n-1] = u[n-1];
}

void NoRecon::reconstructPt(const unsigned int i, const unsigned int n,
                            const double* const RESTRICT u,
                            double* const RESTRICT ul, double* const RESTRICT ur){
  ur[i] = u[i];
  if(i == 0){
    ul[i] = u[i];
    ul[i+1] = u[i];
  }
  else{
    ul[i+1] = u[i];
  }
}
