#include <recon/mp5.h>

MP5::MP5() : Recon(5){
  
}

MP5::~MP5(){

}

void MP5::reconstruct(const unsigned int n, const double* const RESTRICT u,
                      double* const RESTRICT ul, double* const RESTRICT ur){
  #define A(i_) (u[ijk[i_]])
  #define Aplus(i_) (ul[ijk[i_]])
  #define Aminus(i_) (ur[ijk[i_]])

  const double mp5_alpha = 4.0;
  const double mp5_eps = 1.0e-10;

  for(unsigned int i = 2; i < n-2; i++){
    const unsigned int ijk[5] = {i-2, i-1, i, i+1, i+2};

    if(!adaptive){
      Aplus(3) = doMP5(A(0), A(1), A(2), A(3), A(4), 1.0, mp5_eps, mp5_alpha);
      Aminus(2) = doMP5(A(4), A(3), A(2), A(1), A(0), 1.0, mp5_eps, mp5_alpha);
    }
    else{
      const double anorm = sqrt(A(0)*A(0) + A(1)*A(1) + A(2)*A(2) +A(3)*A(3) + A(4)*A(4));

      Aplus(3) = doMP5(A(0), A(1), A(2), A(3), A(4), anorm, mp5_eps, mp5_alpha);
      Aminus(2) = doMP5(A(4), A(3), A(2), A(1), A(0), anorm, mp5_eps, mp5_alpha);
    }
  }

  ul[0] = u[0];
  ul[1] = u[0];
  ul[2] = u[1];
  ur[0] = u[0];
  ur[1] = u[1];

  ur[n-2] = u[n-2];
  ur[n-1] = u[n-1];
  ul[n-1] = u[n-2];

  #undef A
  #undef Aplus
  #undef Aminus
}

void MP5::reconstructPt(const unsigned int i, const unsigned int n,
                        const double* const RESTRICT u,
                        double* const RESTRICT ul, double* const RESTRICT ur){
  const double mp5_alpha = 4.0;
  const double mp5_eps = 1.0e-10;
  // FIXME: Not currently implemented.
  if(i == 0){
    ul[0] = u[0];
    ul[1] = u[0];
    ur[0] = u[0];
  }
  else if(i == 1){
    ul[2] = u[1];
    ur[1] = u[1];
  }
  else if(i == n-2){
    ur[n-2] = u[n-2];
  }
  else if(i == n-1){
    ur[n-1] = u[n-1];
  }
  else{
    if(!adaptive){
      ul[i+1] = doMP5(u[i-2], u[i-1], u[i], u[i+1], u[i+2], 1.0, mp5_eps, mp5_alpha);
      ur[i] = doMP5(u[i+2], u[i+1], u[i], u[i-1], u[i-2], 1.0, mp5_eps, mp5_alpha);
    }
    else{
      const double anorm = sqrt(u[i-2]*u[i-2] + u[i-1]*u[i-1] + u[i]*u[i] + u[i+1]*u[i+1] + u[i+2]*u[i+2]);
      ul[i+1] = doMP5(u[i-2], u[i-1], u[i], u[i+1], u[i+2], 1.0, mp5_eps, mp5_alpha);
      ur[i] = doMP5(u[i+2], u[i+1], u[i], u[i-1], u[i-2], 1.0, mp5_eps, mp5_alpha);
    }
  }
}
