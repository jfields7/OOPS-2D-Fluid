#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

/******************************************************************
 *
 * Class: RiemannSolver
 * Author: Jacob Fields
 * Date Modified: 16-6-2020
 *
 * Description: An abstract class for a Riemann solver.
 *
 *****************************************************************/

class Metric;
class PrimitiveSolver;

class RiemannSolver{
  protected:
    const unsigned int NU    = 5;
    const unsigned int U_D   = 0;
    const unsigned int U_SX  = 1;
    const unsigned int U_SY  = 2;
    const unsigned int U_SZ  = 3;
    const unsigned int U_TAU = 4;

    const unsigned int NV    = 5;
    const unsigned int V_RHO = 0;
    const unsigned int V_VX  = 1;
    const unsigned int V_VY  = 2;
    const unsigned int V_VZ  = 3;
    const unsigned int V_P   = 4;

    // We maintain a pointer to the EOS so that we can get things
    // like the sound speed.
    PrimitiveSolver *primitive;
    Metric* metric;

  public:
    RiemannSolver(PrimitiveSolver* ps, Metric* m);
    virtual ~RiemannSolver();
    
    /**
     * Calculate the flux at a point based on the input values.
     */
    void calcFluxX(double *f, double *u, double *v);
    void calcFluxY(double *f, double *u, double *v);

    /**
     * Solve the Riemann problem at a specific point.
     */
    virtual void solveRiemannProblem(double *F, double *fl, double *fr, 
                                     double *ul, double *ur, double *vl, double *vr,
                                     const int dir)=0;

    inline PrimitiveSolver* getPrimitiveSolver(){
      return primitive;
    }
    inline void setPrimitiveSolver(PrimitiveSolver* ps){
      primitive = ps;
    }

    inline Metric* getMetric(){
      return metric;
    }
    inline void setMetric(Metric* m){
      metric = m;
    }
};

#endif
