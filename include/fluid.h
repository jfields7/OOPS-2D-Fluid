#ifndef FLUID_H
#define FLUID_H
#include <ode.h>

class Metric;
class RiemannSolver;
class PrimitiveSolver;
class Recon;
class FluidICParameters;
class FluidGridParameters;

class Fluid : public ODE{
  private:
    // Variable labels
    // Conserved variables
    static const unsigned int NU    = 5;
    static const unsigned int U_D   = 0;
    static const unsigned int U_SX  = 1;
    static const unsigned int U_SY  = 2;
    static const unsigned int U_SZ  = 3;
    static const unsigned int U_TAU = 4;
    // Primitive variables
    static const unsigned int NV    = 5;
    static const unsigned int V_RHO = 0;
    static const unsigned int V_VX  = 1;
    static const unsigned int V_VY  = 2;
    static const unsigned int V_VZ  = 3;
    static const unsigned int V_P   = 4;
    // Directions
    static const unsigned int DIR_X = 0;
    static const unsigned int DIR_Y = 1;

    // Metric
    Metric* metric;

    // Various parts of the fluid method.
    RiemannSolver *riemann;
    PrimitiveSolver *primitive;
    Recon* recon;

    void reconstruct(double *vl[], double *vr[], double *v1[], double *wv[],
                     const unsigned int dir, const Grid& grid, const unsigned int j);

    void calculateFourVelocity(double *wv[], double *v1[], const unsigned int dir, const Grid& grid,
                               const unsigned int j);

    void calculateSourceTerms(double *dtu[], double *u[], double *v[], const Grid& grid);
    void interpolateBadPoints(std::set<unsigned int>& bad);
    double averageNeighbors(double *vpt, double *v[], std::vector<unsigned int>& neighbors, unsigned int ix,
                            unsigned int jy);

    // Parameters
    FluidICParameters *params;
    FluidGridParameters *gridParams;

    // Initial conditions
    void applyGaussian();
    void applyShockTube1D(const unsigned int dir);

  protected:
    virtual void rhs(std::shared_ptr<FieldMap>& fieldMap);
    virtual void doAfterStage();

  public:
    Fluid(Domain *d, Solver *s);
    virtual ~Fluid();

    virtual void initData();
    virtual void applyBoundaries();
    virtual void doAfterBoundaries();

    // Some inline getters and setters. {{{
    inline Metric* getMetric(){
      return metric;
    }
    inline void setMetric(Metric* m){
      metric = m;
    }

    inline RiemannSolver* getRiemannSolver(){
      return riemann;
    }
    inline void setRiemannSolver(RiemannSolver *rs){
      riemann = rs;
    }

    inline PrimitiveSolver* getPrimitiveSolver(){
      return primitive;
    }
    inline void setPrimitiveSolver(PrimitiveSolver *ps){
      primitive = ps;
    }

    inline Recon* getRecon(){
      return recon;
    }
    inline void setRecon(Recon *r){
      recon = r;
    }

    inline FluidICParameters* getParameters(){
      return params;
    }
    inline void setParameters(FluidICParameters *par){
      params = par;
    }
    inline FluidGridParameters* getGridParameters(){
      return gridParams;
    }
    inline void setGridParameters(FluidGridParameters *par){
      gridParams = par;
    }
    // }}}

    static inline std::string charToString(char *str){
      return std::string(str);
    }
};

#endif
