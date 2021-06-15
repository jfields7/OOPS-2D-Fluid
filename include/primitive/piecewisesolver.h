#ifndef PIECEWISE_SOLVER_H
#define PIECEWISE_SOLVER_H

#include <primitive/primitivesolver.h>

/***************************************************************************
 *
 * Class: PiecewiseSolver
 * Author: Jacob Fields
 * Date Modified: 4-6-2021
 *
 * Description: A solver for an piecewise-polytropic gas.
 *
 **************************************************************************/

class EOSPiece;

class PiecewiseSolver : public PrimitiveSolver{
  protected:
    unsigned int nPieces;
    EOSPiece *pieces;
  
    double tol;
    unsigned int maxIterations;
  public:
    PiecewiseSolver(Metric* m, unsigned int n);
    virtual ~PiecewiseSolver();

    virtual double calcH(double rho, double P);
    
    virtual bool conToPrimPt(double *u, double *v);

    virtual double calcSoundSpeed(double rho, double p);

    void setPiece(unsigned int i, double rho, double kappa, double gamma, double a);

    inline EOSPiece& getPiece(unsigned int i){
      return piece[i];
    }
    
    inline unsigned int getNPieces() const{
      return nPieces;
    }

    inline unsigned int getMaxIterations(){
      return maxIterations;
    }

    inline void setMaxIterations(unsigned int i){
      maxIterations = i;
    }
};

#endif
