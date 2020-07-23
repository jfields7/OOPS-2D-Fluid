#ifndef METRIC_H
#define METRIC_H
#include <domain.h>

/*******************************************************************************
 *
 * Class: Metric
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 *
 * Description: An abstract class for handling basic operations with a metric.
 *
 ******************************************************************************/

class Metric{
  protected:
    double **gd;
    double **gu;
    double ***christoffel;
    const double epsilon;
    double det;
    double sdet;
  public:
    // Metric {{{
    Metric() : epsilon(1e-15){
      gd = new double*[3];
      gu = new double*[3];
      christoffel = new double**[3];
      for(unsigned int i = 0; i < 3; i++){
        gd[i] = new double[3];
        gu[i] = new double[3];
        christoffel[i] = new double*[3];
        for(unsigned int j = 0; j < 3; j++){
          christoffel[i][j] = new double[3];
        }
      }
    }
    // }}}

    // ~Metric {{{
    virtual ~Metric(){
      for(unsigned int i = 0; i < 3; i++){
        delete[] gd[i];
        delete[] gu[i];
        for(unsigned int j = 0; j < 3; j++){
          delete[] christoffel[i][j];
        }
        delete[] christoffel[i];
      }

      delete[] christoffel;
      delete[] gd;
      delete[] gu;
    }
    // }}}
    
    /**
     * Update the metric variables based on the current position.
     */
    virtual void updateMetric(const double pos[2])=0;
    
    /**
     * Update the Christoffel symbols based on the current position.
     */
    virtual void updateChristoffelSymbols(const double pos[2])=0;

    /**
     * Get the length of a position vector measured relative to
     * the origin.
     */
    virtual double getLength(const double pos[2])=0;

    inline double** getMetric() const{
      return gd;
    }

    inline double** getInverseMetric() const{
      return gu;
    }

    inline double*** getChristoffelSymbols() const{
      return christoffel;
    }

    inline double getDeterminant() const{
      return det;
    }

    inline double getDeterminantRoot() const{
      return sdet;
    }

    // squareVector {{{
    inline double squareVector(double vu[3]){
      return gd[0][0]*vu[0]*vu[0] + 2.0*gd[0][1]*vu[0]*vu[1] +
         2.0*gd[0][2]*vu[0]*vu[0] +     gd[1][1]*vu[1]*vu[1] +
         2.0*gd[1][2]*vu[1]*vu[2] +     gd[2][2]*vu[2]*vu[2];
    }
    // }}}

    // squareForm {{{
    inline double squareForm(double vd[3]){
      return gu[0][0]*vd[0]*vd[0] + 2.0*gu[0][1]*vd[0]*vd[1] +
         2.0*gu[0][2]*vd[0]*vd[0] +     gu[1][1]*vd[1]*vd[1] +
         2.0*gu[1][2]*vd[1]*vd[2] +     gu[2][2]*vd[2]*vd[2];
    }
    // }}}

    // lowerVector {{{
    inline void lowerVector(double vd[3], double vu[3]){
      vd[0] = gd[0][0]*vu[0] + gd[0][1]*vu[1] + gd[0][2]*vu[2];
      vd[1] = gd[0][1]*vu[0] + gd[1][1]*vu[1] + gd[1][2]*vu[2];
      vd[2] = gd[0][2]*vu[0] + gd[1][2]*vu[1] + gd[2][2]*vu[2];
    }
    // }}}

    // raiseForm {{{
    inline void raiseForm(double vu[3], double vd[3]){
      vu[0] = gu[0][0]*vd[0] + gu[0][1]*vd[1] + gu[0][2]*vd[2];
      vu[1] = gu[0][1]*vd[0] + gu[1][1]*vd[1] + gu[1][2]*vd[2];
      vu[2] = gu[0][2]*vd[0] + gu[1][2]*vd[1] + gu[2][2]*vd[2];
    }
    // }}}

    /**
     * Transform a position to Cartesian coordinates.
     */
    virtual void toCartesianCoordinates(double out[2], const double in[2])=0;
    
    /**
     * Transform a position from Cartesian coordinates.
     */
    virtual void fromCartesianCoordinates(double out[2], const double in[2])=0;
};

#endif
