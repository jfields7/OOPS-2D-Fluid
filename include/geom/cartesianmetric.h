#ifndef CARTESIAN_METRIC_H
#define CARTESIAN_METRIC_H

#include "metric.h"
#include <cmath>


/*******************************************************************************
 *
 * Class: CartesianMetric
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 *
 * Description: A derived class of Metric that handles Cartesian metrics.
 *              Coordinate conventions are (x,y,z).
 *
 ******************************************************************************/

class CartesianMetric : public Metric{
  public:
    CartesianMetric(){
      gd[0][0] = 1.0;
      gd[0][1] = 0.0;
      gd[0][2] = 0.0;
      gd[1][0] = 0.0;
      gd[1][1] = 1.0;
      gd[1][2] = 0.0;
      gd[2][0] = 0.0;
      gd[2][1] = 0.0;
      gd[2][2] = 1.0;

      gu[0][0] = 1.0;
      gu[0][1] = 0.0;
      gu[0][2] = 0.0;
      gu[1][0] = 0.0;
      gu[1][1] = 1.0;
      gu[1][2] = 0.0;
      gu[2][0] = 0.0;
      gu[2][1] = 0.0;
      gu[2][2] = 1.0;

      det = 1.0;
      sdet = 1.0;
      for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
          for(unsigned int k = 0; k < 3; k++){
            christoffel[i][j][k] = 0;
          }
        }
      }
    }
    virtual ~CartesianMetric(){}

    // There's no need for the Cartesian metric ever to update.
    virtual void updateMetric(const double pos[2]){
      return;
    }
    virtual void updateChristoffelSymbols(const double pos[2]){
      return;
    }

    virtual double getLength(const double pos[2]){
      return sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    }
};


#endif
