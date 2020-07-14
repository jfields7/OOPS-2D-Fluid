#ifndef POLAR_METRIC_H
#define POLAR_METRIC_H

#include "metric.h"
#include <cmath>

/*******************************************************************************
 *
 * Class: PolarMetric
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 *
 * Description: A derived class of Metric that handles a 2d polar metric.
 *              This is best used with polar coordinates.
 *              Coordinate conventions are (r, theta, z);
 *
 ******************************************************************************/

class PolarMetric : public Metric{
  public:
    PolarMetric(const double pos[2]){
      gd[0][0] = 1.0;
      gd[0][1] = 0.0;
      gd[0][2] = 0.0;
      gd[1][0] = 0.0;
      gd[1][1] = pos[0]*pos[0];
      gd[1][2] = 0.0;
      gd[2][0] = 0.0;
      gd[2][1] = 0.0;
      gd[2][2] = 1.0;

      gu[0][0] = 1.0;
      gu[0][1] = 0.0;
      gu[0][2] = 0.0;
      gu[1][0] = 0.0;
      gu[1][1] = 1.0/(pos[0]*pos[0] + epsilon);
      gu[1][2] = 0.0;
      gu[2][0] = 0.0;
      gu[2][1] = 0.0;
      gu[2][2] = 1.0;

      det = gd[1][1];
      sdet = sqrt(det);
      for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
          for(unsigned int k = 0; k < 3; k++){
            christoffel[i][j][k] = 0;
          }
        }
      }
      christoffel[0][1][1] = -pos[0];
      christoffel[1][0][1] = 1.0/(pos[0] + epsilon);
      christoffel[1][1][0] = christoffel[1][0][1];
    }

    virtual ~PolarMetric(){}

    virtual void updateMetric(const double pos[2]){
      gd[1][1] = pos[0]*pos[0];
      gu[1][1] = 1.0/(pos[0]*pos[0] + epsilon);

      det = gd[1][1];
      sdet = sqrt(det);
    }

    virtual void updateChristoffelSymbols(const double pos[2]){
      christoffel[0][1][1] = -pos[0];
      christoffel[1][0][1] = 1.0/(pos[0] + epsilon);
      christoffel[1][1][0] = christoffel[1][0][1];
    }

    virtual double getLength(const double pos[2]){
      return pos[0];
    }
};

 #endif
