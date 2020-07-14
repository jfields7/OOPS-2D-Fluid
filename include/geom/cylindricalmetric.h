#ifndef CYLINDRICAL_METRIC_H
#define CYLINDRICAL_METRIC_H

#include "metric.h"
#include <cmath>

/*******************************************************************************
 *
 * Class: CylindricalMetric
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 *
 * Description: A derived class of Metric that handles a 3d cylindrical metric.
 *              This is best used with Cartesian coordinates.
 *              Coordinate conventions are (r, z, theta);
 *
 ******************************************************************************/

class CylindricalMetric : public Metric{
  public:
    CylindricalMetric(const double pos[2]){
      gd[0][0] = 1.0;
      gd[0][1] = 0.0;
      gd[0][2] = 0.0;
      gd[1][0] = 0.0;
      gd[1][1] = 1.0;
      gd[1][2] = 0.0;
      gd[2][0] = 0.0;
      gd[2][1] = 0.0;
      gd[2][2] = pos[0]*pos[0];

      gu[0][0] = 1.0;
      gu[0][1] = 0.0;
      gu[0][2] = 0.0;
      gu[1][0] = 0.0;
      gu[1][1] = 0.0;
      gu[1][2] = 0.0;
      gu[2][0] = 0.0;
      gu[2][1] = 0.0;
      gu[2][2] = 1.0/(pos[0]*pos[0] + epsilon);

      det = gd[2][2];
      sdet = sqrt(det);
      for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
          for(unsigned int k = 0; k < 3; k++){
            christoffel[i][j][k] = 0;
          }
        }
      }
      christoffel[0][2][2] = -pos[0];
      christoffel[2][0][2] = 1.0/(pos[0] + epsilon);
      christoffel[2][2][0] = christoffel[2][0][2];
    }

    virtual ~CylindricalMetric(){}

    virtual void updateMetric(const double pos[2]){
      gd[2][2] = pos[0]*pos[0];
      gu[2][2] = 1.0/(pos[0]*pos[0] + epsilon);

      det = gd[2][2];
      sdet = sqrt(det);
    }

    virtual void updateChristoffelSymbols(const double pos[2]){
      christoffel[0][2][2] = -pos[0];
      christoffel[2][0][2] = 1.0/(pos[0] + epsilon);
      christoffel[2][2][0] = christoffel[2][0][2];
    }

    virtual double getLength(const double pos[2]){
      return sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    }
};

 #endif
