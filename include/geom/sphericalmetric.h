#ifndef SPHERICAL_METRIC_H
#define SPHERICAL_METRIC_H
#include <cmath>

/*******************************************************************************
 *
 * Class: SphericalMetric 
 * Author: Jacob Fields
 * Date Modified: 15-6-2020
 *
 * Description: A derived class of Metric that handles a 3d spherical metric.
 *              This is best used with polar coordinates.
 *              Coordinate conventions are (r, theta, phi).
 *              Note that theta is actually shifted relative to the rest of the
 *              code; while theta corresponds to the counterclockwise angle 
 *              from the x-axis for purposes of output and finite differencing,
 *              theta here is measured clockwise from the z-axis. The metric
 *              code handles the conversion, so (r, theta) actually treats the
 *              point like (r, pi/2 - theta).
 *
 ******************************************************************************/

class SphericalMetric : public Metric {
  public:
    SphericalMetric(const double pos[2]){
      double rsq = pos[0]*pos[0];
      double sth = cos(pos[1]);
      double cth = sin(pos[1]);
      gd[0][0] = 1.0;
      gd[0][1] = 0.0;
      gd[0][2] = 0.0;
      gd[1][0] = 0.0;
      gd[1][1] = rsq;
      gd[1][2] = 0.0;
      gd[2][0] = 0.0;
      gd[2][1] = 0.0;
      gd[2][2] = rsq*sth*sth;

      gu[0][0] = 1.0;
      gu[0][1] = 0.0;
      gu[0][2] = 0.0;
      gu[1][0] = 0.0;
      gu[1][1] = 1.0/(rsq + epsilon);
      gu[1][2] = 0.0;
      gu[2][0] = 0.0;
      gu[2][1] = 0.0;
      gu[2][2] = 1.0/(rsq*sth*sth + epsilon);

      det = rsq*rsq*sth*sth;
      sdet = sqrt(det);
      for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
          for(unsigned int k = 0; k < 3; k++){
            christoffel[i][j][k] = 0;
          }
        }
      }
      christoffel[0][1][1] = -pos[0];
      christoffel[0][2][2] = -pos[0]*sth*sth;
      christoffel[1][0][1] = 1.0/(pos[0] + epsilon);
      christoffel[1][1][0] = christoffel[1][0][1];
      christoffel[1][2][2] = sth*cth;
      christoffel[2][0][2] = christoffel[1][0][1];
      christoffel[2][2][0] = christoffel[1][0][1];
      christoffel[2][1][2] = cth/(sth + epsilon);
      christoffel[2][2][1] = christoffel[2][1][2];
    }

    virtual ~SphericalMetric(){}

    virtual void updateMetric(const double pos[2]){
      double rsq = pos[0]*pos[0];
      double sth = cos(pos[1]);

      gd[1][1] = rsq;

      gd[2][2] = rsq*sth*sth;

      gu[1][1] = 1.0/(rsq + epsilon);

      gu[2][2] = 1.0/(rsq*sth*sth + epsilon);

      det = rsq*rsq*sth*sth;
      sdet = sqrt(det);
    }

    virtual void updateChristoffelSymbols(const double pos[2]){
      double sth = cos(pos[1]);
      double cth = sin(pos[1]);

      christoffel[0][1][1] = -pos[0];
      christoffel[0][2][2] = -pos[0]*sth*sth;
      christoffel[1][0][1] = 1.0/(pos[0] + epsilon);
      christoffel[1][1][0] = christoffel[1][0][1];
      christoffel[1][2][2] = sth*cth;
      christoffel[2][0][2] = christoffel[1][0][1];
      christoffel[2][2][0] = christoffel[1][0][1];
      christoffel[2][1][2] = cth/(sth + epsilon);
      christoffel[2][2][1] = christoffel[2][1][2];
    }

    virtual double getLength(const double pos[2]){
      return pos[0];
    }
};
#endif
