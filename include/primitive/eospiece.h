#ifndef EOS_PIECE_H
#define EOS_PIECE_H

/**************************************************************************
 *
 * Class: EOSPiece
 * Author: Jacob Fields
 * Date Modified: 4-6-2021
 *
 * Description: A very simple class containing information about a section
 *   of the polytropic EOS.
 *
 *************************************************************************/

class EOSPiece{
  private:
    // The threshold density for this piece.
    const double rho;
    const double kappa;
    const double gamma;
    const double a;

    double pressure;
    double enthalpy;
  public:
    EOSPiece(double r, double k, double g, double c);

    inline double getRho() const{
      return rho;
    }
    inline double getKappa() const{
      return kappa;
    }
    inline double getGamma() const{
      return gamma;
    }
    inline double getA() const{
      return a;
    }

    double getTrialSsq(double D, double tau);

    double calcH(double dens, double P);

    double calcSoundSpeed(double dens, double P);
};

#endif
