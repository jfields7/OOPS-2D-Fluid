#ifndef EOSPOLICY_H
#define EOSPOLICY_H

class EOSPolicy{
  private:
    // Ranges for the policy.
    const double minRho;
    const double maxRho;
    const double minEps;
    const double maxEps;
    
    void applyRhoMaxPolicy(double *u, double *v);
    void applyRhoMinPolicy(double *u, double *v);

    void applyEpsMaxPolicy(double *u, double *v);
    void applyEpsMinPolicy(double *u, double *v);

  public:
    // A constructor for the EOSPolicy.
    EOSPolicy(double minRho_, double maxRho_, double minEps_, double maxEps_);

    void applyPolicy(double *u, double *v);
};

#endif
