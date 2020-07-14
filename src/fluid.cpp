#include <fluid.h>
#include <geom/metric.h>
#include <riemann/riemannsolver.h>
#include <primitive/primitivesolver.h>
#include <recon/recon.h>
#include <iostream>
#include <cmath>
#include <mpicommunicator.h>
#include <fluidicparameters.h>
#include <fluidgridparameters.h>

// Constructor {{{
Fluid::Fluid(Domain *d, Solver *s) : ODE(5, d, s){
  metric = nullptr;
  riemann = nullptr;
  primitive = nullptr;
  recon = nullptr;

  // Add the variables we need.
  addField("Conserved", 5, true, true, 1);
  addField("Primitive", 5, false, false, 4);

  setVariableName("Conserved", U_D  , "U_D");
  setVariableName("Conserved", U_SX , "U_SX");
  setVariableName("Conserved", U_SY , "U_SY");
  setVariableName("Conserved", U_SZ , "U_SZ");
  setVariableName("Conserved", U_TAU, "U_TAU");

  setVariableName("Primitive", V_RHO, "V_RHO");
  setVariableName("Primitive", V_VX , "V_VX");
  setVariableName("Primitive", V_VY , "V_VY");
  setVariableName("Primitive", V_VZ , "V_VZ");
  setVariableName("Primitive", V_P  , "V_P");

  reallocateData();
}
// }}}

// Destructor {{{
Fluid::~Fluid(){
  
}
// }}}

// rhs {{{
void Fluid::rhs(std::shared_ptr<FieldMap>& fieldMap){
  auto conservedVars = fieldMap->getSolverField("Conserved");
  auto primitiveVars = (*fieldMap)["Primitive"];
  const Grid& grid = fieldMap->getGrid();
  unsigned int nx = grid.getSize()[DIR_X];
  unsigned int ny = grid.getSize()[DIR_Y];
  auto points = grid.getPoints();
  unsigned int nb = domain->getGhostPoints();

  // Retrieve the data from our fields for easy access.
  double **u   = conservedVars->getIntermediateData();
  double **dtu = conservedVars->getCurrentRHS();
  double **f   = conservedVars->getLine(0); // Flux line
  double **v   = primitiveVars->getData();
  double **v1  = primitiveVars->getLine(0); // Line in v
  double **vl  = primitiveVars->getLine(1); // Reconstructed left line
  double **vr  = primitiveVars->getLine(2); // Reconstructed right line
  double **wv  = primitiveVars->getLine(3); // Four velocity for reconstruction.

  // Calculate the flux in the x and y directions.
  // FIXME: These should be moved into new functions to make the code cleaner.
  // x-direction
  double dx = grid.getSpacing()[DIR_X];
  double dy = grid.getSpacing()[DIR_Y];
  double dt = domain->getCFL()*fmin(dx, dy);
  
  double offset = 0.5*dx;
  double pos[2];
  for(unsigned int j = nb; j < ny - nb; j++){
    pos[DIR_Y] = points[DIR_Y][j];
    // Get a line in constant y and reconstruct it.
    for(unsigned int m = 0; m < NU; m++){
      for(unsigned int i = 0; i < nx; i++){
        v1[m][i] = v[m][grid.getIndex(i,j)];
      }
    }

    // Apply reconstruction here.
    reconstruct(vl, vr, v1, wv, DIR_X, grid, j);

    // Calculate the conserved variables at the cell interface.
    for(unsigned int i = 0; i < nx; i++){
      double ulpt[NU], urpt[NU], vlpt[NV], vrpt[NV];
      pos[DIR_X] = points[DIR_X][i] - offset;

      for(unsigned int m = 0; m < NV; m++){
        vlpt[m] = vl[m][i];
        vrpt[m] = vr[m][i];
      }

      metric->updateMetric(pos);
      primitive->primToConPt(ulpt, vlpt);
      primitive->primToConPt(urpt, vrpt);

      /*for(unsigned int m = 0; m < NU; m++){
        ul[m][i] = ulpt[m];
        ur[m][i] = urpt[m];
      }*/

      // Solve the Riemann problem at this point.
      double fpt[NU] = {0.0};
      double flpt[NU] = {0.0};
      double frpt[NU] = {0.0};
      riemann->calcFluxX(flpt, ulpt, vlpt);
      riemann->calcFluxX(frpt, urpt, vrpt);

      riemann->solveRiemannProblem(fpt, flpt, frpt, ulpt, urpt, vlpt, vrpt, DIR_X);
      // DEBUGGING ONLY: Lax-Friedrichs flux
      /*for(unsigned int m = 0; m < NU; m++){
        fpt[m] = 0.5*(flpt[m] + frpt[m]) - dx*(urpt[m] - ulpt[m])/(dt*2.0);
      }*/
      for(unsigned int m = 0; m < NU; m++){
        f[m][i] = fpt[m];
      }
    }
    // Calculate the x-derivative of the flux.
    for(unsigned int m = 0; m < NU; m++){
      for(unsigned int i = 0; i < nx-1; i++){
        unsigned int pp = grid.getIndex(i, j);
        dtu[m][pp] = -(f[m][i+1] - f[m][i])/dx;
      }
    }
  }

  // y-direction
  offset = 0.5*dy;
  for(unsigned int i = nb; i < nx-nb; i++){
    pos[DIR_X] = points[DIR_X][i];
    // Get a line in constant x and reconstruct it.
    for(unsigned int m = 0; m < NU; m++){
      for(unsigned int j = 0; j < ny; j++){
        v1[m][j] = v[m][grid.getIndex(i,j)];
      }
    }

    // Apply reconstruction here.
    reconstruct(vl, vr, v1, wv, DIR_Y, grid, i);

    // Calculate the conserved variables at the cell interface.
    for(unsigned int j = 0; j < ny; j++){
      double ulpt[NU], urpt[NU], vlpt[NV], vrpt[NV];
      pos[DIR_Y] = points[DIR_Y][j] - offset;

      for(unsigned int m = 0; m < NV; m++){
        vlpt[m] = vl[m][j];
        vrpt[m] = vr[m][j];
      }

      metric->updateMetric(pos);
      primitive->primToConPt(ulpt, vlpt);
      primitive->primToConPt(urpt, vrpt);

      /*for(unsigned int m = 0; m < NU; m++){
        ul[m][i] = ulpt[m];
        ur[m][i] = urpt[m];
      }*/

      // Solve the Riemann problem at this point.
      double fpt[NU] = {0.0};
      double flpt[NU] = {0.0};
      double frpt[NU] = {0.0};
      riemann->calcFluxY(flpt, ulpt, vlpt);
      riemann->calcFluxY(frpt, urpt, vrpt);

      riemann->solveRiemannProblem(fpt, flpt, frpt, ulpt, urpt, vlpt, vrpt, DIR_Y);
      // DEBUGGING ONLY: Lax-Friedrichs flux
      /*for(unsigned int m = 0; m < NU; m++){
        fpt[m] = 0.5*(flpt[m] + frpt[m]) - dy*(urpt[m] - ulpt[m])/(dt*2.0);
      }*/
      for(unsigned int m = 0; m < NU; m++){
        f[m][j] = fpt[m];
      }
    }
    // Calculate the y-derivative of the flux.
    for(unsigned int m = 0; m < NU; m++){
      for(unsigned int j = 0; j < ny-1; j++){
        unsigned int pp = grid.getIndex(i, j);
        dtu[m][pp] -= (f[m][j+1] - f[m][j])/dy;
      }
    }
  }

  // Apply the source terms.
  calculateSourceTerms(dtu, u, v, grid);
}
// }}}

// initData {{{
void Fluid::initData(){
  switch(params->getInitialConditions()){
    case FluidICParameters::GAUSSIAN:
      applyGaussian();
      break;
    case FluidICParameters::SHOCKTUBE_X:
      applyShockTube1D(DIR_X);
      break;
    case FluidICParameters::SHOCKTUBE_Y:
      applyShockTube1D(DIR_Y);
      break;
    default:
      applyGaussian();
      break;
  }
}
// }}}

// applyBoundaries {{{
void Fluid::applyBoundaries(){
  double **u = fieldData->getSolverField("Conserved")->getIntermediateData();
  unsigned int nb = domain->getGhostPoints();
  const Grid *grid = domain->getGrid();
  unsigned int nx = grid->getSize()[0];
  unsigned int ny = grid->getSize()[1];
  if(domain->hasBoundary(LEFT)){
    // Apply outflow boundaries.
    switch(gridParams->getBoundaryLeft()){
      case FluidGridParameters::LEFT_OUTFLOW:
        for(unsigned int m = 0; m < NU; m++){
          for(unsigned int j = 0; j < ny; j++){
            unsigned int pnb = grid->getIndex(nb,j);
            for(unsigned int i = 0; i < nb; i++){
              unsigned int pp = grid->getIndex(i,j);
              u[m][pp] = u[m][pnb];
            }
          }
        }
        break;
      case FluidGridParameters::LEFT_WALL:
        // In the case of a wall, we treat D, SY, SZ, and TAU
        // as even functions and SX as an odd function.
        for(unsigned int m = 0; m < NU; m++){
          double sign = (m != U_SX) ? 1.0 : -1.0;
          for(unsigned int j = 0; j < ny; j++){
            for(unsigned int p = 0; p < nb; p++){
              u[m][grid->getIndex(nb - p - 1, j)] = sign * u[m][grid->getIndex(nb + p + 1,j)];
            }
          }
        }
        break;
    }
  }
  if(domain->hasBoundary(RIGHT)){
    // Apply outflow boundaries.
    switch(gridParams->getBoundaryRight()){
      case FluidGridParameters::RIGHT_OUTFLOW:
        for(unsigned int m = 0; m < NU; m++){
          for(unsigned int j = 0; j < ny; j++){
            unsigned int pnb = grid->getIndex(nx-nb-1,j);
            for(unsigned int i = nx-nb; i < nx; i++){
              unsigned int pp = grid->getIndex(i,j);
              u[m][pp] = u[m][pnb];
            }
          }
        }
        break;
      case FluidGridParameters::RIGHT_WALL:
        for(unsigned int m = 0; m < NU; m++){
          double sign = (m != U_SX) ? 1.0 : -1.0;
          for(unsigned int j = 0; j < ny; j++){
            for(unsigned int p = 0; p < nb; p++){
              u[m][grid->getIndex(nx - nb + p, j)] = sign * u[m][grid->getIndex(nx - nb - 2 - p, j)];
            }
          }
        }
        break;
    }
  }
  if(domain->hasBoundary(DOWN)){
    // Apply outflow boundaries.
    switch(gridParams->getBoundaryDown()){
      case FluidGridParameters::DOWN_OUTFLOW:
        for(unsigned int m = 0; m < NU; m++){
          for(unsigned int j = 0; j < nb; j++){
            for(unsigned int i = 0; i < nx; i++){
              unsigned int pp = grid->getIndex(i,j);
              unsigned int pnb = grid->getIndex(i,nb);
              u[m][pp] = u[m][pnb];
            }
          }
        }
        break;
      case FluidGridParameters::DOWN_WALL:
        for(unsigned int m = 0; m < NU; m++){
          double sign = (m != U_SY) ? 1.0 : -1.0;
          for(unsigned int p = 0; p < nb; p++){
            for(unsigned int i = 0; i < nx; i++){
              u[m][grid->getIndex(i, nb - p - 1)] = sign * u[m][grid->getIndex(i, nb + p + 1)];
            }
          }
        }
        break;
    }
  }
  if(domain->hasBoundary(UP)){
    // Apply outflow boundaries.
    switch(gridParams->getBoundaryUp()){
      case FluidGridParameters::UP_OUTFLOW:
        for(unsigned int m = 0; m < NU; m++){
          for(unsigned int j = ny-nb; j < ny; j++){
            for(unsigned int i = 0; i < nx; i++){
              unsigned int pp = grid->getIndex(i,j);
              unsigned int pnb = grid->getIndex(i,ny-nb-1);
              u[m][pp] = u[m][pnb];
            }
          }
        }
        break;
      case FluidGridParameters::UP_WALL:
        for(unsigned int m = 0; m < NU; m++){
          double sign = (m != U_SY) ? 1.0 : -1.0;
          for(unsigned int p = 0; p < nb; p++){
            for(unsigned int i = 0; i < nx; i++){
              u[m][grid->getIndex(i, ny - nb + p)] = sign * u[m][grid->getIndex(i, ny - nb - p - 2)];
            }
          }
        }
        break;
    }
  }
}
// }}}

// doAfterBoundaries {{{
void Fluid::doAfterBoundaries(){
  // Retrieve the variables that we need from the available data.
  double **u = fieldData->getSolverField("Conserved")->getIntermediateData();
  double **v = (*fieldData)["Primitive"]->getData();
  const Grid *grid = domain->getGrid();
  unsigned int nx = grid->getSize()[0];
  unsigned int ny = grid->getSize()[1];
  const double *x = grid->getPoints()[0];
  const double *y = grid->getPoints()[1];

  double pos[2];
  double upt[NU], vpt[NV];

  // We need to get the primitive variables from the conserved.
  for(unsigned int j = 0; j < ny; j++){
    pos[1] = y[j];
    for(unsigned int i = 0; i < nx; i++){
      pos[0] = x[i];
      unsigned int pp = grid->getIndex(i,j);

      metric->updateMetric(pos);

      // The conserved variables
      upt[U_D  ] = u[U_D  ][pp];
      upt[U_SX ] = u[U_SX ][pp];
      upt[U_SY ] = u[U_SY ][pp];
      upt[U_SZ ] = u[U_SZ ][pp];
      upt[U_TAU] = u[U_TAU][pp];

      // The primitive variables
      vpt[V_RHO] = v[V_RHO][pp];
      vpt[V_VX ] = v[V_VX ][pp];
      vpt[V_VY ] = v[V_VY ][pp];
      vpt[V_VZ ] = v[V_VZ ][pp];
      vpt[V_P  ] = v[V_P  ][pp];

      bool result = primitive->conToPrimPt(upt, vpt);
      if(!result){
        std::cout << "An error occurred during conToPrimPt! Aborting...\n";
        MPI_Abort(MPI_COMM_WORLD,0);
      }

      // Copy the results back into the data.
      u[U_D  ][pp] = upt[U_D  ];
      u[U_SX ][pp] = upt[U_SX ];
      u[U_SY ][pp] = upt[U_SY ];
      u[U_SZ ][pp] = upt[U_SZ ];
      u[U_TAU][pp] = upt[U_TAU];

      v[V_RHO][pp] = vpt[V_RHO];
      v[V_VX ][pp] = vpt[V_VX ];
      v[V_VY ][pp] = vpt[V_VY ];
      v[V_VZ ][pp] = vpt[V_VZ ];
      v[V_P  ][pp] = vpt[V_P  ];
    }
  }
  #ifdef DEBUG_PRIMITIVE_SOLVER
  std::cout << "Average number of Newton iterations: " << primitive->getNewtonAverage() <<std::endl;
  primitive->resetNewtonAverage();
  #endif

  // Handle axis treatment here.
}
// }}}

// reconstruct {{{
void Fluid::reconstruct(double *vl[], double *vr[], double *v1[], double *wv[],
                        const unsigned int dir, const Grid& grid, const unsigned int j){
  double pos[2];
  unsigned int n;
  const double *x;

  const unsigned int otherdir = 1 - dir;
  n = grid.getSize()[dir];
  x = grid.getPoints()[dir];
  pos[otherdir] = grid.getPoints()[otherdir][j];
  double offset = 0.5*grid.getSpacing()[dir];

  recon->reconstruct(n, v1[V_RHO], vl[V_RHO], vr[V_RHO]);
  recon->reconstruct(n, v1[V_P  ], vl[V_P  ], vr[V_P  ]);

  // Reconstruct the four-velocity instead of the three velocity.
  calculateFourVelocity(wv, v1, dir, grid, j);
  recon->reconstruct(n, wv[0], vl[V_VX], vr[V_VX]);
  recon->reconstruct(n, wv[1], vl[V_VY], vr[V_VY]);
  recon->reconstruct(n, wv[2], vl[V_VZ], vr[V_VZ]);
  for(unsigned int i = 0; i < n; i++){
    pos[dir] = x[i] - offset;
    double wvlpt[3], wvrpt[3];

    wvlpt[0] = vl[V_VX][i];
    wvlpt[1] = vl[V_VY][i];
    wvlpt[2] = vl[V_VZ][i];

    wvrpt[0] = vr[V_VX][i];
    wvrpt[1] = vr[V_VY][i];
    wvrpt[2] = vr[V_VZ][i];

    metric->updateMetric(pos);

    double wl = sqrt(1.0 + metric->squareVector(wvlpt));
    double wr = sqrt(1.0 + metric->squareVector(wvrpt));
    // These are commented out because they should be mathematically unnecessary.
    // However, if there are issues, we can put it back.
    // wl = fmax(1.0, wl);
    // wr = fmax(1.0, wr);

    vl[V_VX][i] = wvlpt[0]/wl;
    vl[V_VY][i] = wvlpt[1]/wl;
    vl[V_VZ][i] = wvlpt[2]/wl;

    vr[V_VX][i] = wvrpt[0]/wr;
    vr[V_VY][i] = wvrpt[1]/wr;
    vr[V_VZ][i] = wvrpt[2]/wr;
  }
  
  double vacuum = primitive->getVacuum();
  for(unsigned int i = 0; i < n; i++){
    if(vl[V_RHO][i] < vacuum){
      vl[V_RHO][i] = vacuum;
      vl[V_VX][i] = 0.0;
      vl[V_VY][i] = 0.0;
      vl[V_VZ][i] = 0.0;
    }
    if(vr[V_RHO][i] < vacuum){
      vr[V_RHO][i] = vacuum;
      vr[V_VX][i] = 0.0;
      vr[V_VY][i] = 0.0;
      vr[V_VZ][i] = 0.0;
    }
    if(vl[V_P  ][i] < vacuum){
      vl[V_P  ][i] = vacuum;
    }
    if(vr[V_P  ][i] < vacuum){
      vr[V_P  ][i] = vacuum;
    }
  }
}
// }}}

// calculateFourVelocity {{{
void Fluid::calculateFourVelocity(double *wv[], double *v[], const unsigned int dir, 
                                  const Grid& grid, const unsigned int j){
  double pos[2];
  unsigned int n;
  const double *x;
  const unsigned int otherdir = 1 -dir;
  const double alpha = 1.0;
  const double beta[3] = {0.0};

  pos[otherdir] = grid.getPoints()[otherdir][j];
  n = grid.getSize()[dir];
  x = grid.getPoints()[dir];
  double offset = 0.5*grid.getSpacing()[dir];

  double vu[3];
  for(unsigned int i = 0; i < n; i++){
    pos[dir] = x[i] - offset;
    metric->updateMetric(pos);

    vu[0] = v[V_VX][i];
    vu[1] = v[V_VY][i];
    vu[2] = v[V_VZ][i];
    double vsq = metric->squareVector(vu);

    // Check if the velocity is superluminal.
    double t1 = vsq - 1.0;
    if(t1 >= 0){
      // Check if it's a roundoff error.
      if(t1 < 1e-15){
        while( t1 >= 0.0 && t1 < 1e-15){
          for(unsigned int m = 0; m < 3; m++){
            vu[m] -= copysign(1.0,vu[m])*5.0e-16;
          }
          vsq = metric->squareVector(vu);
          t1 = vsq - 1.0;
        }
      }
      else{
        // Something has gone horribly wrong.
        // We'll try to fix it.
        printf("calculateFourVelocity: Before reconstruction, vsq = %25.20e.\n",vsq);
        printf("  v = ( %g, %g, %g )\n", vu[0], vu[1], vu[2]);
        double t1 = sqrt(1.0/vsq - 6.0e-16);
        for(unsigned int m = 0; m < 3; m++){
          vu[m] *= t1;
        }
        vsq = metric->squareVector(vu);
        // If that failed, scream and quit.
        if(vsq >= 1.0){
          printf("calculateFourVelocity: Something's still wrong with vsq. Aborting...\n");
          MPI_Abort(MPI_COMM_WORLD, 0);
        }
      }
    }

    double W = sqrt(1.0/(1.0 - vsq));

    wv[0][i] = W*(vu[0] - beta[0]/alpha);
    wv[1][i] = W*(vu[1] - beta[1]/alpha);
    wv[2][i] = W*(vu[2] - beta[2]/alpha);
  }
}
// }}}

// calculateSourceTerms {{{
void Fluid::calculateSourceTerms(double *dtu[], double *u[], double *v[], const Grid& grid){
  unsigned int nx = grid.getSize()[0];
  unsigned int ny = grid.getSize()[1];

  double pos[2] = {0.0};
  double vu[3], Sd[3];
  double P;

  for(unsigned int j = 0; j < ny; j++){
    for(unsigned int i = 0; i < nx; i++){
      unsigned int pp = grid.getIndex(i, j);

      // Get the data for this point.
      vu[0] = v[V_VX][pp];
      vu[1] = v[V_VY][pp];
      vu[2] = v[V_VZ][pp];
      Sd[0] = u[U_SX][pp];
      Sd[1] = u[U_SY][pp];
      Sd[2] = u[U_SZ][pp];
      P = v[V_P][pp];

      // Update the metric based on the position.
      pos[0] = grid.getPoints()[0][i];
      pos[1] = grid.getPoints()[1][j];
      metric->updateMetric(pos);
      metric->updateChristoffelSymbols(pos);
      double sdetg = metric->getDeterminantRoot();
      double ***conn = metric->getChristoffelSymbols();
      double Pdens = sdetg*P;
      // General source term calculation.
      for(unsigned int k = 0; k < 3; k++){
        for(unsigned int l = 0; l < 3; l++){
          double arg = vu[l]*Sd[k] + Pdens*(k==l);
          dtu[U_SX][pp] += conn[k][l][0]*arg;
          dtu[U_SY][pp] += conn[k][l][1]*arg;
          dtu[U_SZ][pp] += conn[k][l][2]*arg;
        }
      }
    }
  }
}
// }}}

// applyGaussian {{{
void Fluid::applyGaussian(){
  double **v = (*fieldData)["Primitive"]->getData();
  double **u = (*fieldData)["Conserved"]->getData();
  const Grid *grid = domain->getGrid();
  unsigned int nx = grid->getSize()[0];
  unsigned int ny = grid->getSize()[1];
  const double *x = grid->getPoints()[0];
  const double *y = grid->getPoints()[1];
  double pos[2];
  
  double rsq;

  double kappa = params->getGaussianKappa();
  double A = params->getGaussianAmplitude();
  double sigma = params->getGaussianSigma();
  double gamma = 5.0/3.0;
  double vacuum = primitive->getVacuum();
  double vpt[NV], upt[NU];
  for(unsigned int j = 0; j < ny; j++){
    for(unsigned int i = 0; i < nx; i++){
      unsigned int pp = grid->getIndex(i,j);
      pos[0] = x[i];
      pos[1] = y[j];
      metric->updateMetric(pos);
      // Calculate r. Right now, assume Cartesian coordinates.
      // FIXME: Generalize this.
      double r = metric->getLength(pos);
      rsq = r*r;
      v[V_RHO][pp] = vacuum + A*exp(-rsq/(sigma*sigma));
      v[V_VX ][pp] = 0.0;
      v[V_VY ][pp] = 0.0;
      v[V_VZ ][pp] = 0.0;
      v[V_P  ][pp] = vacuum + kappa*pow(v[V_RHO][pp], gamma);

      for(unsigned int m = 0; m < NV; m++){
        vpt[m] = v[m][pp];
      }

      // Calculate the conserved variables.
      primitive->primToConPt(upt, vpt);

      for(unsigned int m = 0; m < NU; m++){
        u[m][pp] = upt[m];
      }
    }
  } 
}
// }}}

// applyShockTube1D {{{
void Fluid::applyShockTube1D(const unsigned int dir){
  double **v = (*fieldData)["Primitive"]->getData();
  double **u = (*fieldData)["Conserved"]->getData();
  const Grid *grid = domain->getGrid();
  auto shp = grid->getSize();
  auto points = grid->getPoints();
  const unsigned int otherdir = 1 - dir;
  double pos[2];
  double center = 0.0;
  double rho_left = 1.0;
  double rho_right = 1.0;
  double p_left = 1000.0;
  double p_right = 0.1;
  double v_left = 0.0;
  double v_right = 0.0;
  double vpt[NV], upt[NU];

  for(unsigned int j = 0; j < shp[DIR_Y]; j++){
    pos[DIR_Y] = points[DIR_Y][j];
    for(unsigned int i = 0; i < shp[DIR_X]; i++){
      pos[DIR_X] = points[DIR_X][i];
      const unsigned int pp = grid->getIndex(i,j);
      // Check if we're to the right of the center.
      if(pos[dir] > center){
        v[V_RHO][pp] = rho_right;
        v[V_P  ][pp] = p_right;
        v[V_VX + dir][pp] = v_right;
        v[V_VX + otherdir][pp] = 0.0;
        v[V_VZ][pp] = 0.0;
      }
      else{
        v[V_RHO][pp] = rho_left;
        v[V_P  ][pp] = p_left;
        v[V_VX + dir][pp] = v_left;
        v[V_VX + otherdir][pp] = v_right;
        v[V_VZ][pp] = 0.0;
      }
      for(unsigned int m = 0; m < NV; m++){
        vpt[m] = v[m][pp];
      }
      primitive->primToConPt(upt, vpt);
      for(unsigned int m = 0; m < NU; m++){
        u[m][pp] = upt[m];
      }
    }
  }
}
// }}}