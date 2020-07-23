#include <mpicommunicator.h>
#include <fluid.h>
#include <iostream>
#include <cstdio>
#include <primitive/idealsolver.h>
#include <primitive/enthalpysolver.h>
#include <riemann/hllesolver.h>
#include <geom/cartesianmetric.h>
#include <geom/polarmetric.h>
#include <geom/cylindricalmetric.h>
#include <geom/sphericalmetric.h>
#include <recon/norecon.h>
#include <recon/minmod.h>
#include <recon/mp5.h>
#include <rk4.h>
#include <cmath>
#include <fluidioparameters.h>
#include <fluidioparser.h>
#include <fluidgridparameters.h>
#include <fluidgridparser.h>
#include <fluidmethodparameters.h>
#include <fluidmethodparser.h>
#include <fluidicparameters.h>
#include <fluidicparser.h>

int main(int argc, char *argv[]){
  MPICommunicator *comm = MPICommunicator::getInstance();
  if(comm->init() != SUCCESS){
    std::cout << "There was an error initializing MPICommunicator.\n";
    return 0;
  }

  int rank = comm->getRank();
  int root = comm->getRootRank();

  // Load the parameters.
  FluidIOParameters iopar;
  FluidGridParameters gridpar;
  FluidMethodParameters methodpar;
  FluidICParameters icpar;

  if(rank == root){
    if(argc < 2){
      std::cout << "Usage: ./Fluid <parameter file>\n";
      MPI_Abort(MPI_COMM_WORLD, 0);
      return 0;
    }
    FluidIOParser ioparser;
    bool result = ioparser.updateParameters(argv[1], &iopar);
    if(!result){
      std::cout << "Could not load " << argv[1] << ". Terminating...\n";
      MPI_Abort(MPI_COMM_WORLD, 0);
      return 0;
    }

    // Load the rest of the parameters.
    std::string gridfile = iopar.getGridFile();
    if(gridfile.compare("NULL") != 0){
      FluidGridParser gridparser;
      result = gridparser.updateParameters(gridfile.c_str(), &gridpar);
      if(!result){
        std::cout << "Could not load " << gridfile << ". Terminating...\n";
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 0;
      }
    }
    std::string methodfile = iopar.getMethodFile();
    if(methodfile.compare("NULL") != 0){
      FluidMethodParser methodparser;
      result = methodparser.updateParameters(methodfile.c_str(), &methodpar);
      if(!result){
        std::cout << "Could not load " << methodfile << ". Terminating...\n";
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 0;
      }
    }
    std::string icfile = iopar.getICFile();
    if(icfile.compare("NULL") != 0){
      FluidICParser icparser;
      result = icparser.updateParameters(icfile.c_str(), &icpar);
      if(!result){
        std::cout << "Could not load " << icfile << ". Terminating...\n";
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 0;
      }
    }
  }

  // Broadcast all the parameters.
  iopar.broadcastParameters();
  gridpar.broadcastParameters();
  methodpar.broadcastParameters();
  icpar.broadcastParameters();

  // Assign the right metric to the problem.
  Metric *metric;
  Domain::GridCoordinates coords = Domain::CARTESIAN;
  double periodic = false;
  double pos[2] = {1.0, 1.0};
  switch(gridpar.getCoordinates()){
    case FluidGridParameters::CARTESIAN:
      metric = new CartesianMetric();
      if(rank == root){
        std::cout << "Cartesian metric selected.\n";
      }
      break;
    case FluidGridParameters::POLAR:
      metric = new PolarMetric(pos);
      coords = Domain::POLAR;
      if(gridpar.getBoundaryDown() == FluidGridParameters::DOWN_PERIODIC || 
         gridpar.getBoundaryUp() == FluidGridParameters::UP_PERIODIC){
        if(gridpar.getBoundaryDown() != FluidGridParameters::DOWN_PERIODIC ||
           gridpar.getBoundaryUp() != FluidGridParameters::UP_PERIODIC){
          if(rank == root){
            std::cout << "Warning! Periodic boundaries mismatched!\n";
            std::cout << "Adjusting mismatched boundary.\n";
          }
          gridpar.setBoundaryDown(FluidGridParameters::DOWN_PERIODIC);
          gridpar.setBoundaryUp(FluidGridParameters::UP_PERIODIC);
        }
        periodic = true;
      }
      if(rank == root){
        std::cout << "Polar metric selected.\n";
      }
      break;
    case FluidGridParameters::CYLINDRICAL:
      metric = new CylindricalMetric(pos);
      if(rank == root){
        std::cout << "Cylindrical metric selected.\n";
      }
      break;
    case FluidGridParameters::SPHERICAL:
      metric = new SphericalMetric(pos);
      coords = Domain::POLAR;
      if(gridpar.getBoundaryDown() == FluidGridParameters::DOWN_PERIODIC || 
         gridpar.getBoundaryUp() == FluidGridParameters::UP_PERIODIC){
        if(gridpar.getBoundaryDown() != FluidGridParameters::DOWN_PERIODIC ||
           gridpar.getBoundaryUp() != FluidGridParameters::UP_PERIODIC){
          if(rank == root){
            std::cout << "Warning! Periodic boundaries mismatched!\n";
            std::cout << "Adjusting mismatched boundary.\n";
          }
          gridpar.setBoundaryDown(FluidGridParameters::DOWN_PERIODIC);
          gridpar.setBoundaryUp(FluidGridParameters::UP_PERIODIC);
        }
        periodic = true;
      }
      if(rank == root){
        std::cout << "Spherical metric selected.\n";
      }
      break;
    default:
      if(rank == root){
        std::cout << "Unknown coordinate system requested. Defaulting to Cartesian.\n";
      }
      break;
  }
  // Set up the Domain.
  Domain domain = Domain();
  pair2<double> bounds;
  bounds[0][0] = gridpar.getDomainMinX();
  bounds[0][1] = gridpar.getDomainMaxX();
  bounds[1][0] = gridpar.getDomainMinY();
  bounds[1][1] = gridpar.getDomainMaxY();
  domain.setCFL(gridpar.getCFL());
  domain.setBounds(bounds);
  domain.setGhostPoints(3);

  unsigned int shp[2] = {(unsigned int)gridpar.getGridSizeX(), 
                         (unsigned int)gridpar.getGridSizeY()};
  // FIXME: Set up periodic conditions for polar coordinates.
  bool result = domain.buildMesh(shp,coords,periodic);

  // Set up the ODE.
  RK4 rk4 = RK4();
  // Get the reconstruction method.
  Recon *recon;
  switch(methodpar.getReconstruction()){
    case FluidMethodParameters::NONE:
      recon = new NoRecon();
      if(rank == root){
        std::cout << "No reconstruction enabled.\n";
      }
      break;
    case FluidMethodParameters::MINMOD:
      recon = new Minmod();
      if(rank == root){
        std::cout << "Minmod reconstruction selected.\n";
      }
      break;
    case FluidMethodParameters::WENO5:
      recon = new Minmod();
      if(rank == root){
        std::cout << "WENO5 not currently implemented. Reverting to minmod.\n";
      }
      break;
    case FluidMethodParameters::MP5:
      recon = new MP5();
      if(rank == root){
        std::cout << "MP5 reconstruction selected.\n";
      }
      break;
    default:
      recon = new NoRecon();
      if(rank == root){
        std::cout << "Unknown reconstruction method requested. Disabling reconstruction.\n";
      }
      break;
  }
  // Get the EOS solver.
  PrimitiveSolver *primitive;
  switch(methodpar.getEOS()){
    case FluidMethodParameters::FIXED:
      primitive = new IdealSolver(metric);
      ((IdealSolver*)primitive)->setGamma(methodpar.getGamma());
      if(rank == root){
        std::cout << "Fixed gamma EOS selected.\n";
      }
      break;
    case FluidMethodParameters::RYU:
      if(rank == root){
        std::cout << "Ryu EOS not currently implemented. Reverting to fixed gamma EOS.\n";
      }
      primitive = new IdealSolver(metric);
      ((IdealSolver*)primitive)->setGamma(methodpar.getGamma());
      break;
    case FluidMethodParameters::FIXED_ENTHALPY:
      primitive = new EnthalpySolver(metric);
      ((IdealSolver*)primitive)->setGamma(methodpar.getGamma());
      if(rank == root){
        std::cout << "Fixed gamma EOS with enthalpy solver selected.\n";
      }
      break;
    default:
      if(rank == root){
        std::cout << "Unknown EOS method requested. Reverting to fixed gamma EOS.\n";
      }
      primitive = new IdealSolver(metric);
      ((IdealSolver*)primitive)->setGamma(methodpar.getGamma());
  }
  // Get the Riemann solver.
  RiemannSolver *riemann;
  switch(methodpar.getRiemannSolver()){
    case FluidMethodParameters::HLLE:
      riemann = new HLLESolver(primitive,metric);
      if(rank == root){
        std::cout << "HLLE Riemann solver selected.\n";
      }
      break;
    default:
      if(rank == root){
        std::cout << "Unknown Riemann solver requested. Reverting to HLLE.\n";
      }
      riemann = new HLLESolver(primitive,metric);
      break;
  }
  Fluid ode = Fluid(&domain, &rk4);
  ode.setMetric(metric);
  ode.setParameters(&icpar);
  ode.setGridParameters(&gridpar);
  ode.setPrimitiveSolver(primitive);
  ode.setRiemannSolver(riemann);
  ode.setRecon(recon);
  ode.initData();
  ode.setVariableOutput("Primitive",0,true);
  ode.setVariableOutput("Primitive",1,true);
  ode.setVariableOutput("Primitive",2,true);
  ode.setVariableOutput("Primitive",3,true);
  ode.setVariableOutput("Primitive",4,true);
  ode.setVariableOutput("Conserved",0,true);
  ode.setVariableOutput("Conserved",1,true);
  ode.setVariableOutput("Conserved",2,true);
  ode.setVariableOutput("Conserved",3,true);
  ode.setVariableOutput("Conserved",4,true);

  double ti = 0;
  double tf = gridpar.getTime();

  auto dx = domain.getGrid()->getSpacing();
  double spacing = fmin(dx[0], dx[1]);
  double dt = domain.getCFL()*spacing;
  unsigned int M = (tf - ti)/dt;
  unsigned int output = iopar.getOutputFrequency();
  unsigned int print = iopar.getPrintFrequency();

  if(output != 0){
    ode.outputVTK("fluid00000", 0);
  }
  if(print != 0 && rank == root){
    std::cout << "Step: 0 Time: 0.0\n";
  }
  for(unsigned int i = 0; i < M; i++){
    double t = (i + 1)*dt;
    ode.evolveStep(dt);

    char buffer[20];
    sprintf(buffer,"fluid%05d",i+1);
    if(output != 0 && (i+1) % output == 0){
      ode.outputVTK(buffer, t);
    }
    
    if(print != 0 && (i+1) % print == 0 && rank == root){
      std::cout << "Step: " << i+1 << " Time: " << t << std::endl;
    }
  }

  delete primitive;
  delete riemann;
  delete recon;
  delete metric;

  if(comm->cleanup() != SUCCESS){
    std::cout << "There was an error cleaning up MPICommunicator.\n";
  }
  return 0;
}
