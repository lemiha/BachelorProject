#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<random>
#include "ownfuncs.h"
#include "forces.h"
#include "integrators.h"
#include "algorithms.h"

int main(int argc, char** argv) {
    if (argc != 11) {
		std::cerr << "usage: " << argv[0] <<" <algo> <integrator> <D> <L> <niter> <ntherm> <nstep> <nmeas> <betaorkappa> <fname>" << std::endl;
		exit(1);
	}

    std::string algoname; std::stringstream(argv[1]) >> algoname;
    std::string nameofintegrator; std::stringstream(argv[2]) >> nameofintegrator;    
	int D;         std::stringstream(argv[3]) >> D;
    int L;         std::stringstream(argv[4]) >> L;
	int niter;    std::stringstream(argv[5]) >> niter;
	int ntherm; std::stringstream(argv[6]) >> ntherm;
	int nstep;  std::stringstream(argv[7]) >> nstep;
	int nmeas;  std::stringstream(argv[8]) >> nmeas;
    double betaorkappa; std::stringstream(argv[9]) >> betaorkappa;
    std::string fname; std::stringstream(argv[10]) >> fname;

    BasicAlgo * algo = NULL;
    if (algoname == "hmc") {
        algo = new HMC();
    } else if (algoname == "hmcn") {
        algo = new HMC_NoMetro();
    } else if (algoname == "smd") {
        algo = new SMD();
    } else if (algoname == "smdn") {
        algo = new SMD_NoMetro();
    } else if (algoname == "xy") {
        algo = new XY();
    } else if (algoname == "xyn") {
        algo = new XY_NoMetro();
    } else if (algoname == "xysmd") {
        algo = new XY_SMD();
    } else if (algoname == "xysmdn") {
        algo = new XY_SMD_NoMetro();
    } else {
        std::cerr << "No algorithm with the given name" << std::endl;
    }
    BasicAlgo &Algo = *algo;

    Force * f = NULL;
    Integrator * integrate = NULL;

    if (algoname == "hmc" || "hmcn" || "smd" || "smdn") {
        f = new ForceScalarField();
        if (nameofintegrator == "omf2") {
            integrate = new OMF2Integrator();            
        } else {
            integrate = new OMF4Integrator();            
        }
    }

    if (algoname == "xy" || "xyn" || "xysmd" || "xysmdn") {
        f = new ForceXY(); 
        if (nameofintegrator == "omf2") {
            integrate = new OMF2IntegratorXY();
        } else {
            integrate = new OMF4IntegratorXY();
        }
    }

    Force * Fo = f;
    Integrator * integrator = integrate; 

    // Timer start
    std::clock_t starttime;
    starttime = std::clock();     

    Algo.setVariables(nameofintegrator, starttime, D, L, niter, ntherm, nstep, nmeas, fname, Fo, integrator, betaorkappa);
    Algo.simulate();
    
    return 0;
}