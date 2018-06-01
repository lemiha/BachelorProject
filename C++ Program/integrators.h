#ifndef __integrators_h__
#define __integrators_h__

#include <vector>
#include "forces.h"
#include "ownfuncs.h"

class Integrator {
    protected:
	    Force* F;
        double epsilon;
        int ndof;
        std::vector<double> qt; 
	    std::vector<double> pt; 

	public: 
        void setVariables(Force* F_in, double& eps_in) {
            F = F_in;
            epsilon = eps_in;
			ndof = F -> getNDof();
        }

		void setInitialConditions(const std::vector<double>& q0, 
				const std::vector<double>& p0) {
			qt=q0; pt=p0;
		}

        virtual void integrate(std::vector<double>& qtp1, std::vector<double>& ptp1) =0;
};

class OMF2Integrator: public Integrator {
	static constexpr double xi=0.1931833275037836;

	public: 
		void integrate(std::vector<double>& qtp1, std::vector<double>& ptp1) {
			// std::cout << "asdasd" << std::endl;
			
			int ndof = F -> getNDof(); 

			for (int i=0; i<ndof; i++) 
				qt[i] += xi*epsilon*pt[i];
			for (int i=0; i<ndof; i++) 
				pt[i] += 0.5*epsilon*(F -> forceEnergy(i,qt)); 	
			for (int i=0; i<ndof; i++) 
				qt[i] += (1.0-2.0*xi)*epsilon*pt[i];
			for (int i=0; i<ndof; i++) 
				pt[i] += 0.5*epsilon*(F -> forceEnergy(i,qt));
			for (int i=0; i<ndof; i++) 
				qt[i] += xi*epsilon*pt[i];

			qtp1=qt; ptp1=pt;
		}
};

class OMF2IntegratorXY: public Integrator {
	static constexpr double xi=0.1931833275037836;

	public: 
		void integrate(std::vector<double>& qtp1, std::vector<double>& ptp1) {
			int ndof = F -> getNDof(); 
			// std::cout << "asdasd" << std::endl;

			for (int i=0; i<ndof; i++) {
				qt[i] += xi*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			}
			for (int i=0; i<ndof; i++)
				pt[i] += 0.5*epsilon*(F -> forceEnergy(i,qt)); 	
			for (int i=0; i<ndof; i++) {
				qt[i] += (1.0-2.0*xi)*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			}
			for (int i=0; i<ndof; i++) 
				pt[i] += 0.5*epsilon*(F -> forceEnergy(i,qt));
			for (int i=0; i<ndof; i++) {
				qt[i] += xi*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			}

			qtp1=qt; ptp1=pt;
		}
};

class OMF4Integrator: public Integrator {
	static constexpr double xi=0.1644986515575760;
	static constexpr double lambda=-0.02094333910398989;
	static constexpr double chi=1.235692651138917;

	public: 
		void integrate(std::vector<double>& qtp1, std::vector<double>& ptp1) {
			int ndof = F -> getNDof();

			for (int i=0; i<ndof; i++) 
				qt[i] += xi*epsilon*pt[i];
			for (int i=0; i<ndof; i++)
				pt[i] += 0.5*(1.0-2.0*lambda)*epsilon*(F -> forceEnergy(i,qt)); 	
			for (int i=0; i<ndof; i++) 
				qt[i] += chi*epsilon*pt[i];
			for (int i=0; i<ndof; i++) 
				pt[i] += lambda*epsilon*(F -> forceEnergy(i,qt));
			for (int i=0; i<ndof; i++) 
				qt[i] += (1.0-2.0*(chi+xi))*epsilon*pt[i];
			for (int i=0; i<ndof; i++) 
				pt[i] += lambda*epsilon*(F -> forceEnergy(i,qt));
			for (int i=0; i<ndof; i++) 
				qt[i] += chi*epsilon*pt[i];
			for (int i=0; i<ndof; i++) 
				pt[i] += 0.5*(1.0-2.0*lambda)*epsilon*(F -> forceEnergy(i,qt)); 	
			for (int i=0; i<ndof; i++) 
				qt[i] += xi*epsilon*pt[i];

			qtp1=qt; ptp1=pt;
		}
};

class OMF4IntegratorXY: public Integrator {
	static constexpr double xi=0.1644986515575760;
	static constexpr double lambda=-0.02094333910398989;
	static constexpr double chi=1.235692651138917;

	public: 
		void integrate(std::vector<double>& qtp1, std::vector<double>& ptp1) {
			int ndof = F -> getNDof();

			for (int i=0; i<ndof; i++) {
				qt[i] += xi*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			}
			for (int i=0; i<ndof; i++) 
				pt[i] += 0.5*(1.0-2.0*lambda)*epsilon*(F -> forceEnergy(i,qt)); 	
			for (int i=0; i<ndof; i++) {
				qt[i] += chi*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			}
			for (int i=0; i<ndof; i++) 
				pt[i] += lambda*epsilon*(F -> forceEnergy(i,qt));
			for (int i=0; i<ndof; i++) {
				qt[i] += (1.0-2.0*(chi+xi))*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			} 
			for (int i=0; i<ndof; i++) 
				pt[i] += lambda*epsilon*(F -> forceEnergy(i,qt));
			for (int i=0; i<ndof; i++) {
				qt[i] += chi*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			} 
			for (int i=0; i<ndof; i++) 
				pt[i] += 0.5*(1.0-2.0*lambda)*epsilon*(F -> forceEnergy(i,qt)); 	
			for (int i=0; i<ndof; i++) {
				qt[i] += xi*epsilon*pt[i];
				if (qt[i]<0)
					qt[i]+=2.0*M_PI;
				else if (qt[i]>(2.0*M_PI))
					qt[i]-=2.0*M_PI;
			} 
			qtp1=qt; ptp1=pt;
		}
};

#endif 