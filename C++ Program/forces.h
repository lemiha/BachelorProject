#ifndef __forces_h__
#define __forces_h__

#include <iostream>
#include <iomanip>
#include <sstream>
#include <random>
#include "ownfuncs.h"

class Force
{
  public:
	// Constructor
	void setVariables(int &L_in, int &D_in, double &lambda_in, double &kappa_in, double &beta_in, std::vector<std::vector<int>> &hop_in)
	{
		L = L_in;
		D = D_in;
		lambda = lambda_in;
		kappa = kappa_in;
		ndof = pow(L, D);
		beta = beta_in;
		hop = hop_in;
	}

	// Nr. of degrees of freedom
	int getNDof()
	{
		return ndof;
	}

	virtual double forceEnergy(int i, std::vector<double> &sigma) = 0;
	virtual double calcPotential(std::vector<double> &sigma) = 0;

  protected:
	int L;
	int D;
	double kappa;
	double lambda;
	double beta;
	int ndof;
	std::vector<std::vector<int>> hop;
};

class ForceScalarField : public Force
{
  public:
	double forceEnergy(int i, std::vector<double> &sigma)
	{
		double sigman;

		// Calculate the neighbours of each site
		sigman = 0;
		for (int mu = 0; mu < 2 * D; mu++)
			sigman += sigma[hop[i][mu]];

		// Potential for each site combined into a total potential
		return 2*kappa*sigman - 2*sigma[i] - lambda*4*(sigma[i]*sigma[i] - 1)*sigma[i];
	}

	double calcPotential(std::vector<double> &sigma)
	{
		double sigman, en, sigma2;

		en = 0;

		// Calculate the total potential
		for (int i = 0; i < pow(L, D); i++)
		{
			// Calculate the forward neighbours of each site
			sigman = 0;
			for (int mu = 0; mu < D; mu++)
				sigman += sigma[hop[i][mu]];

			// Potential for each site combined into a total potential
			sigma2 = sigma[i] * sigma[i];
			en += -2*kappa*sigman*sigma[i] + sigma2 + lambda*(sigma2-1.0)*(sigma2-1.0);
		}

		return en;
	}
};

class ForceXY : public Force
{
	double forceEnergy(int i, std::vector<double> &sigma)
	{
		// Calculate the neighbours of each site
		double en = 0.0;
		for (int mu = 0; mu < 2 * D; mu++)
			en += sin(sigma[i] - sigma[hop[i][mu]]);

		return -beta * en;
	}

	double calcPotential(std::vector<double> &sigma)
	{
		double pot = 0.0;

		for (int i = 0; i < L * L; i++)
		{
			// Calculate the forward neighbours of each site
			double en = 0.0;
			for (int mu = 0; mu < D; mu++)
				en += cos(sigma[i] - sigma[hop[i][mu]]);

			pot += en;
		}

		return -beta * pot;
	}
};

#endif