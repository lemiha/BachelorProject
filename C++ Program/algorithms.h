#ifndef __algorithms_h__
#define __algorithms_h__

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <random>
#include <algorithm>
#include "ownfuncs.h"
#include "forces.h"
#include "integrators.h"

// Basic skeleton for the algorithms
class BasicAlgo
{
    public:
      void setVariables(std::string nameofintegrator_in, std::clock_t& starttime_in, int &D_in, int &L_in, int &niter_in, int &ntherm_in,
                        int &nstep_in, int &nmeas_in, std::string &fname_in,
                        Force *&F_in, Integrator *int_in, double betaorkappa_in)
      {
            nameofintegrator = nameofintegrator_in;
            starttime = starttime_in;
            D = D_in;
            L = L_in;
            niter = niter_in;
            ntherm = ntherm_in;
            nstep = nstep_in;
            nmeas = nmeas_in;
            fname = fname_in;
            F = F_in;
            integrator = int_in;

            sigma.assign(pow(L, D), 0.0);
            pi.assign(pow(L, D), 0.0);
            r = {0,2,4,6,8,10};

            hop.resize(pow(L, D), std::vector<int>(D * 2, 0));
            hopping(L, D, hop);

            lambda = 0.5;
            kappa = betaorkappa_in;
            beta = betaorkappa_in;
            // beta = 1.0;
            // kappa = 1.3;
            // lambda = 1.1689;
            // kappa = 0.185825;
            tau = 1.0;
            eps = tau / double(niter);

            acc = 0.0;
            seed = 8475565;
            gen.seed(seed);
            myfile.open(fname);
            corrfile.open(std::string("corr") + fname);
            myfile << "#"
                   << " D=" << D << " L=" << L << " niter=" << niter << " ntherm=" << ntherm << " nstep=" << nstep << " nmeas=" << nmeas << std::endl;
            myfile << "#"
                   << " lambda=" << lambda << " betaorkappa=" << kappa << " tau=" << tau << " eps=" << eps << " integrator=" << nameofintegrator << std::endl;
      }

      void writeToFile(int t, double Eprint)
      {
            if ((t >= ntherm) && (((t - ntherm) % nmeas) == 0))
            {
                  std::vector<double> rsum(r.size()-1);

                  for (int x = 0; x < pow(L, D); x++)
                  {
                        std::vector<double> rneig(D, x);

                        for (int i = 1; i < r.size(); i++)
                        {
                              for (int j = r[i - 1]; j < r[i]; j++)
                              {
                                    for (int k = 0; k < rneig.size(); k++)
                                    {
                                          // std::cout << rneig[k] << std::endl;
                                          rneig[k] = hop[rneig[k]][k];
                                    }
                              }

                              for (int k = 0; k < rneig.size(); k++)
                              {
                                    rsum[i-1] += sigma[x]*sigma[rneig[k]];
                              }
                        }
                  }
                  
                  for (int i = 0; i < rsum.size(); i++) {
                        if (i == rsum.size()-1)
                              corrfile << rsum[i] / (double)(D*pow(L, D)) << std::endl;
                        else {
                              corrfile << rsum[i] / (double)(D*pow(L, D)) << "  ";
                        }
                  }
                  myfile << Eprint << std::endl;
            }

      }

      void writeToFileXY(int t, double Eprint)
      {
            if ((t >= ntherm) && (((t - ntherm) % nmeas) == 0))
            {
                  std::vector<double> rsum(r.size()-1);

                  for (int x = 0; x < pow(L, D); x++)
                  {
                        std::vector<double> rneig(D, x);

                        for (int i = 1; i < r.size(); i++)
                        {
                              for (int j = r[i - 1]; j < r[i]; j++)
                              {
                                    for (int k = 0; k < rneig.size(); k++)
                                    {
                                          // std::cout << rneig[k] << std::endl;
                                          rneig[k] = hop[rneig[k]][k];
                                    }
                              }

                              for (int k = 0; k < rneig.size(); k++)
                              {
                                    rsum[i-1] += cos(sigma[x] - sigma[rneig[k]]);
                              }
                        }
                  }
                  
                  for (int i = 0; i < rsum.size(); i++) {
                        if (i == rsum.size()-1)
                              corrfile << rsum[i] / (double)(D*pow(L, D)) << std::endl;
                        else {
                              corrfile << rsum[i] / (double)(D*pow(L, D)) << "  ";
                        }
                  }

                  myfile << Eprint << std::endl;
            }

      }

      void writeAccTime(double& acc) {
            double duration;
            duration = ( std::clock() - starttime ) / (double) CLOCKS_PER_SEC;
            myfile << "#"
                   << " acc=" << acc / (nstep + ntherm)
                   << " time=" << duration                   
                   << std::endl;
            myfile.close();
      }

      virtual void simulate() = 0;

    protected:
      std::clock_t starttime;
      int D;      // Dimensions
      int L;      // Lattice size
      int niter;  // Nr. of integration steps
      int ntherm; // Nr. of thermalization steps
      int nstep;  // Nr. of updates
      int nmeas;  //Frequency of measurements
      std::string nameofintegrator;
      Force *F;
      Integrator *integrator;

      std::vector<double> sigma;
      std::vector<double> pi;
      std::vector<int> r;
      int nb; // Nr. of neighbours
      std::vector<std::vector<int>> hop;

      double lambda;
      double kappa;
      double beta;
      double tau;
      double eps;
      unsigned seed;
      double acc; // Acceptance ratio
      std::mt19937 gen;
      std::ofstream myfile;
      std::ofstream corrfile;
      std::string fname;
};

class HMC : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);
            double Etempnew, Eprint;
            double Etempold = F->calcPotential(sigma);

            myfile << "#"
                   << " algo=hmc" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  std::vector<double> sigma_new = sigma;

                  double Eold = 0;
                  for (int x = 0; x < pow(L, D); x++) {
                        pi[x] = sdist(gen);
                        Eold += 0.5 * pi[x] * pi[x];
                  }

                  Eold += Etempold;

                  integrator->setInitialConditions(sigma_new, pi);
                  for (int t = 0; t < niter; t++)
                        integrator->integrate(sigma_new, pi);

                  double Enew = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        Enew += 0.5 * pi[x] * pi[x];

                  Etempnew = F->calcPotential(sigma_new);
                  Enew += Etempnew;

                  double deltaE = Enew - Eold;
                  double p = exp(-1.0 * (deltaE));

                  if (adist(gen) < p)
                  {
                        acc++;
                        sigma = sigma_new;
                        Etempold = Etempnew;
                        Eprint = Etempnew;
                  }
                  else
                  {
                        Eprint = Etempold;
                  }

                  writeToFile(t, Eprint);
            }

            writeAccTime(acc);
      }
};

class HMC_NoMetro : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);

            double Etemp;

            myfile << "#"
                   << " algo=hmcn" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        pi[x] = sdist(gen);
                  }

                  integrator->setInitialConditions(sigma, pi);
                  for (int t = 0; t < niter; t++)
                        integrator->integrate(sigma, pi);

                  double E = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        E += 0.5 * pi[x] * pi[x];

                  Etemp = F->calcPotential(sigma);
                  E += Etemp;

                  writeToFile(t, Etemp);
            }

            writeAccTime(acc);
      }
};

class SMD : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);
            double gamma = 0.6;
            double c1 = exp(-gamma * eps);
            double c2 = sqrt(1 - c1 * c1);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);

            double Etempnew, Eprint;
            double Etempold = F->calcPotential(sigma);

            for (int x = 0; x < pow(L, D); x++)
            {
                  pi[x] = sdist(gen);
            }

            myfile << "#"
                   << " algo=smd" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  std::vector<double> sigma_new = sigma;

                  double Eold = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        pi[x] = c1 * pi[x] + c2 * sdist(gen);
                        Eold += 0.5 * pi[x] * pi[x];
                  }

                  std::vector<double> pi_new = pi;
                  Eold += Etempold;

                  integrator->setInitialConditions(sigma_new, pi_new);
                  // for (int t=0;t<niter;t++)
                  integrator->integrate(sigma_new, pi_new);

                  double Enew = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        Enew += 0.5 * pi_new[x] * pi_new[x];

                  Etempnew = F->calcPotential(sigma_new);
                  Enew += Etempnew;

                  double deltaE = Enew - Eold;
                  double p = exp(-1.0 * (deltaE));

                  if (adist(gen) < p)
                  {
                        acc++;
                        sigma = sigma_new;
                        pi = pi_new;
                        Etempold = Etempnew;
                        Eprint = Etempnew;
                  }
                  else
                  {
                        for (int x = 0; x < pow(L, D); x++)
                        {
                              pi[x] = -pi[x];
                        }
                        Eprint = Etempold;
                  }

                  writeToFile(t, Eprint);
            }

            writeAccTime(acc);
      }
};

class SMD_NoMetro : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);
            double gamma = 0.6;
            double c1 = exp(-gamma * eps);
            double c2 = sqrt(1 - c1 * c1);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);

            for (int x = 0; x < pow(L, D); x++)
            {
                  pi[x] = sdist(gen);
            }

            double Etemp;

            myfile << "#"
                   << " algo=smdn" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  // double v = sdist(gen);
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        // pi[x]=c1*pi[x] + c2*v;
                        pi[x] = c1 * pi[x] + c2 * sdist(gen);
                  }

                  integrator->setInitialConditions(sigma, pi);
                  // for (int t=0;t<niter;t++)
                  integrator->integrate(sigma, pi);

                  double E = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        E += 0.5 * pi[x] * pi[x];

                  double Etemp = F->calcPotential(sigma);
                  E += Etemp;

                  writeToFile(t, Etemp);
            }

            writeAccTime(acc);
      }
};

class XY : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);
            double Etempnew, Eprint;
            double Etempold = F->calcPotential(sigma);

            myfile << "#"
                   << " algo=xy" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  std::vector<double> sigma_new = sigma;

                  double Eold = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        pi[x] = sdist(gen);
                        Eold += 0.5 * pi[x] * pi[x];
                  }

                  Eold += Etempold;

                  integrator->setInitialConditions(sigma_new, pi);
                  for (int t = 0; t < niter; t++)
                        integrator->integrate(sigma_new, pi);

                  double Enew = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        Enew += 0.5 * pi[x] * pi[x];

                  Etempnew = F->calcPotential(sigma_new);
                  Enew += Etempnew;

                  double deltaE = Enew - Eold;
                  double p = exp(-1.0 * (deltaE));

                  if (adist(gen) < p)
                  {
                        acc++;
                        sigma = sigma_new;
                        Etempold = Etempnew;
                        Eprint = Etempnew;
                  }
                  else
                  {
                        Eprint = Etempold;
                  }

                  writeToFileXY(t, Eprint);
            }

            writeAccTime(acc);
      }
};

class XY_NoMetro : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);

            double Etemp;

            myfile << "#"
                   << " algo=xyn" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        pi[x] = sdist(gen);
                  }

                  integrator->setInitialConditions(sigma, pi);
                  for (int t = 0; t < niter; t++)
                        integrator->integrate(sigma, pi);

                  double E = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        E += 0.5 * pi[x] * pi[x];

                  Etemp = F->calcPotential(sigma);
                  E += Etemp;

                  writeToFileXY(t, Etemp);
            }

            writeAccTime(acc);
      }
};

class XY_SMD : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);
            double gamma = 0.6;
            double c1 = exp(-gamma * eps);
            double c2 = sqrt(1 - c1 * c1);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);

            double Etempnew, Eprint;
            double Etempold = F->calcPotential(sigma);

            for (int x = 0; x < pow(L, D); x++)
            {
                  pi[x] = sdist(gen);
            }

            myfile << "#"
                   << " algo=xysmd" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  std::vector<double> sigma_new = sigma;

                  double Eold = 0.0;
                  // double v = sdist(gen);
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        pi[x] = c1 * pi[x] + c2 * sdist(gen);
                        Eold += 0.5 * pi[x] * pi[x];
                        // pi[x]=c1*pi_new[x] + c2*v;
                  }

                  std::vector<double> pi_new = pi;
                  Eold += Etempold;

                  integrator->setInitialConditions(sigma_new, pi_new);
                  // for (int t=0;t<niter;t++)
                  integrator->integrate(sigma_new, pi_new);

                  double Enew = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        Enew += 0.5 * pi_new[x] * pi_new[x];

                  Etempnew = F->calcPotential(sigma_new);
                  Enew += Etempnew;

                  double deltaE = Enew - Eold;
                  double p = exp(-1.0 * (deltaE));
                  double pEn;

                  if (adist(gen) < p)
                  {
                        acc++;
                        sigma = sigma_new;
                        pi = pi_new;
                        Etempold = Etempnew;
                        Eprint = Etempnew;
                  }
                  else
                  {
                        for (int x = 0; x < pow(L, D); x++)
                        {
                              pi[x] = -pi[x];
                        }
                        Eprint = Etempold;
                  }

                  writeToFileXY(t, Eprint);
            }

            writeAccTime(acc);
      }
};

class XY_SMD_NoMetro : public BasicAlgo
{
    public:
      void simulate()
      {
            std::uniform_real_distribution<double> adist(0, 1.0);
            std::normal_distribution<double> sdist(0, 1.0);
            double gamma = 0.6;
            double c1 = exp(-gamma * eps);
            double c2 = sqrt(1 - c1 * c1);

            F->setVariables(L, D, lambda, kappa, beta, hop);
            integrator->setVariables(F, eps);

            for (int x = 0; x < pow(L, D); x++)
            {
                  pi[x] = sdist(gen);
            }

            double Etemp;

            myfile << "#"
                   << " algo=xysmdn" << std::endl;

            for (int t = 0; t < (ntherm + nstep); t++)
            {
                  // double v = sdist(gen);
                  for (int x = 0; x < pow(L, D); x++)
                  {
                        // pi[x]=c1*pi[x] + c2*v;
                        pi[x] = c1 * pi[x] + c2 * sdist(gen);
                  }

                  integrator->setInitialConditions(sigma, pi);
                  // for (int t=0;t<niter;t++)
                  integrator->integrate(sigma, pi);

                  double E = 0.0;
                  for (int x = 0; x < pow(L, D); x++)
                        E += 0.5 * pi[x] * pi[x];

                  double Etemp = F->calcPotential(sigma);
                  E += Etemp;

                  writeToFileXY(t, Etemp);
            }

            writeAccTime(acc);
      }
};

#endif