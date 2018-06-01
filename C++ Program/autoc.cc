#include <iterator>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include<random>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <libgen.h>
using namespace std;

double rho(vector<double> x, double xave, int ndat, int t)
{
    int n, t0;
    double r = 0.0;
    n = ndat - t;
    
    // Auto correlation
    for (t0 = 0; t0 < n; t0++)
        r += (x[t0] - xave) * (x[t0 + t] - xave);
    r /= n;
    return r;
}

int main(int argc, char **argv)
{
    string fname; stringstream(argv[1]) >> fname;
    cout << fname << endl;
    int tmax; stringstream(argv[2]) >> tmax;
    ifstream is(fname);
    vector<double> x;
    string line;
    while (getline(is, line))
    {
        if (line.find("#") != string::npos){}
        else {
            string tmp = line.substr(line.find("  ") + 1);
            x.push_back(std::stod(tmp));
        }
    }

    ofstream myfile;
    string autofname = string("auto") + fname;
    myfile.open(autofname);

    double norm;
    int i, ndat, t, tcut;
    double O, dO, chi, dchi;
    ndat = x.size();
    double xave = 0;

    for (int t0 = 0; t0 < x.size(); t0++)
    {
        xave += x[t0];
    }
    xave /= (double)x.size();
    cout << xave << endl;
    vector<double> r(tmax);
    for (t = 0; t < tmax; t++)
        r[t] = rho(x, xave, ndat, t);
    norm = 1.0 / r[0];
    for (t = 0; t < tmax; t++)
        r[t] *= norm; 
    //tau[t] stores integrated autocorrelation times with tcut=t
    vector<double> tau(tmax);
    tau[0] = 0.5;
    for (tcut = 1; tcut < tmax; tcut++)
    {
        for (t = 0; t <= tcut; t++)
            tau[tcut] += r[t]; //sum of r[t]
        std::uniform_real_distribution<double> adist(0,tcut);
    }

    for (t = 0; t < tmax; t++)
        myfile << r[t] << "  " << tau[t] << '\n';

    return 0;
}