#include <iostream>
#include <vector>
#include <string>
using namespace std;

int main() {
    string algon = "xy";
    string integrator = " omf4";
    string D = " 2";
    string L = " 20";
    vector<int> niter = {2, 3, 4, 5, 6, 7, 10, 15, 20};    
    int ntherminit = 200000;
    int nstepinit = 400000;
    int nmeasinit = 1;
    double betaorkappa = 0.90;
    string fname;

    for(int i =0; i < niter.size(); i++) {
        if (niter[i] < 10) {
            fname = string(" ") + algon + string("0") + to_string(niter[i]) + string("leapf.dat");
        } else {
            fname = string(" ") + algon + to_string(niter[i]) + string("leapf.dat");
        }
        
        string command = string("./simu ") + algon + integrator + D + L + string(" ") + to_string(niter[i]) 
                + string(" ") + to_string(ntherminit) + string(" ") + to_string(nstepinit) 
                + string(" ") + to_string(nmeasinit) + string(" ") + to_string(betaorkappa) + fname;
        system(command.c_str());
    }

    return 0;
}