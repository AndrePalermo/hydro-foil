#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <iostream>
#include <array>
using namespace std;

//meant to read a the beta.dat file from the output of vhlle (polarization branches)

struct element {
    double tau, x, y, eta;
    vector<double> u;
    vector<double> dsigma;
    double T, mub, muq, mus;
    vector<vector<double>> dbeta;
    element() : u(4, 0), dsigma(4, 0), dbeta(4, vector<double>(4, 0)){};
    void print();
};

void read_hypersrface(std::string filename, vector<element> &hypersurface);

#endif