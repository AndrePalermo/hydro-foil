#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <iostream>
#include <array>
#include "utils.h"
using namespace std;

// meant to read a the beta.dat file from the output of vhlle (polarization branches)

struct element
{
    double tau, x, y, eta;
    vector<double> u;
    vector<double> dsigma;
    double T, mub, muq, mus;
    vector<vector<double>> dbeta;
    element() : u(4, 0), dsigma(4, 0), dbeta(4, vector<double>(4, 0)){};
    void print();
    double get_normal_size();

    bool is_timelike()
    {
        return get_normal_size() > 0;
    }

private:
    double _normal_size = 0;
};

void read_hypersrface(std::string filename, vector<element> &hypersurface);

#endif