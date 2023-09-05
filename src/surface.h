#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <iostream>
#include <array>
#include "utils.h"

// meant to read a the beta.dat file from the output of vhlle (polarization branches)
namespace surface
{
    struct element
    {
        double tau, x, y, eta;
        std::vector<double> u;
        std::vector<double> dsigma;
        double T, mub, muq, mus;
        std::vector<std::vector<double>> dbeta;
        element() : u(4, 0), dsigma(4, 0), dbeta(4, std::vector<double>(4, 0)){};
        void print();
        double get_normal_size();

        bool is_timelike()
        {
            return get_normal_size() > 0;
        }

        friend std::istream &operator>>(std::istream &stream, element &cell);
        friend std::ostream &operator<<(std::ostream &stream, element &cell);

    private:
        double _normal_size = 0;
    };

    // void read_hypersrface(std::string filename, std::vector<element> &hypersurface);
    void read_hypersrface(std::ifstream &file, std::vector<element> &hypersurface);

    std::array<double, 4> get_mins(const std::vector<element> &hypersup);
    std::array<double, 4> get_maxs(const std::vector<element> &hypersup);
}
#endif