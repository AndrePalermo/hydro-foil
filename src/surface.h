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
        double t();
        double z();
        double get_normal_size_sq();
        utils::four_vec get_pos();
        utils::four_vec get_pos_mink();

        bool is_timelike()
        {
            return get_normal_size_sq() > 0;
        }

        friend std::istream &operator>>(std::istream &stream, element &cell);
        friend std::ostream &operator<<(std::ostream &stream, element &cell);
        friend bool operator==(const element &lhs, const element &rhs);
        friend bool operator!=(const element &lhs, const element &rhs);

        double distance(element other);

        static double distance(element first, element second);

    private:
        double _normal_size = 0;
    };

    // void read_hypersrface(std::string filename, std::vector<element> &hypersurface);
    void read_hypersrface(std::ifstream &file, std::vector<element> &hypersurface);

    std::array<double, 4> get_mins(const std::vector<element> &hypersup);
    std::array<double, 4> get_maxs(const std::vector<element> &hypersup);
    // Find smallest delta_x^\mu to a reference point
    std::array<double, 4> get_deltas(const std::vector<element> &hypersup, element refpoint, std::array<element, 4> &points);
    element get_random_point(const std::vector<element> &hypersup, int& index);
    int get_random_path(const std::vector<element> &hypersup, int start, std::array<double, 4> deltas, int count, std::vector<element> &path);
    std::vector<element> get_diff_tau(const std::vector<element> &hypersup, const element &refpoint);
    std::vector<element> get_diff_x(const std::vector<element> &hypersup, const element &refpoint);
    std::vector<element> get_diff_y(const std::vector<element> &hypersup, const element &refpoint);
    std::vector<element> get_diff_eta(const std::vector<element> &hypersup, const element &refpoint);
}
#endif