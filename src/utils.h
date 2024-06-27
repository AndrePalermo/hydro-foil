#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <cmath>
#include <vector>
#include <chrono>


const std::array<int,4> gmumu = {1,-1,-1,-1};
const std::array<int,4> t_vector = {1,0,0,0};
int g(int mu, int nu);

const double Gevtofm = 5.067728853;
const double hbarC = 1. / 5.067728853; //=0.197 Gev*fm
const double PI = std::acos(-1);

int levi(int i, int j, int k, int l);
std::vector<double> linspace(double min, double max, int size);
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> file_to_4_column(std::string filename, int f_column);

//class for interpolation
class interpolator{
    private:
    public:
    std::vector<double> x_, y_, z_, f_; //x,y,z vector of a table of the values x. All quantities have the same lenght
    int nx, ny, nz;
    double xmaxG;
    double ymaxG;
    double zmaxG;
    double xminG;
    double yminG;
    double zminG;
    double dxG;
    double dyG;
    double dzG;
        interpolator(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> f);
        interpolator(std::string filename, int f_column);
        ~interpolator();

        double trilinear_interpolation(double x, double y, double z);
        double trilinear_interpol(double x, double y, double z);
        double f_index(int ix, int iy, int iz);

};

constexpr std::array<std::array<int,4>, 24> NonZeroIndicesLeviCivita() {
    ///List of arrays of indices for which the Levi Civita symbol is non-zero
    return {{
        {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}, {0, 1, 3, 2}, {0, 3, 2, 1}, {0, 2, 3, 1},
        {1, 0, 2, 3}, {1, 2, 0, 3}, {1, 3, 0, 2}, {1, 0, 3, 2}, {1, 2, 3, 0}, {1, 3, 2, 0},
        {2, 0, 1, 3}, {2, 1, 0, 3}, {2, 3, 0, 1}, {2, 0, 3, 1}, {2, 1, 3, 0}, {2, 3, 1, 0},
        {3, 0, 1, 2}, {3, 1, 0, 2}, {3, 2, 0, 1}, {3, 0, 2, 1}, {3, 1, 2, 0}, {3, 2, 1, 0},       
    }};
}

#endif
