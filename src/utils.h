#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <cmath>
#include <vector>


const std::array<int,4> gmumu = {1,-1,-1,-1};
const std::array<int,4> t_vector = {1,0,0,0};

const double Gevtofm = 5.067728853;
const double hbarC = 1. / 5.067728853; //=0.197...
const double PI = std::acos(-1);


//std::vector<std::vector<int> gmatrix = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
int levi(int i, int j, int k, int l);
std::vector<double> linspace(double min, double max, int size);

template<size_t ndim>
class interpolation{
    private:
        std::array<size_t, ndim> dim;
        std::array<double, ndim> cordinate_max, cordinate_min;
        std::vector<double> grid_values;

    public:
        interpolation(std::array<double, ndim> cordinate_max, std::array<double, ndim> cordinate_min, std::vector<double> &grid_values);       
        ~interpolation();

        // interpolation read_grid(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> grid_values);
        double interpolate(double x_point, double y_point, double z_point);
}

#endif
