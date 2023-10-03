#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <cmath>
#include <vector>


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
        std::vector<double> x_, y_, z_, f_; //x,y,z vector of a table of the values x. All quantities have the same lenght

    public:
        interpolator(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> f);
        interpolator(std::string filename, int f_column);
        ~interpolator();

        double trilinear_interpolation(double x, double y, double z);
};

#endif
