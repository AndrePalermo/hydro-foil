#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <cmath>
#include <vector>
#include <random>

template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N> const& vec)
{
	os << "[";
    std::for_each(begin(vec), end(vec) - 1, [&os](T val) {os << val << ","; });
	os << vec[vec.size() - 1] << "]";
    return os;
}

namespace utils
{
    // Random engine type
    typedef std::mt19937 randomizer;
    // Generates random integer
    int rand_int(int min = 0, int max = 10);

    typedef std::array<double, 4> four_vec;

    const std::array<int, 4> gmumu = {1, -1, -1, -1};
    const std::array<int, 4> t_vector = {1, 0, 0, 0};
    int g(int mu, int nu);

    const double Gevtofm = 5.067728853;
    const double hbarC = 1. / 5.067728853; //=0.197 Gev*fm
    const double PI = std::acos(-1);

    double get_norm_sq(four_vec vec);

    int levi(int i, int j, int k, int l);
    std::vector<double> linspace(double min, double max, int size);

    // class for interpolation
    class interpolator
    {
    private:
        std::vector<double> x_, y_, z_, f_; // x,y,z vector of a table of the values x. All quantities have the same lenght

        double evaluate_f_lattice(double x, double y, double z);
        std::tuple<double, double> adjacent_points(double x, std::vector<double> direction);

    public:
        interpolator(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> f);
        ~interpolator();

        double trilinear_interpolation(double x, double y, double z);
    };

    template <typename T>
    T absolute_error(const T approx, const T exact);

    template <typename T>
    T relative_error(const T approx, const T exact);

}

#endif
