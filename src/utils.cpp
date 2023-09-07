#include <vector>
#include <tuple>
#include <iostream>
#include "utils.h"
#include "particle_data_group.h"

double utils::get_norm_sq(four_vec vec)
{
    double norm = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        norm += vec[mu] * vec[mu] * utils::gmumu[mu];
    }
    return norm;
}

int utils::levi(int i, int j, int k, int l)
{
    // Levi-Civita symbols
    // i,j,k,l = 0...3 i.e. upper indices
    if ((i == j) || (i == k) || (i == l) || (j == k) || (j == l) || (k == l))
        return 0;
    else
        return ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l) / 12);
}

int utils::rand_int(int min, int max)
{
    static randomizer rng;
    static bool initialized = false;
    if (!initialized)
    {
        rng.seed(time(0));
        initialized = true;
    }
    std::uniform_int_distribution<uint32_t> dist(0, 10);
    return dist(rng);
}

int utils::g(int mu, int nu)
{
    if (mu > 3 || mu < 0 || nu > 3 || nu < 0)
    {
        std::cout << "error with the indices of the metric" << std::endl;
        exit(1);
    }
    if (mu == nu)
    {
        if (mu == 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
    return 0;
}

std::vector<double> utils::linspace(double min, double max, int size)
{
    // returns a vector of doubles from min to max (included), equally spaced with spacing max-min/size
    std::vector<double> result;
    double interval = (max - min) / size;
    double tmp_min = min;
    for (int i = 0; i <= size; i++)
    {
        result.push_back(tmp_min);
        tmp_min += interval;
    }
    return result;
}

utils::interpolator::interpolator(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> f) : x_(x), y_(y), z_(z), f_(f) {}

utils::interpolator::~interpolator() {}

std::tuple<double, double> utils::interpolator::adjacent_points(double x, std::vector<double> direction)
{
    //"x" is some point that we want to locate in the vector "direction"
    // gives for x in the interval [x0,x1] returns the tuple (x0,x1). If x coincides with one of
    // the grid points it is displaced by an infinitesimal amount to make the algorithm work.
    double dummy_x = x;
    if (dummy_x == direction[direction.size() - 1])
    {
        dummy_x -= 1e-8;
    }
    for (int i = 0; i < direction.size(); i++)
    {
        if (dummy_x > direction[i] & dummy_x < direction[i + 1])
        {
            return {direction[i], direction[i + 1]};
        }
        if (dummy_x == direction[i])
        {
            dummy_x += 1e-8;
        }
    }
    std::cout << "Out of bounds! ERROR" << std::endl;
    exit(1);
}

double utils::interpolator::evaluate_f_lattice(double x, double y, double z)
{
    // Gives the value of f at the point x,y,z in the table. If the point is not on the table it returns an error
    // Supposedly, the sizes of the vectors x,y and z are the same (e.g. they are read as columns of a file)
    for (int i = 0; i < x_.size(); i++)
    {
        if (x == x_[i] & y == y_[i] & z == z_[i])
        {
            return f_[i];
        }
    }
    std::cout << "Not on the lattice!" << std::endl;
    exit(1);
}

double utils::interpolator::trilinear_interpolation(double x, double y, double z)
{
    // algorithm taken from https://en.wikipedia.org/wiki/Trilinear_interpolation
    auto [x0, x1] = adjacent_points(x, x_);
    auto [y0, y1] = adjacent_points(y, y_);
    auto [z0, z1] = adjacent_points(z, z_);

    double c000 = evaluate_f_lattice(x0, y0, z0);
    double c001 = evaluate_f_lattice(x0, y0, z1);
    double c010 = evaluate_f_lattice(x0, y1, z0);
    double c100 = evaluate_f_lattice(x1, y0, z0);
    double c011 = evaluate_f_lattice(x0, y1, z1);
    double c110 = evaluate_f_lattice(x1, y1, z0);
    double c101 = evaluate_f_lattice(x1, y0, z1);
    double c111 = evaluate_f_lattice(x1, y1, z1);

    // std::cout<<c000;
    double dx = x0 - x1 + 1e-8;
    double dy = y0 - y1 + 1e-8;
    double dz = z0 - z1 + 1e-8;
    double denominator = dx * dy * dz;

    double a0 = (-c000 * x1 * y1 * z1 + c001 * x1 * y1 * z0 + c010 * x1 * y0 * z1 - c011 * x1 * y0 * z0) / denominator +
                (c100 * x0 * y1 * z1 - c101 * x0 * y1 * z0 - c110 * x0 * y0 * z1 + c111 * x0 * y0 * z0) / denominator;
    double a1 = (c000 * y1 * z1 - c001 * y1 * z0 - c010 * y0 * z1 + c011 * y0 * z0) / denominator +
                (-c100 * y1 * z1 + c101 * y1 * z0 + c110 * y0 * z1 - c111 * y0 * z0) / denominator;
    double a2 = (c000 * x1 * z1 - c001 * x1 * z0 - c010 * x1 * z1 + c011 * x1 * z0) / denominator +
                (-c100 * x0 * z1 + c101 * x0 * z0 + c110 * x0 * z1 - c111 * x0 * z0) / denominator;
    double a3 = (c000 * x1 * y1 - c001 * x1 * y1 - c010 * x1 * y0 + c011 * x1 * y0) / denominator +
                (-c100 * x0 * y1 + c101 * x0 * y1 + c110 * x0 * y0 - c111 * x0 * y0) / denominator;
    double a4 = (-c000 * z1 + c001 * z0 + c010 * z1 - c011 * z0 + c100 * z1 - c101 * z0 - c110 * z1 + c111 * z0) / denominator;
    double a5 = (-c000 * y1 + c001 * y1 + c010 * y0 - c011 * y0 + c100 * y1 - c101 * y1 - c110 * y0 + c111 * y0) / denominator;
    double a6 = (-c000 * x1 + c001 * x1 + c010 * x1 - c011 * x1 + c100 * x0 - c101 * x0 - c110 * x0 + c111 * x0) / denominator;
    double a7 = (c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) / denominator;

    return a0 + a1 * x + a2 * y + a3 * z + a4 * x * y + a5 * x * z + a6 * y * z + a7 * x * y * z;
}

template <typename T>
inline T utils::absolute_error(const T approx, const T exact)
{
    return abs(approx - exact);
}

template <typename T>
T utils::relative_error(const T approx, const T exact)
{
    return utils::absolute_error(approx, exact) / exact;
}
