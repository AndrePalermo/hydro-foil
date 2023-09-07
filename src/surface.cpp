#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "surface.h"

using namespace surface;

std::istream &surface::operator>>(std::istream &stream, element &cell)
{
    stream >> cell.tau >> cell.x >> cell.y >> cell.eta;
    for (int mu = 0; mu < 4; mu++)
    {
        stream >> cell.dsigma[mu];
    }
    for (int mu = 0; mu < 4; mu++)
    {
        stream >> cell.u[mu];
    }
    stream >> cell.T >> cell.mub >> cell.muq >> cell.mus;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream >> cell.dbeta[i][j];
        }
    }
    return stream;
}

std::ostream &surface::operator<<(std::ostream &stream, element &cell)
{
    stream << cell.tau << "\t" << cell.x << "\t" << cell.y << "\t"
           << cell.eta << "\t";
    for (int mu = 0; mu < 4; mu++)
    {
        stream << cell.dsigma[mu] << "\t";
    }
    for (int mu = 0; mu < 4; mu++)
    {
        stream << cell.u[mu] << "\t";
    }

    stream << cell.T << "\t" << cell.mub << "\t" << cell.muq << "\t" << cell.mus << "\t";

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            stream << cell.dbeta[i][j] << "\t";
        }
    }
    return stream;
}

bool surface::operator==(const element &lhs, const element &rhs)
{
    return lhs.tau == rhs.tau && lhs.x == rhs.x && lhs.y == rhs.y && lhs.eta == rhs.eta;
}

bool surface::operator!=(const element &lhs, const element &rhs)
{
    return lhs.tau != rhs.tau && lhs.x != rhs.x && lhs.y != rhs.y && lhs.eta != rhs.eta;
}

void surface::read_hypersrface(std::ifstream &input_file, std::vector<element> &hypersurface)
{
    std::string line;

    while (std::getline(input_file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue; // skip empty or comment lines
        }

        std::istringstream iss(line);
        element cell;

        iss >> cell;

        hypersurface.push_back(cell);
    }
}

std::array<double, 4> surface::get_mins(const std::vector<element> &hypersup)
{
    auto min_tau = std::min_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.tau < second.tau;
                                    })
                       ->tau;

    auto min_x = std::min_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.x < second.x;
                                  })
                     ->x;

    auto min_y = std::min_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.y < second.y;
                                  })
                     ->y;

    auto min_eta = std::min_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.eta < second.eta;
                                    })
                       ->eta;
    return {min_tau, min_x, min_y, min_eta};
}

std::array<double, 4> surface::get_maxs(const std::vector<element> &hypersup)
{
    auto max_tau = std::max_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.tau < second.tau;
                                    })
                       ->tau;

    auto max_x = std::max_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.x < second.x;
                                  })
                     ->x;

    auto max_y = std::max_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.y < second.y;
                                  })
                     ->y;

    auto max_eta = std::max_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.eta < second.eta;
                                    })
                       ->eta;
    return {max_tau, max_x, max_y, max_eta};
}

std::array<double, 4> surface::get_deltas(const std::vector<element> &hypersup, element refpoint, std::array<element, 4> &points)
{
    auto diff = surface::get_diff_tau(hypersup, refpoint);

    auto nearest_tau = std::min_element(diff.begin(), diff.end(),
                                        [refpoint](const element &first, const element &second)
                                        {
                                            return (abs(first.tau - refpoint.tau) < abs(second.tau - refpoint.tau));
                                        });

    diff = surface::get_diff_x(hypersup, refpoint);

    auto nearest_x = std::min_element(diff.begin(), diff.end(), [refpoint](const element &first, const element &second)
                                      { return abs(first.x - refpoint.x) < abs(second.x - refpoint.x); });

    diff = surface::get_diff_y(hypersup, refpoint);

    auto nearest_y = std::min_element(diff.begin(), diff.end(), [refpoint](const element &first, const element &second)
                                      { return abs(first.y - refpoint.y) < abs(second.y - refpoint.y); });

    diff = surface::get_diff_eta(hypersup, refpoint);

    auto nearest_eta = std::min_element(diff.begin(), diff.end(), [refpoint](const element &first, const element &second)
                                        { return abs(first.eta - refpoint.eta) < abs(second.eta - refpoint.eta); });

    points = {*nearest_tau, *nearest_x, *nearest_y, *nearest_eta};

    return {abs(refpoint.tau - nearest_tau->tau),
            abs(refpoint.x - nearest_x->x),
            abs(refpoint.y - nearest_y->y),
            abs(refpoint.eta - nearest_eta->eta)};
}

element surface::get_random_point(const std::vector<element> &hypersup, int &index)
{
    index = utils::rand_int(0, hypersup.size());
    return hypersup[index];
}

int surface::get_random_path(const std::vector<element> &hypersup, int start,
                             std::array<double, 4> deltas, int count, std::vector<element> &path)
{
    int actual = 0;
    element current = hypersup[start];
    for (size_t i = 0; i < count; i++)
    {
        /* code */
    }

    return 0;
}

std::vector<element> surface::get_diff_tau(const std::vector<element> &hypersup, const element &refpoint)
{
    std::vector<element> diff_tau;
    std::copy_if(hypersup.begin(), hypersup.end(), std::back_inserter(diff_tau), 
    [refpoint] (element e) { return e.tau != refpoint.tau && e.x == refpoint.x && e.y == refpoint.y && e.eta == refpoint.eta; });
    return diff_tau;
}

std::vector<element> surface::get_diff_x(const std::vector<element> &hypersup, const element &refpoint)
{
    std::vector<element> diff;
    std::copy_if(hypersup.begin(), hypersup.end(), std::back_inserter(diff), 
    [refpoint] (element e) { return e.tau == refpoint.tau && e.x != refpoint.x && e.y == refpoint.y && e.eta == refpoint.eta; });
    return diff;
}

std::vector<element> surface::get_diff_y(const std::vector<element> &hypersup, const element &refpoint)
{
    std::vector<element> diff;
    std::copy_if(hypersup.begin(), hypersup.end(), std::back_inserter(diff), 
    [refpoint] (element e) { return e.tau == refpoint.tau && e.x == refpoint.x && e.y != refpoint.y && e.eta == refpoint.eta; });
    return diff;
}


std::vector<element> surface::get_diff_eta(const std::vector<element> &hypersup, const element &refpoint)
{
    std::vector<element> diff;
    std::copy_if(hypersup.begin(), hypersup.end(), std::back_inserter(diff), 
    [refpoint] (element e) { return e.tau == refpoint.tau && e.x == refpoint.x && e.y == refpoint.y && e.eta != refpoint.eta; });
    return diff;
}

void element::print()
{
    std::cout << "Printing hypersurface element:" << std::endl
              << *this << std::endl;
}

double surface::element::t()
{
    return tau * cosh(eta);
}

double surface::element::z()
{
    return tau * sinh(eta);
}

double element::get_normal_size_sq()
{
    if (_normal_size == 0)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            _normal_size += dsigma[mu] * dsigma[mu] * utils::gmumu[mu];
        }
    }
    return _normal_size;
}

utils::four_vec surface::element::get_pos()
{
    return {tau, x, y, eta};
}

utils::four_vec surface::element::get_pos_mink()
{
    return {t(), x, y, z()};
}

double surface::element::distance(element second)
{
    return distance(*this, second);
}

double surface::element::distance(surface::element first, surface::element second)
{
    utils::four_vec vec = {second.t() - first.t(), second.x - first.x, second.y - first.y, second.z() - first.z()};
    return utils::get_norm_sq(vec);
}
