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
                                    })->tau;

    auto min_x = std::min_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.x < second.x;
                                  })->x;

    auto min_y = std::min_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.y < second.y;
                                  })->y;

    auto min_eta = std::min_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.eta < second.eta;
                                    })->eta;
     return {min_tau, min_x, min_y, min_eta};
}


std::array<double, 4> surface::get_maxs(const std::vector<element> &hypersup)
{
     auto max_tau = std::max_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.tau < second.tau;
                                    })->tau;

    auto max_x = std::max_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.x < second.x;
                                  })->x;

    auto max_y = std::max_element(hypersup.begin(),
                                  hypersup.end(),
                                  [](const element &first, const element &second)
                                  {
                                      return first.y < second.y;
                                  })->y;

    auto max_eta = std::max_element(hypersup.begin(),
                                    hypersup.end(),
                                    [](const element &first, const element &second)
                                    {
                                        return first.eta < second.eta;
                                    })->eta;
    return {max_tau, max_x, max_y, max_eta};
}

void element::print()
{
    std::cout << "Printing hypersurface element:" << std::endl
    << this << std::endl;
}

double element::get_normal_size()
{
    if (_normal_size == 0)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            _normal_size += dsigma[mu] * dsigma[mu] * gmumu[mu];
        }
    }
    return _normal_size;
}


