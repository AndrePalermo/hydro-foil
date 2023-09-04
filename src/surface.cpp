#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "surface.h"

using namespace std;

istream &operator>>(istream &stream, element &cell)
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

void read_hypersrface(string filename, vector<element> &hypersurface)
{
    ifstream input_file(filename);
    if (!input_file.is_open())
    {
        cout << "Failed to open " << filename << endl;
        exit(1);
    }
    cout << "Reading hypersurface from " << filename << endl;

    string line;
    while (getline(input_file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue; // skip empty or comment lines
        }

        std::istringstream iss(line);
        element cell;

        iss >> cell;

        // iss >> cell.tau >> cell.x >> cell.y >> cell.eta;
        // for (int mu = 0; mu < 4; mu++)
        // {
        //     iss >> cell.dsigma[mu];
        // }
        // for (int mu = 0; mu < 4; mu++)
        // {
        //     iss >> cell.u[mu];
        // }
        // iss >> cell.T >> cell.mub >> cell.muq >> cell.mus;
        // for (int i = 0; i < 4; i++)
        // {
        //     for (int j = 0; j < 4; j++)
        //     {
        //         iss >> cell.dbeta[i][j];
        //     }
        // }

        hypersurface.push_back(cell);
    }
    cout << "Reading successful!" << endl;
}

std::array<double, 4> get_mins(vector<element> &hypersup)
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

std::array<double, 4> get_maxs(vector<element> &hypersup)
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
    cout << "Printing hypersurface element:" << endl;
    cout << tau << "    " << x << "    " << y << "    " << eta << "    ";
    for (int mu = 0; mu < 4; mu++)
    {
        cout << dsigma[mu] << "    ";
    }
    for (int mu = 0; mu < 4; mu++)
    {
        cout << u[mu] << "    ";
    }
    cout << T << "    " << mub << "    " << muq << "    " << mus << "    ";
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << dbeta[i][j] << "    ";
        }
    }
    cout << endl;
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
