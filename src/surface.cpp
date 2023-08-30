#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "surface.h"

using namespace std;

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

        iss >> cell.tau >> cell.x >> cell.y >> cell.eta;
        for (int mu = 0; mu < 4; mu++)
        {
            iss >> cell.dsigma[mu];
        }
        for (int mu = 0; mu < 4; mu++)
        {
            iss >> cell.u[mu];
        }
        iss >> cell.T >> cell.mub >> cell.muq >> cell.mus;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                iss >> cell.dbeta[i][j];
            }
        }

        hypersurface.push_back(cell);
    }
    cout << "Reading successful!" << endl;
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
