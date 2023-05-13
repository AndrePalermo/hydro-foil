#ifndef INTEGRALS_H
#define INTEGRALS_H

#include <vector>
#include <fstream>

#include "surface.h"
#include "utils.h"
#include "particle_data_group.h"

using namespace std;


//computes the mean spin vector for particles emitted at midrapidity for the particle "particle" using
//the freezeout file "freezeout_sup" and writing the output in "fileout" as a table:
//pt phi denominator numerator_varpi numerator_xi
//refer to 2009.13449
void polarization_midrapidity(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

void shear_correction_to_spectra(double pT, double phi, double Dt, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

#endif