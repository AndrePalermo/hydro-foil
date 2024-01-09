#ifndef INTEGRALS_H
#define INTEGRALS_H

#include <vector>
#include <fstream>

#include "surface.h"
#include "utils.h"
#include "particle_data_group.h"

using namespace std;


//computes the mean spin vector for particles emitted at midrapidity for the particle "particle" using
//the freezeout file "freezeout_sup". For the thermal vorticity, eq. 64 of 2304.02276 is used
//The output is written in "fileout" as a table:
//pt phi denominator numerator_varpi numerator_xi
//refer to 2103.10917
void polarization_midrapidity(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

//same as the previous function but uses the linear approximation for the vorticity induced polarization. This function is faster.
void polarization_midrapidity_linear(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

//Similar to "polarization_midrapidity_linear" but computes polarization projecting the gradients 
//tangent to the hypersurface at point x before integrating  
void polarization_projected(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

//Same as "polarization_midrapidity", but the table now includes the rapidity "y":
//pt phi y denominator numerator_varpi numerator_xi
void polarization_exact_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

//auxiliary function to compute polarization for any spin
double aux_exact_polarization(double spin, double pu, double T, double mutot, double theta_sq);

//auxiliary function used for distribution functions
int statistics(double spin);

//integrates the thermal spectrum of "particle"
void spectrum_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);

void Lambda_polarization_FeedDown(double pT, double phi, double y_rap, pdg_particle mother, 
    interpolator &spectrum_interpolator, array<interpolator,4> &S_vorticity_interpolator, array<interpolator,4> &S_shear_interpolator, ofstream &fileout);

void polarization_components(double pT, double phi, double y_rap, pdg_particle particle, 
        std::array<vector<element>,5> components, std::array<ofstream,5> fileout_list);

void Lambda_FeedDown_nointerpolation(double pT, double phi, double y_rap, pdg_particle mother, 
    vector<element> &freeze_out_sup, ofstream &fileout);

#endif