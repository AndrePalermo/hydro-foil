#include "integrals.h"
#include "particle_data_group.h"
#include "utils.h"
#include "surface.h"

#define CALC_PROJ true
// TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		std::cout << "INVALID SINTAX!"
				  << std::endl
				  << "use './foil <surface_file> <output_file>' to compute Lambda polarization at decoupling."
				  << std::endl;
		exit(1);
	}

	std::string surface_file = argv[1];
	std::string output_file = argv[2];

	std::vector<element> hypersup = {};
	read_hypersrface(surface_file, hypersup);

	auto min_points = get_mins(hypersup);
	auto max_points = get_maxs(hypersup);

	std::cout << "tau range: [" << min_points[0] << "," << max_points[0] << "]\n"
			  << "x range: [" << min_points[1] << "," << max_points[1] << "]\n"
			  << "y range: [" << min_points[2] << "," << max_points[2] << "]\n"
			  << "eta range: [" << min_points[3] << "," << max_points[3] << "]"
			  << std::endl;

	int timelikes = 0;
	for (auto el : hypersup)
	{
		if (el.is_timelike())
		{
			timelikes++;
		}
	}

	std::cout << 100 * timelikes / hypersup.size() << "% timelike elements." << std::endl;

	// assume a reference point x0
	// for x0 - n * d to x0 + n * d
	// {
	// 	b1 = (taylor expand u(x+d)) / T(x0)
	// 	b2 = (taylor expand b(x+d))
	// 	write x+d error(b1, b) error(b2, b)
	// }

	// 	std::ofstream fout(output_file);

	// #if CALC_PROJ
	// 	std::ofstream fout_proj(output_file + "_projected");
	// #endif

	// 	if (!fout)
	// 	{
	// 		std::cout << "I/O error with " << output_file << std::endl;
	// 		exit(1);
	// 	}

	// #if CALC_PROJ
	// 	if (!fout_proj)
	// 	{
	// 		std::cout << "I/O error with " << output_file + "_projected" << std::endl;
	// 		exit(1);
	// 	}
	// #endif

	// 	int size_pt = 31;
	// 	int size_phi = 40;
	// 	int size_y = 20;
	// 	std::vector<double> pT = linspace(0, 6.2, size_pt);
	// 	std::vector<double> phi = linspace(0, 2 * PI, size_phi);
	// 	// std::vector<double> y_rap =  linspace(-1,1,size_y);

	// 	pdg_particle part(3122);
	// 	part.print();
	// 	// for(double iy : y_rap)
	// 	for (double ipt : pT)
	// 	{
	// 		for (double iphi : phi)
	// 		{
	// 			polarization_midrapidity_linear(ipt, iphi, part, hypersup, fout);
	// #if CALC_PROJ
	// 			polarization_projected(ipt, iphi, part, hypersup, fout_proj);
	// #endif
	// 		}
	// 	}
	// 	std::cout << "The calculation is done!" << std::endl;

	// 	fout.close();

	// #if CALC_PROJ
	// 	fout_proj.close();
	// #endif

	return 0;
}
