#include "integrals.h"
#include "particle_data_group.h"
#include "utils.h"
#include "surface.h"
#include <filesystem>

//TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

using namespace std;

int main(int argc, char** argv){

if(argc!=3){
    cout<< "INVALID SINTAX!"<<endl;
	cout<<"use './foil <surface_file> <output_folder>' to compute Lambda polarization at decoupling."<<endl;
	exit(1);
}

string surface_file = argv[1];
string output_folder = argv[2];
filesystem::create_directories("./"+output_folder);

std::array<string,5> components_names={output_folder+"/vort_acc",
									   output_folder+"/shear_acc",
									   output_folder+"/ang_vel",
									   output_folder+"/proper_shear",
									   output_folder+"/expansion"};

for(auto names : components_names){
	filesystem::create_directories("./"+names);
}

vector<element> hypersup = {};
read_hypersrface(surface_file, hypersup);

std::array<vector<element>,5> comp_FO = components_freeze_out(hypersup);



	pdg_particle Lambda(3122);
	Lambda.print();
	polarization_components(Lambda, comp_FO, components_names);
	
cout<<"The calculation is done!"<<endl;
return 0;
}
