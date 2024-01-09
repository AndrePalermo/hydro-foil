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

vector<element> hypersup = {};
read_hypersrface(surface_file, hypersup);

std::array<vector<element>,5> comp_FO = components_freeze_out(hypersup);

// for(int iii=0;iii<5;iii++){
// 	cout<<"check "<<iii<<endl;
// 	cout<<comp_FO[iii].size()<<" "<<hypersup.size()<<endl;
// for(element cell : comp_FO[iii]){
// 	for(int mu=0; mu<4;mu++){
// 		for(int nu=0; nu<4;nu++){
// 			if(std::isnan(cell.dbeta[mu][nu])){
// 				cout<<"NAN ecnounter!"<<endl;
// 			}
// 		}	
// 	}
// }
// }

std::array<string,5> components_names={output_folder+"/vort_acc",
									   output_folder+"/shear_acc",
									   output_folder+"/ang_vel",
									   output_folder+"/proper_shear",
									   output_folder+"/expansion"};
	

	pdg_particle Lambda(3122);
	Lambda.print();
	polarization_components(Lambda, comp_FO, components_names);
	
cout<<"The calculation is done!"<<endl;
return 0;
}
