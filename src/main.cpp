#include "integrals.h"
#include "particle_data_group.h"
#include "utils.h"
#include "surface.h"

//TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

using namespace std;

int main(int argc, char** argv){

if(argc!=3){
    cout<< "INVALID SINTAX!"<<endl;
	cout<<"use './calc <surface_file> <output_file>' to compute Lambda polarization at decoupling."<<endl;
	exit(1);
}

string surface_file = argv[1];
string output_file = argv[2];

vector<element> hypersup = {};
read_hypersrface(surface_file, hypersup);

ofstream fout(output_file);
	 if (!fout) {
		 cout << "I/O error with " << output_file << endl;
		 exit(1);
	 }

int size_pt = 11;
int size_phi = 20;
int size_y = 20;
vector<double> pT = linspace(0,3.2,size_pt);
vector<double> phi =  linspace(0,2*PI,size_phi);
// vector<double> y_rap =  linspace(-1,1,size_y);

pdg_particle lambda(3122);
lambda.print();
// for(double iy : y_rap)
	for(double ipt : pT){
		for(double iphi : phi){
			polarization_projected(ipt, iphi, lambda, hypersup, fout);
		}
	}
cout<<"The calculation is done!"<<endl;
return 0;
}
