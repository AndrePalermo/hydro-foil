#include "integrals.h"
#include "particle_data_group.h"
#include "utils.h"
#include "surface.h"

//TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

using namespace std;

int main(int argc, char** argv){

bool decay = false;

if(argc<3){
    cout<< "INVALID SINTAX!"<<endl;
	cout<<"use './foil <surface_file> <output_folder> <flags>' to compute Lambda polarization at decoupling."<<endl;
	exit(1);
}

if(argv[3]=="-D"s){
	decay = true;
	cout<<"Including calculations for the feed-down corrections!"<<endl;
}

string surface_file = argv[1];
string output_folder = argv[2];

vector<element> hypersup = {};
read_hypersrface(surface_file, hypersup);

string name_file_primary = output_folder+"/primary";

ofstream fout(name_file_primary);
	 if (!fout) {
		 cout << "I/O error with " << name_file_primary << endl;
		 exit(1);
	 }
int size_pt = 30;
int size_phi = 40;
int size_y = 20;
vector<double> pT = linspace(0,6.2,size_pt);
vector<double> phi =  linspace(0,2*PI,size_phi);
vector<double> y_rap =  linspace(-1,1,size_y);

pdg_particle Lambda(3122);
Lambda.print();
for(double ipt : pT){
	for(double iphi : phi){
		for(double iy : y_rap){
			polarization_exact_rapidity(ipt, iphi, iy, Lambda, hypersup, fout);
		}
	}
}


////////////////////////////////////DECAYS///////////////////////////////	 
if(decay){
	string table_file_sigma0 = output_folder + "/TableSigma0";
	string table_file_sigmastar = output_folder + "/TableSigmastar";
	ofstream fout_sigma0(table_file_sigma0);
		if (!fout_sigma0) {
			cout << "I/O error with " << table_file_sigma0 << endl;
			exit(1);
		}
	ofstream fout_sigmastar(table_file_sigmastar);
		if (!fout_sigmastar) {
			cout << "I/O error with " << table_file_sigmastar << endl;
			exit(1);
		}

	vector<double> pT_table = linspace(0,6.2,2*size_pt);
	vector<double> phi_table =  linspace(0,2*PI,2*size_phi);
	vector<double> y_rap_table =  linspace(-1,1,2*size_y);
	pdg_particle Sigma0(3212);
	pdg_particle SigmaStar(3224); //NB: there are three sigma* decaying to lambdas
	Sigma0.print();
	SigmaStar.print();
	for(double ipt : pT_table){
		for(double iphi : phi_table){
			for(double iy : y_rap_table){
				polarization_exact_rapidity(ipt, iphi, iy, Sigma0, hypersup, fout_sigma0);
				polarization_exact_rapidity(ipt, iphi, iy, SigmaStar, hypersup, fout_sigmastar);
			}
		}
	}

	string FD_file_sigma0 = output_folder + "/FeedDown_Sigma0";
	string FD_file_sigmastar = output_folder + "/FeedDown_Sigmastar";
	ofstream FDoutsigma0(FD_file_sigma0);
		if (!FDoutsigma0) {
			cout << "I/O error with " << FD_file_sigma0 << endl;
			exit(1);
		}
	ofstream FDoutsigmastar(FD_file_sigmastar);
		if (!FDoutsigmastar) {
			cout << "I/O error with " << FD_file_sigmastar << endl;
			exit(1);
		}

	for(double ipt : pT){
		for(double iphi : phi){
			for(double iy : y_rap){
				Lambda_polarization_FeedDown(ipt, iphi, iy, Sigma0, table_file_sigma0, FDoutsigma0);
				Lambda_polarization_FeedDown(ipt, iphi, iy, SigmaStar, table_file_sigmastar, FDoutsigmastar);
			}
		}
	}
}

cout<<"The calculation is done!"<<endl;
return 0;
}
