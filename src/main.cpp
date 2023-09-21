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
string output_string = argv[2];

vector<element> hypersup = {};
read_hypersrface(surface_file, hypersup);

string table_file_sigma0 = output_string + "_TableSigma0";
string table_file_sigmastar = output_string + "_TableSigmastar";
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


int size_pt = 30;
int size_phi = 60;
int size_y = 40;
vector<double> pTtable = linspace(0,6.2,size_pt);
vector<double> phitable =  linspace(0,2*PI,size_phi);
vector<double> y_raptable =  linspace(-2,2,size_y);

pdg_particle Sigma0(3212);
pdg_particle SigmaStar(3224); //NB: there are three sigma* decaying to lambdas
Sigma0.print();
SigmaStar.print();
for(double ipt : pTtable){
	for(double iphi : phitable){
		for(double iy : y_raptable){
			polarization_exact_rapidity(ipt, iphi, iy, Sigma0, hypersup, fout_sigma0);
			polarization_exact_rapidity(ipt, iphi, iy, SigmaStar, hypersup, fout_sigmastar);
		}
	}
}

string FD_file_sigma0 = output_string + "_FeedDown_Sigma0";
string FD_file_sigmastar = output_string + "_FeedDown_Sigmastar";
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

size_pt = 30;
size_phi = 60;
size_y = 20;
vector<double> pT = linspace(0,6.2,size_pt);
vector<double> phi =  linspace(0,2*PI,size_phi);
vector<double> y_rap =  linspace(-1,1,size_y);
for(double ipt : pTtable){
	for(double iphi : phitable){
		for(double iy : y_raptable){
			Lambda_polarization_FeedDown(ipt, iphi, iy, Sigma0, table_file_sigma0, FDoutsigma0);
			Lambda_polarization_FeedDown(ipt, iphi, iy, SigmaStar, table_file_sigmastar, FDoutsigmastar);
		}
	}
}

cout<<"The calculation is done!"<<endl;
return 0;
}
