#include "integrals.h"
#include "particle_data_group.h"
#include "utils.h"
#include "surface.h"
#include<tuple>

//TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

using namespace std;

int main(int argc, char** argv){
	std::string filename = "./polarization_primary";

// int size_pt = 101;
// int size_phi = 120;
// int size_y = 120;
// auto [pT,phi,y_rap,table] = file_to_4_column(filename,5);// =  linspace(-1,1,size_y);
// int total_size = size(pT)*size(phi)*size(y_rap);


// vector<double> table(total_size,0);
// // vector<double> tablept(total_size,0);
// // vector<double> tablephi(total_size,0);
// // vector<double> tabley(total_size,0);

// int count = 0; 
// for(double iy : y_rap){
// 	for(double iphi : phi){
// 		for(double ipt : pT){
// 			table[count] = (ipt*iy*sin(iphi));
// 			// tablept[count] = ipt;
// 			// tablephi[count] = iphi;
// 			// tabley[count] = iy;;
// 			count++;
// 		}
// 	}
// }

interpolator pippo(filename,5);
int size_pttest = 40;
int size_phitest = 31;
int size_ytest = 30;
vector<double> pTtest = linspace(0,3.113,size_pttest);
vector<double> phitest =  linspace(0.1,2*PI,size_phitest);
vector<double> y_raptest =  linspace(0.3,0.76,size_ytest);
double errorinf = 0;
double errorsq = 0;
for(double iy : y_raptest){
	for(double iphi : phitest){
		for(double ipt : pTtest){
			// auto [x0,x1] = pippo.adjacent_points(ipt,pT);
			double error = pippo.trilinear_interpolation(ipt,iphi,iy);
			// cout<<"differenza fra funzioni: "<<pippo.trilinear_interpolation(ipt,iphi,iy)<<"  "<<(ipt*iy*sin(iphi))<<endl;
			// cout << error << endl; 
			errorinf += fabs(error);
			errorsq += error*error;

		}
	}
}
cout<<"errorinf: "<<errorinf<<"    errorsq: "<<sqrt(errorsq)<<endl;

cout<<"The calculation is done!"<<endl;
return 0;
}
