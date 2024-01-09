#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "surface.h"
#include "utils.h"


using namespace std;

void read_hypersrface(string filename, vector<element> &hypersurface){
    ifstream input_file(filename);
    if(!input_file.is_open()){
        cout << "Failed to open " << filename << endl;
        exit(1);
    }
    cout << "Reading hypersurface from " << filename << endl;

    string line;
    while (getline(input_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // skip empty or comment lines
        }
        
        std::istringstream iss(line);
        element cell;
        
        iss >> cell.tau >> cell.x >> cell.y >> cell.eta;
        for(int mu=0;mu<4;mu++){
            iss >> cell.dsigma[mu];
        }        
        for(int mu=0;mu<4;mu++){
            iss >> cell.u[mu];
        }
        iss >> cell.T >> cell.mub >> cell.muq >> cell.mus;
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                iss >> cell.dbeta[i][j];
            }
        }
        
        hypersurface.push_back(cell);
    }
    cout << "Reading successful!" << endl;
}

void element::print(){
    cout<<"Printing hypersurface element:"<<endl;
    cout << tau << "    " << x << "    " << y << "    " << eta << "    ";
        for(int mu=0;mu<4;mu++){
            cout << dsigma[mu] << "    ";
        }        
        for(int mu=0;mu<4;mu++){
            cout << u[mu] << "    ";
        }
        cout << T << "    " << mub << "    " << muq << "    " << mus << "    ";
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                cout << dbeta[i][j] << "    ";
            }
        }
        cout << endl;
}

element new_dbeta(element surf_old, int tag){
    std::cout<<"REMEMBER: this function should only be used in the isothermal freeze-out case!"<<std::endl;
    
    element new_elem = surf_old; 

    double theta_over_T=0;
    std::array<double,4> u_up={0}, u_low = {0}, acc_over_T_shear_low={0},
                         acc_over_T_vort_low={0};
    std::array<std::array<double,4>,4> vorticity_low={0}, fullshear_low={0}, deltamunu_low={0},
                        sigma_over_T={0},ang_vel_part_of_vorticity_over_T={0};
    
    for(int mu=0; mu<4; mu++){
        u_up[mu] = surf_old.u[mu];
        u_low[mu] = surf_old.u[mu]*gmumu[mu];   
        for(int nu=0;nu<4;nu++){
            vorticity_low[mu][nu]=0.5*( surf_old.dbeta[nu][mu] - surf_old.dbeta[mu][nu]); // dbeta[mu][nu] is partial_mu\beta_nu!
            fullshear_low[mu][nu]=0.5*( surf_old.dbeta[nu][mu] + surf_old.dbeta[mu][nu]);
            }
    }

    for(int mu=0;mu<4;mu++){
		for(int nu=0;nu<4;nu++){
			acc_over_T_shear_low[mu] += 2*fullshear_low[mu][nu]*u_up[nu];
			acc_over_T_vort_low[mu] += 2*vorticity_low[mu][nu]*u_up[nu];
            if(mu==nu){
                deltamunu_low[mu][nu] = gmumu[mu]-u_low[mu]*u_low[nu];
            }
            else{
                deltamunu_low[mu][nu] = -u_low[mu]*u_low[nu];
            }	
		}
	}
    if(acc_over_T_shear_low[0]!=acc_over_T_vort_low[0]){
        std::cout<<"Acceleration is different if computed from shear and from vorticity! Error!"<<std::endl;
        exit(1);
        }

    for(int mu=0; mu<4; mu++){
			theta_over_T += fullshear_low[mu][mu]*gmumu[mu];
    }
    
    for(int mu=0;mu<4;mu++){
		for(int nu=0;nu<4;nu++){
            sigma_over_T[mu][nu] = fullshear_low[mu][nu]
                    -0.5*(acc_over_T_shear_low[mu]*u_low[nu]+acc_over_T_shear_low[nu]*u_low[mu])
                    -theta_over_T/3 * deltamunu_low[mu][nu];
        }
    }

    for(int mu=0;mu<4;mu++){
		for(int nu=0;nu<4;nu++){
            ang_vel_part_of_vorticity_over_T[mu][nu] = vorticity_low[mu][nu]
                    -0.5*(acc_over_T_vort_low[mu]*u_low[nu]-acc_over_T_vort_low[nu]*u_low[mu]);
        }
    }
    
    switch(tag){
        case 0:
        //acceleration vorticity
        for(int mu=0;mu<4;mu++){
		    for(int nu=0;nu<4;nu++){
                new_elem.dbeta[mu][nu] = 0.5*(acc_over_T_vort_low[mu]*u_low[nu]-acc_over_T_vort_low[nu]*u_low[mu]);
            }
        }
        break;

        case 1:
        //acceleration thermal shear
        for(int mu=0;mu<4;mu++){
		    for(int nu=0;nu<4;nu++){
                new_elem.dbeta[mu][nu] = 0.5*(acc_over_T_shear_low[mu]*u_low[nu]+acc_over_T_shear_low[nu]*u_low[mu]);            
            }
        }
        break;

        case 2:
        //angular velocity vorticity
        for(int mu=0;mu<4;mu++){
		    for(int nu=0;nu<4;nu++){
                new_elem.dbeta[mu][nu] = ang_vel_part_of_vorticity_over_T[mu][nu];
            }
        }
        break;

        case 3:
        //proper shear
        for(int mu=0;mu<4;mu++){
		    for(int nu=0;nu<4;nu++){
                new_elem.dbeta[mu][nu] = sigma_over_T[mu][nu];
            }
        }
        break;
        
        case 4:
        //expansion thermal shear
        for(int mu=0;mu<4;mu++){
		    for(int nu=0;nu<4;nu++){
                new_elem.dbeta[mu][nu] = theta_over_T/3 * deltamunu_low[mu][nu];
            }
        }
        break;

        default:
        std::cout<< "Tag Unkown! Error!"<<std::endl;
        exit(1);
    }
    return new_elem;
}

std::array<vector<element>,5> components_freeze_out(vector<element> &freeze_out_sup){
    std::array<vector<element>,5> components = {};
    for(element cell : freeze_out_sup){
        components[0].push_back(new_dbeta(cell,0)); //acceleration vorticity
        components[2].push_back(new_dbeta(cell,1)); //acceleration shear
        components[3].push_back(new_dbeta(cell,2)); //angular velocity
        components[4].push_back(new_dbeta(cell,3)); //proper shear
        components[5].push_back(new_dbeta(cell,4)); //expansion
    }
    
    return components;
}



