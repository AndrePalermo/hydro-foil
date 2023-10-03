#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "particle_data_group.h"

#define path_to_database "pdg_database/baryons_mesons.txt"

using namespace std;

pdg_particle::pdg_particle(int id_number){
    ifstream input_file(path_to_database);
    if(!input_file.is_open()){
        cout << "Failed to PDG database at " << path_to_database << endl;;
        exit(1);
    }
    
    int sign = 1; //used to choose between particle and antiparticle if pdg<0
    if(id_number < 0){
        id_number *= -1;
        sign = -1;
    }

    int _pdg_id,_q,_b,_s;
    float _spin;
    double _m; //GeV
    string _name;

    string line;
    while (getline(input_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // skip empty or comment lines
        }
        
        std::istringstream iss(line);
        
        iss >> _pdg_id >> _m >> _name >> _q >> _spin >> _b >> _s;
        if(_pdg_id == id_number){
            mass = _m;
            if(sign == -1){
                _name.insert(0,"anti-");
            }
            name = _name;
            pdg_ID = sign*_pdg_id;
            q = sign*_q;
            spin = _spin;
            b = sign*_b;
            s = sign*_s;
            return ;
        }
    }
    cout << "ERROR: No particle with PDG " << id_number << " is known. Check the id given or add it to the database at " << path_to_database << endl;
    exit(1);        
}

pdg_particle::~pdg_particle(){}

double pdg_particle::get_mass(){
    return mass;
}
int pdg_particle::get_id(){
    return pdg_ID;
}
int pdg_particle::get_q(){
    return q;
}
int pdg_particle::get_b(){
    return b;
}
int pdg_particle::get_s(){
    return s;
}
float pdg_particle::get_spin(){
    return spin;
}
std::string pdg_particle::get_name(){
    return name;
}

void pdg_particle::print(){
cout<<"Name: "<< name << ", Mass " << mass << ", Q_charge " << q << ", B_charge " << b << ", S_charge " << s << endl;
}

pdg_particle& pdg_particle::operator=(pdg_particle other){
    mass = other.get_mass();
    pdg_ID = other.get_id();
    q = other.get_q();
    b = other.get_b();
    s = other.get_s();
    name = other.get_name();

    return *this;
}

