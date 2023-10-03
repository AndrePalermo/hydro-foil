#ifndef PARTICLE_DATA_GROUP_H //if it has not been defined already executes the if-endif block
#define PARTICLE_DATA_GROUP_H
#include <string>

class pdg_particle{
    private:
        int pdg_ID;
        double mass;
        std::string name;
        float spin;
        //charges
        int q; //electric charge
        int b; //baryon number
        int s; //strangeness
    public:
    pdg_particle(int id_number);
    ~pdg_particle();
    
    double get_mass();
    int get_id();
    std::string get_name();
    int get_q();
    int get_b();
    int get_s();
    float get_spin();
    void print();
    
    pdg_particle& operator=(pdg_particle other);
};

#endif

