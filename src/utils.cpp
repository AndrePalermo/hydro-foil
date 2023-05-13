#include<vector>

#include "utils.h"
#include "particle_data_group.h"


int levi(int i, int j, int k, int l){
// Levi-Civita symbols
// i,j,k,l = 0...3 i.e. upper indices
 if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) return 0;
 else return ( (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12 );
}

std::vector<double> linspace(double min, double max, int size){
    std::vector<double> result;
    double interval = (max-min)/size;
    double tmp_min = min;
    for( int i=0;i<size;i++){
        result.push_back(tmp_min);
        tmp_min+=interval;
    }
    return result;
}