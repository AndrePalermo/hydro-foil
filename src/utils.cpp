#include<vector>
#include <iostream>
#include "utils.h"
#include "particle_data_group.h"


int levi(int i, int j, int k, int l){
// Levi-Civita symbols
// i,j,k,l = 0...3 i.e. upper indices
 if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) return 0;
 else return ( (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12 );
}

int g(int mu, int nu){
    if(mu > 3 || mu < 0 || nu > 3 || nu < 0){
        std::cout<<"error with the indices of the metric"<<std::endl;
        exit(1);
    }
    if(mu==nu){
        if(mu==0){
            return 1;
        }
        else{
            return -1;
        }
    }
    return 0;
}

std::vector<double> linspace(double min, double max, int size){
    std::vector<double> result;
    double interval = (max-min)/size;
    double tmp_min = min;
    for( int i=0;i<=size;i++){
        result.push_back(tmp_min);
        tmp_min+=interval;
    }
    return result;
}

