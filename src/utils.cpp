#include<vector>
#include<tuple>
#include<iostream>
#include <algorithm> 
#include <fstream>
#include <sstream>
#include <string>
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


std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> file_to_4_column(std::string filename, int f_column){
    std::ifstream input_file(filename);
    if(!input_file.is_open()){
        std::cout << "Failed to open " << filename << std::endl;
        exit(1);
    }
    std::vector<double> x_vector,y_vector,z_vector,fxyz_vector;
    std::string line;
    while (getline(input_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // skip empty or comment lines
        }
        double x,y,z,fxyz;
        std::istringstream iss(line);

        iss>>x>>y>>z;
        int i=3;
        double dummy = 0;
        while(i<f_column){
            iss>>dummy;
            i++;
        }
        iss>>fxyz;

        x_vector.push_back(x);
        y_vector.push_back(y);
        z_vector.push_back(z);
        fxyz_vector.push_back(fxyz);
        
    }

    std::vector<double>::iterator it;
    std::sort(x_vector.begin(), x_vector.end());
    it = std::unique (x_vector.begin(), x_vector.end());            
    x_vector.resize( std::distance(x_vector.begin(),it) );

    std::sort(y_vector.begin(), y_vector.end());
    it = std::unique (y_vector.begin(), y_vector.end());            
    y_vector.resize( std::distance(y_vector.begin(),it) );
    
    std::sort(z_vector.begin(), z_vector.end());
    it = std::unique (z_vector.begin(), z_vector.end());            
    z_vector.resize( std::distance(z_vector.begin(),it) );

    if (x_vector.size()*y_vector.size()*z_vector.size()!=fxyz_vector.size()){
        std::cout << "Size mismatch :("<<std::endl; 
    }

    return {x_vector,y_vector,z_vector,fxyz_vector};
}

std::vector<double> linspace(double min, double max, int size){
    //returns a vector of doubles from min to max (included), equally spaced with spacing max-min/size
    std::vector<double> result;
    double interval = (max-min)/size;
    double tmp_min = min;
    for( int i=0;i<=size;i++){
        result.push_back(tmp_min);
        tmp_min+=interval;
    }
    return result;
}

interpolator::interpolator(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> f):
    x_(x), y_(y), z_(z), f_(f) {}

interpolator::interpolator(std::string filename, int f_column){
    std::ifstream input_file(filename);
    if(!input_file.is_open()){
        std::cout << "Failed to open " << filename << std::endl;
        exit(1);
    }
    std::vector<double> x_vector,y_vector,z_vector,fxyz_vector;
    std::string line;
    while (getline(input_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // skip empty or comment lines
        }
        double x,y,z,fxyz;
        std::istringstream iss(line);

        iss>>x>>y>>z;
        int i=3;
        double dummy = 0;
        while(i<f_column){
            iss>>dummy;
            i++;
        }
        iss>>fxyz;

        x_vector.push_back(x);
        y_vector.push_back(y);
        z_vector.push_back(z);
        fxyz_vector.push_back(fxyz);
        
    }
    
    std::vector<double>::iterator it;
    std::sort(x_vector.begin(), x_vector.end());
    it = std::unique (x_vector.begin(), x_vector.end());            
    x_vector.resize( std::distance(x_vector.begin(),it) );

    std::sort(y_vector.begin(), y_vector.end());
    it = std::unique (y_vector.begin(), y_vector.end());            
    y_vector.resize( std::distance(y_vector.begin(),it) );
    
    std::sort(z_vector.begin(), z_vector.end());
    it = std::unique (z_vector.begin(), z_vector.end());            
    z_vector.resize( std::distance(z_vector.begin(),it) );

    if (x_vector.size()*y_vector.size()*z_vector.size()!=fxyz_vector.size()){
        std::cout << "Size mismatch :("<<std::endl; 
        exit(1);
    }

    x_ = x_vector;
    y_ = y_vector;
    z_ = z_vector;
    f_ = fxyz_vector;
}

interpolator::~interpolator(){}

double interpolator::trilinear_interpolation(double x, double y, double z){
    const int  nx_grid = x_.size();
    const int  ny_grid = y_.size();
    const int  nz_grid = z_.size();
    const double xmaxG= x_[nx_grid-1];
    const double xminG= x_[0];

    const double ymaxG= y_[ny_grid-1];
    const double yminG= y_[0];

    const double zmaxG= z_[ny_grid-1];
    const double zminG= z_[0];


    const double dxG = (xmaxG - xminG) / (nx_grid - 1);
    const double dyG = (ymaxG - yminG) / (ny_grid - 1);
    const double dzG = (zmaxG - zminG) / (nz_grid - 1);
    int ix = (int)((x - xminG) / dxG);
    int iy = (int)((y - yminG) / dyG);
    int iz = (int)((z - zminG) / dzG);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix > nx_grid - 2) ix = nx_grid - 2;
    if (iy > ny_grid - 2) iy = ny_grid - 2;
    if (iz > nz_grid - 2) iz = nz_grid - 2;
    const double xm = x - xminG - ix * dxG;
    const double ym = y - yminG - iy * dyG;
    const double zm = z - zminG - iz * dzG;
    double wx[2] = {1. - xm / dxG, xm / dxG};
    double wy[2] = {1. - ym / dyG, ym / dyG};
    double wz[2] = {1. - zm / dzG, zm / dzG};
    double return_val = 0.;
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++){
    
       return_val += wx[jx] * wy[jy] * wz[jz] * f_[ix + jx +nx_grid*(iy + jy +ny_grid*(iz + jz))];
     }

    return return_val;
}


