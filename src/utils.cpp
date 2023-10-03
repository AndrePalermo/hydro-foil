#include<vector>
#include<tuple>
#include<iostream>
#include <algorithm> 
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
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
    x_(x), y_(y), z_(z), f_(f), nx(x_.size()), ny(y_.size()), nz(z_.size()),
    xmaxG(x_[nx-1]),ymaxG(y_[ny-1]),zmaxG(z_[nz-1]),
    xminG(x_[0]),yminG(y_[0]),zminG(z_[0]),
    dxG((xmaxG - xminG) / (nx - 1)),dyG((ymaxG - yminG) / (ny - 1)),dzG((zmaxG - zminG) / (nz - 1)) {}

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
    nx = x_.size();
    ny = y_.size();
    nz = z_.size();
    xmaxG = x_[nx -1];
    ymaxG = y_[ny -1];
    zmaxG = z_[nz -1];
    xminG = x_[0];
    yminG = y_[0];
    zminG = z_[0];
    dxG = (xmaxG - xminG) / (nx - 1);
    dyG = (ymaxG - yminG) / (ny - 1);
    dzG = (zmaxG - zminG) / (nz - 1);
}

interpolator::~interpolator(){}

double interpolator::trilinear_interpolation(double x, double y, double z){
    const double xmaxG= x_[nx-1];
    const double xminG= x_[0];

    const double ymaxG= y_[ny-1];
    const double yminG= y_[0];

    const double zmaxG= z_[nz-1];
    const double zminG= z_[0];


    const double dxG = (xmaxG - xminG) / (nx - 1);
    const double dyG = (ymaxG - yminG) / (ny - 1);
    const double dzG = (zmaxG - zminG) / (nz - 1);
    int ix = (int)((x - xminG) / dxG);
    int iy = (int)((y - yminG) / dyG);
    int iz = (int)((z - zminG) / dzG);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix > nx - 2) ix = nx - 2;
    if (iy > ny - 2) iy = ny - 2;
    if (iz > nz - 2) iz = nz - 2;
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
    
       return_val += wx[jx] * wy[jy] * wz[jz] * f_index(ix+jx,iy+jy,iz+jx);//f_[ix + jx +nx_grid*(iy + jy +ny_grid*(iz + jz))];
     }

    return return_val;
}

double interpolator::f_index(int ix, int iy, int iz){
    return f_[iz+nz*iy+nz*ny*ix];
}

double interpolator::trilinear_interpol(double x, double y, double z){    
    int ix = (int)((x - xminG) / dxG);
    int iy = (int)((y - yminG) / dyG);
    int iz = (int)((z - zminG) / dzG);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix > nx - 2) ix = nx - 2;
    if (iy > ny - 2) iy = ny - 2;
    if (iz > nz - 2) iz = nz - 2;
    
    const double f_x00 = std::lerp(f_index(ix,iy,iz),f_index(ix+1,iy,iz),(x-x_[ix])/dxG);
    const double f_x10 = std::lerp(f_index(ix,iy+1,iz),f_index(ix+1,iy+1,iz),(x-x_[ix])/dxG);
    const double f_xy0 =  std::lerp(f_x00,f_x10,(y-y_[iy])/dyG);
    const double f_x01 = std::lerp(f_index(ix,iy,iz+1),f_index(ix+1,iy,iz+1),(x-x_[ix])/dxG);
    const double f_x11 = std::lerp(f_index(ix,iy+1,iz+1),f_index(ix+1,iy+1,iz+1),(x-x_[ix])/dxG);
    const double f_xy1 =  std::lerp(f_x01,f_x11,(y-y_[iy])/dyG);
    const double f_xyz = std::lerp(f_xy0,f_xy1,(z-z_[iz])/dzG);

    return f_xyz; 
}



