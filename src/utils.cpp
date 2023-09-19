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


std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> file_to_4_column(std::string filename)
{
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

        iss>>x>>y>>z>>fxyz;

        x_vector.push_back(x);
        y_vector.push_back(y);
        z_vector.push_back(z);
        fxyz_vector.push_back(fxyz);
        
    }

    std::vector<double>::iterator it;
    it = std::unique (x_vector.begin(), x_vector.end());            
    x_vector.resize( std::distance(x_vector.begin(),it) );

    it = std::unique (y_vector.begin(), y_vector.end());            
    y_vector.resize( std::distance(y_vector.begin(),it) );

    it = std::unique (z_vector.begin(), z_vector.end());            
    z_vector.resize( std::distance(z_vector.begin(),it) );

    if (x_vector.size()*y_vector.size()*z_vector.size()!=fxyz_vector.size()){
        std::cout << "Size miss_match "<<std::endl; 
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

interpolator::~interpolator(){}

std::tuple<int, int > interpolator::adjacent_points(double x, std::vector<double> direction){
    //"x" is some point that we want to locate in the vector "direction"
    //gives for x in the interval [x0,x1] returns the tuple (x0,x1).

    //xi = xmin+Deltax*i where Deltax = xmax-xmin/(size(x)-1)
    // i = (xi-xmin)/Deltax= (size(x)-1)*(xi-xmin)/(xmax-xmin)
    double xmin = direction[0];
    double xmax = direction[std::size(direction)-1];
    double deltax = (xmax-xmin)/(std::size(direction)-1.0);
    double i = (x-xmin)/deltax;
    int idx = floor(i);
    if(idx>std::size(direction)-1 || idx < 0){
        std::cout<<"Out of bounds! ERROR"<<std::endl;
        exit(1);
    }

    if(fabs(i-idx)<1e-8){
        return {idx,idx};
    }
    return {idx,idx+1};

}

double interpolator::trilinear_interpolation(double x, double y, double z){
    const int  nx_grid = x_.size();
    const int  ny_grid = y_.size();
    const int  nz_grid = z_.size();
    const double xmaxG= x_[nx_grid-1];
    const double xminG= x_[0];

    const double ymaxG= y_[nx_grid-1];
    const double yminG= y_[0];

    const double zmaxG= z_[nx_grid-1];
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
    double weta[2] = {1. - zm / dzG, zm / dzG};
    double return_val = 0.;
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++){
    
       return_val += wx[jx] * wy[jy] * weta[jz] * f_[ix + jx +nx_grid*(iy + jy +ny_grid*(iz + jz))];
     }

    return return_val;
}
//     //algorithm taken from https://en.wikipedia.org/wiki/Trilinear_interpolation
//     auto [ix0,ix1] = adjacent_points(x,x_);
//     auto [iy0,iy1] = adjacent_points(y,y_);
//     auto [iz0,iz1] = adjacent_points(z,z_);
    
//     double c000 = f_[ix0+size(x_)*iy0+size(x_)*size(y_)*iz0];
//     double c001 = f_[ix0+size(x_)*iy0+size(x_)*size(y_)*iz1];
//     double c010 = f_[ix0+size(x_)*iy1+size(x_)*size(y_)*iz0];
//     double c100 = f_[ix1+size(x_)*iy0+size(x_)*size(y_)*iz0];
//     double c011 = f_[ix0+size(x_)*iy1+size(x_)*size(y_)*iz1];
//     double c110 = f_[ix1+size(x_)*iy1+size(x_)*size(y_)*iz0];
//     double c101 = f_[ix1+size(x_)*iy0+size(x_)*size(y_)*iz1];
//     double c111 = f_[ix1+size(x_)*iy1+size(x_)*size(y_)*iz1];

//     // std::cout<<c000;
//     double x0 = x_[ix0];
//     double x1 = x_[ix1];
//     double y0 = x_[iy0];
//     double y1 = x_[iy1];
//     double z0 = x_[iz0];
//     double z1 = x_[iz1];
//     double dx = x_[ix0]-x_[ix1]+1e-8;
//     double dy = y_[iy0]-y_[iy1]+1e-8;
//     double dz = z_[iz0]-z_[iz1]+1e-8;
//     double denominator = dx * dy * dz;

//     double a0 = (-c000 * x1 * y1 * z1 + c001 * x1 * y1 * z0 + c010 * x1 * y0 * z1 - c011 * x1 * y0 * z0) / denominator +
//                 (c100 * x0 * y1 * z1 - c101 * x0 * y1 * z0 - c110 * x0 * y0 * z1 + c111 * x0 * y0 * z0) / denominator;
//     double a1 = (c000 * y1 * z1 - c001 * y1 * z0 - c010 * y0 * z1 + c011 * y0 * z0) / denominator +
//                 (-c100 * y1 * z1 + c101 * y1 * z0 + c110 * y0 * z1 - c111 * y0 * z0) / denominator;
//     double a2 = (c000 * x1 * z1 - c001 * x1 * z0 - c010 * x1 * z1 + c011 * x1 * z0) / denominator +
//                 (-c100 * x0 * z1 + c101 * x0 * z0 + c110 * x0 * z1 - c111 * x0 * z0) / denominator;
//     double a3 = (c000 * x1 * y1 - c001 * x1 * y1 - c010 * x1 * y0 + c011 * x1 * y0) / denominator +
//                 (-c100 * x0 * y1 + c101 * x0 * y1 + c110 * x0 * y0 - c111 * x0 * y0) / denominator;
//     double a4 = (-c000 * z1 + c001 * z0 + c010 * z1 - c011 * z0 + c100 * z1 - c101 * z0 - c110 * z1 + c111 * z0) / denominator;
//     double a5 = (-c000 * y1 + c001 * y1 + c010 * y0 - c011 * y0 + c100 * y1 - c101 * y1 - c110 * y0 + c111 * y0) / denominator;
//     double a6 = (-c000 * x1 + c001 * x1 + c010 * x1 - c011 * x1 + c100 * x0 - c101 * x0 - c110 * x0 + c111 * x0) / denominator;
//     double a7 = (c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) / denominator;
//     std::cout<<"cose: "<<c000<<" "<<c001<<" "<<c010<<" "<<c100<<std::endl;
//     std::cout<<"cose2: "<<c011<<" "<<c110<<" "<<c101<<" "<<c111<<std::endl;
//     std::cout<<"coseA: "<<a0<<" "<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<std::endl;
//     std::cout<<"result:"<<a0 + a1*x + a2*y + a3*z + a4*x*y + a5*x*z + a6*y*z + a7*x*y*z<<std::endl;
//     return a0 + a1*x + a2*y + a3*z + a4*x*y + a5*x*z + a6*y*z + a7*x*y*z;
}


