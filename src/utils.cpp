#include<vector>
#include<tuple>
#include<iostream>
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

std::tuple<double, double > interpolator::adjacent_points(double x, std::vector<double> direction){
    //"x" is some point that we want to locate in the vector "direction"
    //gives for x in the interval [x0,x1] returns the tuple (x0,x1).

    //xi = xmin+Deltax*i where Deltax = xmax-xmin/(size(x)-1)
    // i = (xi-xmin)/Deltax= (size(x)-1)*(xi-xmin)/(xmax-xmin)
    double xmin = direction[0];
    double xmax = direction[std::size(direction)-1];
    double i = (std::size(direction)-1)*(x-xmin)/(xmax-xmin);
    int idx = floor(i);
    if(i>std::size(direction)-1 || i < 0){
        std::cout<<"Out of bounds! ERROR"<<std::endl;
        exit(1);
    }
    if(abs(i-idx)<1e-8){
        return {idx,idx};
    }
    return {idx,idx+1};

}

double interpolator::trilinear_interpolation(double x, double y, double z){
    //algorithm taken from https://en.wikipedia.org/wiki/Trilinear_interpolation
    auto [x0,x1,boolx] = adjacent_points(x,x_);
    auto [y0,y1,booly] = adjacent_points(y,y_);
    auto [z0,z1,boolz] = adjacent_points(z,z_);

    double c000 = f[x0+size(x_)*y0+size(x_)*size(y_)*z0];
    double c001 = f[x0+size(x_)*y0+size(x_)*size(y_)*z1];
    double c010 = f[x0+size(x_)*y1+size(x_)*size(y_)*z0];
    double c100 = f[x1+size(x_)*y0+size(x_)*size(y_)*z0];
    double c011 = f[x0+size(x_)*y1+size(x_)*size(y_)*z1];
    double c110 = f[x1+size(x_)*y1+size(x_)*size(y_)*z0];
    double c101 = f[x1+size(x_)*y0+size(x_)*size(y_)*z1];
    double c111 = f[x1+size(x_)*y1+size(x_)*size(y_)*z1];

    // std::cout<<c000;
    double dx = x0 - x1+1e-8;
    double dy = y0 - y1+1e-8;
    double dz = z0 - z1+1e-8;
    double denominator = dx * dy * dz;

    double a0 = (-c000 * x1 * y1 * z1 + c001 * x1 * y1 * z0 + c010 * x1 * y0 * z1 - c011 * x1 * y0 * z0) / denominator +
                (c100 * x0 * y1 * z1 - c101 * x0 * y1 * z0 - c110 * x0 * y0 * z1 + c111 * x0 * y0 * z0) / denominator;
    double a1 = (c000 * y1 * z1 - c001 * y1 * z0 - c010 * y0 * z1 + c011 * y0 * z0) / denominator +
                (-c100 * y1 * z1 + c101 * y1 * z0 + c110 * y0 * z1 - c111 * y0 * z0) / denominator;
    double a2 = (c000 * x1 * z1 - c001 * x1 * z0 - c010 * x1 * z1 + c011 * x1 * z0) / denominator +
                (-c100 * x0 * z1 + c101 * x0 * z0 + c110 * x0 * z1 - c111 * x0 * z0) / denominator;
    double a3 = (c000 * x1 * y1 - c001 * x1 * y1 - c010 * x1 * y0 + c011 * x1 * y0) / denominator +
                (-c100 * x0 * y1 + c101 * x0 * y1 + c110 * x0 * y0 - c111 * x0 * y0) / denominator;
    double a4 = (-c000 * z1 + c001 * z0 + c010 * z1 - c011 * z0 + c100 * z1 - c101 * z0 - c110 * z1 + c111 * z0) / denominator;
    double a5 = (-c000 * y1 + c001 * y1 + c010 * y0 - c011 * y0 + c100 * y1 - c101 * y1 - c110 * y0 + c111 * y0) / denominator;
    double a6 = (-c000 * x1 + c001 * x1 + c010 * x1 - c011 * x1 + c100 * x0 - c101 * x0 - c110 * x0 + c111 * x0) / denominator;
    double a7 = (c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) / denominator;

    return a0 + a1*x + a2*y + a3*z + a4*x*y + a5*x*z + a6*y*z + a7*x*y*z;
}


