#include<iostream>
#include<fstream>
#include<vector>
#include "Spline1D.h"

class TensorProductSpline {
    private:
        std::vector<Spline1D> splines1D;
        int dim;
        std::vector<double> coeff;
        std::vector<double> point;
        int degree_x, degree_y;
        int nx, ny;
    public:
        TensorProductSpline(int d=2){
            dim = d;
            splines1D.resize(dim);
            point.resize(dim);
            degree_x = 0;
            degree_y = 0;
            nx = 0;
            ny = 0;
        }
        double eval(double x, double y);
        void deriv(double x, double y, double& dx, double& dy);
        void second_deriv(double x, double y, double& d2x, double& d2y, double& dxy);
        void set_splines(std::vector<Spline1D>& s1D);
        void read_file(std::string fname);
};
