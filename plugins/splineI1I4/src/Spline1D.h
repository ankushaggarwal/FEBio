#include<iostream>
#include<fstream>
#include<vector>

class Spline1D {
    private:
        std::vector<double> knot_v;
        std::vector<double> coeff;
        int nKnots, degree, nNodes;
        
    public:
        Spline1D(){
            nKnots = degree = nNodes = 0;
        }
        void read_file(std::string fname);
        std::vector<double> cal_shape(double x);
        void set_degree(int d){degree = d;};
        void set_knot(std::vector<double> knot);
        int get_degree(){return degree;};
        int get_nNodes(){return nNodes;};
        const std::vector<double>& get_knot() const{return knot_v;};
        int find_interval(double x) const;
        double eval(double x);
        double deriv(double x);
        double second_deriv(double x);
        int eval(double x, std::vector<double> &N);
        int deriv(double x, std::vector<double> &dN);
        int second_deriv(double x, std::vector<double>& d2N);
};

