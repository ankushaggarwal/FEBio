#include "TensorProductSpline.h"

void TensorProductSpline::read_file(std::string fname){
    std::ifstream indata(fname, std::ios::in | std::ios::binary); // opens the file
    
    if(!indata) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    indata.read((char*)&degree_x,sizeof(int));
    int kx;
    indata.read((char*)&kx,sizeof(int));
    std::vector<double> knot_x(kx,0);
    for(int i=0;i<kx;i++)
        indata.read((char*)&knot_x[i],sizeof(double));
    
    indata.read((char*)&degree_y,sizeof(int));
    int ky;
    indata.read((char*)&ky,sizeof(int));
    std::vector<double> knot_y(ky,0);
    for(int i=0;i<ky;i++)
        indata.read((char*)&knot_y[i],sizeof(double));

    int N;
    indata.read((char*)&N,sizeof(int));
    coeff.resize(N);
    for(int i=0;i<N;i++)
        indata.read((char*)&coeff[i],sizeof(double));

    indata.close();

    splines1D[0].set_degree(degree_x);
    splines1D[1].set_degree(degree_y);
    splines1D[0].set_knot(knot_x);
    splines1D[1].set_knot(knot_y);
    nx = splines1D[0].get_nNodes();
    ny = splines1D[1].get_nNodes();
    if (nx*ny != N)
        std::cerr << "Number of points in the input spline is not consistent" << std::endl;
    std::cout << "Read the spline from file " << fname << std::endl;
    /*
    std::cout << degree_x << "\t" << kx << std::endl;
    for(int i=0; i<kx; i++)
        std::cout << knot_x[i] << std::endl;
    std::cout << degree_y << "\t" << ky << std::endl;
    for(int i=0; i<ky; i++)
        std::cout << knot_y[i] << std::endl;
    std::cout << N << "\t" << N << std::endl;
    for(int i=0; i<N; i++)
        std::cout << coeff[i] << std::endl;
    */
}

double TensorProductSpline::eval(double x, double y){
    std::vector<double> Nx(splines1D[0].get_degree());
    std::vector<double> Ny(splines1D[1].get_degree());
    const int i_x = splines1D[0].eval(x,Nx);
    const int j_y = splines1D[1].eval(y,Ny);
    double val=0.;
    for(int i=0; i<=degree_x; i++){
        for(int j=0; j<=degree_y; j++){
            val += Nx[i]*Ny[j]*coeff[(i+i_x)*ny+(j+j_y)];
        }
    }
    return val;
}

void TensorProductSpline::deriv(double x, double y, double& dx, double& dy){
    std::vector<double>  Nx(degree_x+1,0.);
    std::vector<double> dNx(degree_x+1,0.);
    std::vector<double>  Ny(degree_y+1,0.);
    std::vector<double> dNy(degree_y+1,0.);
    const int i_x = splines1D[0].eval(x,Nx);
    const int j_y = splines1D[1].eval(y,Ny);
    splines1D[0].deriv(x,dNx);
    splines1D[1].deriv(y,dNy);
    dx=0.;dy=0.;
    for(int i=0; i<=degree_x; i++){
        for(int j=0; j<=degree_y; j++){
            dx += dNx[i]* Ny[j]*coeff[(i+i_x)*ny+(j+j_y)];
            dy +=  Nx[i]*dNy[j]*coeff[(i+i_x)*ny+(j+j_y)];
        }
    }
    return;
}

void TensorProductSpline::second_deriv(double x, double y, double& d2x, double& d2y, double& dxy){
    std::vector<double>   Nx(degree_x+1,0.);
    std::vector<double>  dNx(degree_x+1,0.);
    std::vector<double> d2Nx(degree_x+1,0.);
    std::vector<double>   Ny(degree_y+1,0.);
    std::vector<double>  dNy(degree_y+1,0.);
    std::vector<double> d2Ny(degree_y+1,0.);
    const int i_x = splines1D[0].eval(x,Nx);
    const int j_y = splines1D[1].eval(y,Ny);
    splines1D[0].deriv(x,dNx);
    splines1D[1].deriv(y,dNy);
    splines1D[0].second_deriv(x,d2Nx);
    splines1D[1].second_deriv(y,d2Ny);
    d2x=0.; d2y=0.; dxy=0.;
    for(int i=0; i<=degree_x; i++){
        for(int j=0; j<=degree_y; j++){
            d2x += d2Nx[i]* Ny[j]*coeff[(i+i_x)*ny+(j+j_y)];
            d2y +=  Nx[i]*d2Ny[j]*coeff[(i+i_x)*ny+(j+j_y)];
            dxy += dNx[i]* dNy[j]*coeff[(i+i_x)*ny+(j+j_y)];
        }
    }
    return;
}

void TensorProductSpline::set_splines(std::vector<Spline1D>& s1D){
    if (dim != s1D.size())
        std::cerr << "Size of input 1D splines is not the same as the dimension of TensorProductSpline set when initializing" << std::endl;
    for(int i=0;i<dim;i++){
        splines1D[i].set_degree(s1D[i].get_degree());
        splines1D[i].set_knot(s1D[i].get_knot());
    }
    degree_x = splines1D[0].get_degree();
    degree_y = splines1D[1].get_degree();
    nx = splines1D[0].get_nNodes();
    ny = splines1D[1].get_nNodes();
    coeff.resize(nx*ny);
}
