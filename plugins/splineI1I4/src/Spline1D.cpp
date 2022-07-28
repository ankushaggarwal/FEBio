#include "Spline1D.h"

void Spline1D::read_file(std::string fname){
    std::ifstream indata;
    indata.open(fname); // opens the file
    if(!indata) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    indata >> degree;
    indata >> nKnots;
    knot_v.resize(nKnots);
    for(int i=0;i<nKnots;i++)
       indata >> knot_v[i];

    indata >> nNodes;
    coeff.resize(nNodes);
    for(int i=0;i<nNodes;i++)
       indata >> coeff[i];

    indata.close();
    std::cout << "File read." << std::endl;
    std::cout << "Degree is " << degree << std::endl;
    std::cout << "Number of control points is " << nNodes << std::endl;
    std::cout << "Size of the knot vector is " << nKnots << std::endl;
    if(nKnots != nNodes + degree + 1)
        std::cerr << "Size of knot vector, coefficients, and degree are not consistent" << std::endl;
}

std::vector<double> Spline1D::cal_shape(double x){
    //An inefficient implmenetation to calculate shape function and its first and second derivatives with the entire length of the shape function
    std::vector<double> shape_f(nNodes,0);
    std::vector<double> dN(nNodes,0);
    std::vector<double> dN_temp(nNodes,0);
    std::vector<double> d2N(nNodes,0);
    const int in_k=find_interval(x);
    shape_f[in_k] = 1.;
    for(int d=1;d<degree+1;d++){
        shape_f[in_k-d] = (knot_v[in_k+1]-x)/(knot_v[in_k+1]-knot_v[in_k-d+1])*shape_f[in_k-d+1];
        for(int j=in_k-d+1;j<in_k;j++)
            shape_f[j]  = (x-knot_v[j])/(knot_v[j+d]-knot_v[j])*shape_f[j] 
                        + (knot_v[j+d+1]-x)/(knot_v[j+d+1]-knot_v[j+1])*shape_f[j+1];
        shape_f[in_k] = (x-knot_v[in_k])/(knot_v[in_k+d]-knot_v[in_k])*shape_f[in_k];

        if(d==degree-1){
            //derivatives use shape functions from previous degree, i.e. d-1
            dN[in_k-degree] = -degree*shape_f[in_k-degree+1]/(knot_v[in_k+1]-knot_v[in_k-degree+1]); //left end special case
            for(int j=in_k-degree+1;j<in_k;j++)
                dN[j] = degree*(shape_f[j]/(knot_v[j+degree]-knot_v[j]) - shape_f[j+1]/(knot_v[j+1+degree]-knot_v[j+1]));
            dN[in_k] = degree*shape_f[in_k]/(knot_v[in_k+degree]-knot_v[in_k]); //right end special case
        }

        if(d==degree-2){
            //calculation of second derivative
            //first calculate first derivatives (of degree d-1)
            dN_temp[in_k-degree+1] = -(degree-1)*shape_f[in_k-(degree-1)+1]/(knot_v[in_k+1]-knot_v[in_k-(degree-1)+1]); //left end special case
            for(int j=in_k-(degree-1)+1;j<in_k;j++)
                dN_temp[j] = (degree-1)*(shape_f[j]/(knot_v[j+degree-1]-knot_v[j]) - shape_f[j+1]/(knot_v[j+1+(degree-1)]-knot_v[j+1]));
            dN_temp[in_k] = (degree-1)*shape_f[in_k]/(knot_v[in_k+degree-1]-knot_v[in_k]); //right end special case
           
            //now calculate the second derivative
            d2N[in_k-degree] = -degree*dN_temp[in_k-degree+1]/(knot_v[in_k+1]-knot_v[in_k-degree+1]); //left end special case
            for(int j=in_k-degree+1; j<in_k; j++)
                d2N[j] = degree*(dN_temp[j]/(knot_v[j+degree]-knot_v[j]) - dN_temp[j+1]/(knot_v[j+1+degree]-knot_v[j+1]));
            d2N[in_k] = degree*dN_temp[in_k]/(knot_v[in_k+degree]-knot_v[in_k]); //right end special case
        }
    }
    return shape_f;
}

double Spline1D::eval(double x){
    std::vector<double> shape_f(degree+1,0);
    const int in_k=eval(x,shape_f);
    double val=0.;
    for(int i=0; i<=degree; i++)
        val += shape_f[i]*coeff[in_k+i];
    return val;
}

int Spline1D::eval(double x, std::vector<double> &shape_f){
    const int in_k=find_interval(x);
    shape_f[in_k-in_k+degree] = 1.;
    for(int d=1;d<degree+1;d++){
        shape_f[degree-d] = (knot_v[in_k+1]-x)/(knot_v[in_k+1]-knot_v[in_k-d+1])*shape_f[degree-d+1]; //left end special case
        for(int j=in_k-d+1;j<in_k;j++)
            shape_f[j-in_k+degree]  = (x-knot_v[j])/(knot_v[j+d]-knot_v[j])*shape_f[j-in_k+degree] 
                        + (knot_v[j+d+1]-x)/(knot_v[j+d+1]-knot_v[j+1])*shape_f[j+1-in_k+degree];
        shape_f[degree] = (x-knot_v[in_k])/(knot_v[in_k+d]-knot_v[in_k])*shape_f[degree]; //right end special case
    }
    return in_k-degree>0?in_k-degree:0;
}

double Spline1D::deriv(double x){
    std::vector<double> dN(degree+1,0);
    const int in_k = deriv(x,dN);
    double val=0.;
    for(int i=0; i<=degree; i++){
        //std::cout << i << "\t" << dN[i] << std::endl;
        val += dN[i]*coeff[in_k+i];
    }
    return val;
}

int Spline1D::deriv(double x, std::vector<double> &dN){
    std::vector<double> shape_f(degree+1,0);
    const int in_k=find_interval(x);
    shape_f[in_k-in_k+degree] = 1.;
    for(int d=1;d<degree;d++){
        shape_f[degree-d] = (knot_v[in_k+1]-x)/(knot_v[in_k+1]-knot_v[in_k-d+1])*shape_f[degree-d+1]; //left end special case
        for(int j=in_k-d+1;j<in_k;j++)
            shape_f[j-in_k+degree]  = (x-knot_v[j])/(knot_v[j+d]-knot_v[j])*shape_f[j-in_k+degree] 
                        + (knot_v[j+d+1]-x)/(knot_v[j+d+1]-knot_v[j+1])*shape_f[j+1-in_k+degree];
        shape_f[degree] = (x-knot_v[in_k])/(knot_v[in_k+d]-knot_v[in_k])*shape_f[degree]; //right end special case
    }

    //derivatives use shape functions from previous degree, i.e. d-1
    dN[0] = -degree*shape_f[1]/(knot_v[in_k+1]-knot_v[in_k-degree+1]); //left end special case
    for(int j=in_k-degree+1;j<in_k;j++)
        dN[j-in_k+degree] = degree*(shape_f[j-in_k+degree]/(knot_v[j+degree]-knot_v[j]) - shape_f[j+1-in_k+degree]/(knot_v[j+1+degree]-knot_v[j+1]));
    dN[degree] = degree*shape_f[degree]/(knot_v[in_k+degree]-knot_v[in_k]); //right end special case

    return in_k-degree>0?in_k-degree:0;
}

double Spline1D::second_deriv(double x){
    std::vector<double> d2N(degree+1,0);
    const int in_k = second_deriv(x,d2N);
    double val=0.;
    for(int i=0; i<=degree; i++){
        //std::cout << i << "\t" << d2N[i] << std::endl;
        val += d2N[i]*coeff[in_k+i];
    }
    return val;
}

int Spline1D::second_deriv(double x, std::vector<double>& d2N){
    std::vector<double> shape_f(degree+1,0);
    std::vector<double> dN_temp(degree+1,0);
    const int in_k=find_interval(x);
    shape_f[degree] = 1.;
    for(int d=1;d<degree-1;d++){
        shape_f[degree-d] = (knot_v[in_k+1]-x)/(knot_v[in_k+1]-knot_v[in_k-d+1])*shape_f[degree-d+1]; //left end special case
        for(int j=in_k-d+1;j<in_k;j++)
            shape_f[j-in_k+degree]  = (x-knot_v[j])/(knot_v[j+d]-knot_v[j])*shape_f[j-in_k+degree] 
                        + (knot_v[j+d+1]-x)/(knot_v[j+d+1]-knot_v[j+1])*shape_f[j+1-in_k+degree];
        shape_f[degree] = (x-knot_v[in_k])/(knot_v[in_k+d]-knot_v[in_k])*shape_f[degree]; //right end special case
    }

    //calculation of second derivative
    //first calculate first derivatives (of degree d-1)
    dN_temp[1] = -(degree-1)*shape_f[2]/(knot_v[in_k+1]-knot_v[in_k-(degree-1)+1]); //left end special case
    for(int j=in_k-(degree-1)+1;j<in_k;j++)
        dN_temp[j-in_k+degree] = (degree-1)*(shape_f[j-in_k+degree]/(knot_v[j+degree-1]-knot_v[j]) - shape_f[j+1-in_k+degree]/(knot_v[j+1+(degree-1)]-knot_v[j+1]));
    dN_temp[degree] = (degree-1)*shape_f[degree]/(knot_v[in_k+degree-1]-knot_v[in_k]); //right end special case
   
    //now calculate the second derivative
    d2N[0] = -degree*dN_temp[1]/(knot_v[in_k+1]-knot_v[in_k-degree+1]); //left end special case
    for(int j=in_k-degree+1; j<in_k; j++)
        d2N[j-in_k+degree] = degree*(dN_temp[j-in_k+degree]/(knot_v[j+degree]-knot_v[j]) - dN_temp[j+1-in_k+degree]/(knot_v[j+1+degree]-knot_v[j+1]));
    d2N[degree] = degree*dN_temp[degree]/(knot_v[in_k+degree]-knot_v[in_k]); //right end special case

    return in_k-degree>0?in_k-degree:0; 
}

void Spline1D::set_knot(std::vector<double> knot){
    knot_v = knot;
    nKnots = knot_v.size();
    if (degree==0)
        std::cerr << "Degree has not been set. Set degree before setting the knot vector " << std::endl;
    nNodes = nKnots - degree - 1;
    return;
}

int Spline1D::find_interval(double x) const{
    if (x<knot_v[0]) return degree; //extrapolate to the left
    int in_k;
    for(int i=0;i<nKnots-2;i++){
        if(x<knot_v[i+1]){
           // std::cout << "Found the interval " << i << " " << x << " " << knot_v[i+1] << std::endl;
           return i;
           in_k = i;
        }
    }
    //if (x>=knot_v[nKnots-1]) 
    return nKnots-degree-2; //extrapolates to the right
}
