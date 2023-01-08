#include <iostream>
//#include "stdafx.h"
#include "FEElastinMaterial.h"

#ifdef WIN32

#define M_PI 3.14159265359
//the following two functions are defined in math_err_gamma.cpp
double erf(double);
double tgamma(double);
#endif // WIN32

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEElastinMaterial, FEElasticMaterial)
	ADD_PARAMETER(vol_fraction, "phi");
	ADD_PARAMETER(fiber_modulus, "kappaf");
	ADD_PARAMETER(mean_theta, "thetabar");
	ADD_PARAMETER(theta_std, "thetasigma");
	ADD_PARAMETER(aniso_fraction, "de");
END_FECORE_CLASS();
 
//-----------------------------------------------------------------------------
FEElastinMaterial::FEElastinMaterial(FEModel* pfem) : FEElasticMaterial(pfem){
	//define the integration rule
	nseg=3;
	order_quad=5;
	nint=nseg*order_quad;
	theta_iterate.resize(nint,0.);
	theta_weight_iterate.resize(nint,0.);
	
	update_integrationTheta(-M_PI/2.,M_PI/2.);
	ErrorF = erf(M_PI/(2.*sqrt(2)*theta_std));
}
//-----------------------------------------------------------------------------
bool FEElastinMaterial::Init()
{
	// Don't forget the base class initialization first.
	if (FEElasticMaterial::Init() == false) return false;
	
	//add conditions for ranges on the parameters
    return true;
}

//-----------------------------------------------------------------------------
mat3ds FEElastinMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	
	// calculate left Cauchy-Green tensor
	//mat3ds b = pt.LeftCauchyGreen();

	// calculate Green-Lagrange strain tensor
	mat3ds E = pt.Strain();

	mat3ds s(0.);

	// integrate over the theta
	for(int i=0;i<nint;i++){
		double thetai = theta_iterate[i]; //ith iteration point
		double weighti = theta_weight_iterate[i]; //ith integration point weight
		vec3d N = vec3d(cos(thetai),sin(thetai),0);
		double gammai = Gamma(thetai-mean_theta,theta_std,aniso_fraction);
		double Ef = N*(E*N);
		s += weighti*vol_fraction*fiber_modulus*gammai*Ef* dyad(N); 
	}

	return pt.push_forward(s);//1/J*F*s*F.transpose();
}

//-----------------------------------------------------------------------------
tens4ds FEElastinMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	
	// calculate left Cauchy-Green tensor
	//mat3ds b = pt.LeftCauchyGreen();

	// calculate Green-Lagrange strain tensor
	mat3ds E = pt.Strain();

	// integrate over the theta
	tens4ds D(0.);

	// integrate over the theta
	for(int i=0;i<nint;i++){
		double thetai = theta_iterate[i]; //ith iteration point
		double weighti = theta_weight_iterate[i]; //ith integration point weight
		vec3d N = vec3d(cos(thetai),sin(thetai),0);
		double gammai = Gamma(thetai-mean_theta,theta_std,aniso_fraction);
		//double Ef = N*(E*N);
		mat3ds dyadN = dyad(N);
		D += weighti*vol_fraction*fiber_modulus*gammai* dyad1s(dyadN); 
	}
	return pt.push_forward(D); //push forward
}

//-----------------------------------------------------------------------------
double FEElastinMaterial::Gamma(double theta, double sigma, double aniso_fraction){
	double normalize=ErrorF;
	return aniso_fraction*exp(-theta*theta/(2.*sigma*sigma))/normalize + (1-aniso_fraction)/M_PI;
}
void FEElastinMaterial::update_integrationTheta(double T1, double T2){
	double DeltaTheta=(T2-T1)/nseg;
	std::vector<double> temp_point(order_quad,0.), temp_weight(order_quad,0.);
	for(int i=0;i<nseg;i++){
		mapGauss(T1+i*DeltaTheta,T1+(i+1)*DeltaTheta,temp_point,temp_weight);
		for(int j=0;j<order_quad;j++){
			theta_iterate[i*order_quad+j]=temp_point[j];
			theta_weight_iterate[i*order_quad+j]=temp_weight[j];
		}
	}
}
void FEElastinMaterial::mapGauss(double x1, double x2, std::vector<double>& xp, std::vector<double>& wp){
	double jac=(x2-x1)/2.;
	double dsize=(x2+x1)/2.;
	if(order_quad==5){
		//in the parameteric coordinates
		double zetap[] = {
			-0.9061798459386640,
			-0.5384693101056831,
			0.0000000000000000,
			0.5384693101056831,
			0.9061798459386640};
		double weightp[] = {
			0.2369268850561891,
			0.4786286704993665,
			0.5688888888888889,
			0.4786286704993665,
			0.2369268850561891};	
		//transform them into the physical coordinates
		for(int i=0;i<order_quad;i++){
		xp[i]=jac*zetap[i]+dsize;
		wp[i]=jac*weightp[i];
		}
	}
	else
		std::cout<<"Order other than 5 not implemented yet"<<std::endl;
	return;
}
