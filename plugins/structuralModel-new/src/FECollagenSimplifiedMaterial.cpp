#include <iostream>
//#include "stdafx.h"
#include "FECollagenSimplifiedMaterial.h"

#ifdef WIN32

#define M_PI 3.14159265359
//the following two functions are defined in math_err_gamma.cpp
double erf(double);
double tgamma(double);
#endif // WIN32

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECollagenSimplifiedMaterial, FEElasticMaterial)
	ADD_PARAMETER(A, "A");
	ADD_PARAMETER(B, "B");
	ADD_PARAMETER(mean_theta, "thetabar");
	ADD_PARAMETER(theta_std, "thetasigma");
	ADD_PARAMETER(aniso_fraction, "de");
END_FECORE_CLASS();
 
//-----------------------------------------------------------------------------
FECollagenSimplifiedMaterial::FECollagenSimplifiedMaterial(FEModel* pfem) : FEElasticMaterial(pfem){
	//define the integration rule
	nseg=1;
	order_quad=10;
	nint=nseg*order_quad;
	theta_iterate.resize(nint,0.);
	theta_weight_iterate.resize(nint,0.);
	Gamma_iterate.resize(nint,0.);
	vec3d xaxis = vec3d(1.,0.,0.);
	nvector_iterate.resize(nint,xaxis);
	
	update_integrationTheta(-M_PI/2.,M_PI/2.);
	
	// calculate the Nvector_iterate
	for(int i=0;i<nint;i++){
		double thetai = theta_iterate[i]; 
		nvector_iterate[i] = vec3d(cos(thetai),sin(thetai),0);
	}
}
//-----------------------------------------------------------------------------
bool FECollagenSimplifiedMaterial::Init()
{
	// Don't forget the base class initialization first.
	if (FEElasticMaterial::Init() == false) return false;
	
	// evaluate the additional parameters needed for recruitment function which do not depend upon deformation
	double factor=sqrt(2.*M_PI)*theta_std;
	ErrorF = erf(M_PI/(2.*sqrt(2)*theta_std))*factor;
	//calculate the gamma function beforehand
	for(int i=0;i<nint;i++){
		Gamma_iterate[i] = Gamma(theta_iterate[i]-mean_theta,theta_std,aniso_fraction);
	}

	//add conditions for ranges on the parameters
    return true;
}

//-----------------------------------------------------------------------------
mat3ds FECollagenSimplifiedMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	// mat3d &F = pt.m_F;
	
	// calculate left Cauchy-Green tensor
	//mat3ds b = pt.LeftCauchyGreen();

	// get the element's local coordinate system
	mat3d Q = GetLocalCS(mp);

	// calculate Green-Lagrange strain tensor
	mat3ds E = pt.Strain();

	mat3ds s(0.);

	// integrate over the theta
	for(int i=0;i<nint;i++){
		//double thetai = theta_iterate[i]; //ith iteration point
		double weighti = theta_weight_iterate[i]; //ith integration point weight
		vec3d n = nvector_iterate[i]; //vec3d(cos(thetai),sin(thetai),0);
		// rotate the vector to material coordinate frame
		vec3d N = Q*n;
		double gammai = Gamma_iterate[i];
		double Ef = N*(E*N);
		s += weighti*A*gammai*(exp(B*Ef)-1)*dyad(N); 
	}

	return pt.push_forward(s);//1/J*F*s*F.transpose();
}

//-----------------------------------------------------------------------------
tens4ds FECollagenSimplifiedMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	// mat3d &F = pt.m_F;
	
	// calculate left Cauchy-Green tensor
	// mat3ds b = pt.LeftCauchyGreen();

	// get the element's local coordinate system
	mat3d Q = GetLocalCS(mp);

	// calculate Green-Lagrange strain tensor
	mat3ds E = pt.Strain();

	// integrate over the theta
	tens4ds D(0.);

	// integrate over the theta
	for(int i=0;i<nint;i++){
		//double thetai = theta_iterate[i]; //ith iteration point
		double weighti = theta_weight_iterate[i]; //ith integration point weight
		vec3d n = nvector_iterate[i]; //vec3d(cos(thetai),sin(thetai),0);
		// rotate the vector to material coordinate frame
		vec3d N = Q*n;
		double gammai = Gamma_iterate[i];
		double Ef = N*(E*N);
		mat3ds dyadN = dyad(N);
		D += weighti*A*gammai*B*exp(B*Ef)*dyad1s(dyadN); 
	}
	return pt.push_forward(D); //push forward
}

//-----------------------------------------------------------------------------
double FECollagenSimplifiedMaterial::Gamma(double theta, double sigma, double aniso_fraction){
	double normalize=ErrorF;
	return aniso_fraction*exp(-theta*theta/(2.*sigma*sigma))/normalize + (1-aniso_fraction)/M_PI;
}
//-----------------------------------------------------------------------------
void FECollagenSimplifiedMaterial::update_integrationTheta(double T1, double T2){
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
//-----------------------------------------------------------------------------
void FECollagenSimplifiedMaterial::mapGauss(double x1, double x2, std::vector<double>& xp, std::vector<double>& wp){
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
	else if(order_quad==9){
		//in the parameteric coordinates
		double zetap[] = {
			0.0000000000000000,
			-0.8360311073266358,
			0.8360311073266358,
			-0.9681602395076261,
			0.9681602395076261,
			-0.3242534234038089,
			0.3242534234038089,
			-0.6133714327005904,
			0.6133714327005904};
		double weightp[] = {
			0.3302393550012598,	
			0.1806481606948574,	
			0.1806481606948574,	
			0.0812743883615744,	
			0.0812743883615744,	
			0.3123470770400029,	
			0.3123470770400029,	
			0.2606106964029354,	
			0.2606106964029354};	
		//transform them into the physical coordinates
		for(int i=0;i<order_quad;i++){
		xp[i]=jac*zetap[i]+dsize;
		wp[i]=jac*weightp[i];
		}
	}
	else if(order_quad==10){
		//in the parameteric coordinates
		double zetap[] = {
			-0.9739065285171717,
			-0.8650633666889845,
			-0.6794095682990244,
			-0.4333953941292472,
			-0.1488743389816312,
			0.1488743389816312,
			0.4333953941292472,
			0.6794095682990244,
			0.8650633666889845,
			0.9739065285171717};
		double weightp[] = {
			0.0666713443086881,	
			0.1494513491505806,	
			0.2190863625159820,	
			0.2692667193099963,	
			0.2955242247147529,	
			0.2955242247147529,	
			0.2692667193099963,	
			0.2190863625159820,	
			0.1494513491505806,	
			0.0666713443086881};	
		//transform them into the physical coordinates
		for(int i=0;i<order_quad;i++){
		xp[i]=jac*zetap[i]+dsize;
		wp[i]=jac*weightp[i];
		}
	}
	else
		std::cout<<"Order other than 5 and 10 not implemented yet"<<std::endl;
	return;
}
