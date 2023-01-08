#pragma once
#include <vector>
#include "FEBioMech/FEElasticMaterial.h"

//
//! Elastin Material

//! Implementation of a structural model for the elastin fibers.
class FEElastinMaterial : public FEElasticMaterial
{
public:
	FEElastinMaterial(FEModel* pfem);

public:
	double vol_fraction;      //!< Volume fraction
	double fiber_modulus;     //!< Modulus for a single fiber
	double mean_theta;        //!< Preferred direction of the elastin distribution
	double theta_std;         //!< Standard deviation of the elastin distribution
	double aniso_fraction;    //!< Fraction of the anisotropic elastin

public:
	//! calculation the stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calcualte the tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! data initialization and checking
	bool Init() override;

private:
	//! for calculating the Gamma function - the fiber angular distribution 
	double Gamma(double theta, double sigma, double aniso_fraction);
	void update_integrationTheta(double T1, double T2);
	void mapGauss(double x1, double x2, std::vector<double>& xp, std::vector<double>& wp);

	//! private variables
	int nint, nseg, order_quad;
	std::vector<double> theta_iterate, theta_weight_iterate;
	double ErrorF; //For normalizing the Gaussian distribution correctly

protected:

	// declare parameter list
	DECLARE_FECORE_CLASS();
};
