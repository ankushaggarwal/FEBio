#pragma once
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEUncoupledFiberExpLinear.h>
#include "TensorProductSpline.h"

//-----------------------------------------------------------------------------
// Constitutive formulation from:
// Holzapfel, et.a., "Determination of layer-specific mechanical properties of human
// coronary arteries with nonatherosclerotic intimal thickening
// and related constitutive modeling", Am J Physiol Heart Circ Physiol 289
class FESplineI1I4 : public FEUncoupledMaterial
{
public:
	FESplineI1I4(FEModel* pfem);

public:
	double		m_rho;
	double		m_mu;
    std::string inp_fname;
	FEParamVec3 m_fiber;
    TensorProductSpline spline;

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;

	// This (optional) function can be used to do one-time material initialization.
	virtual bool Init() override;

protected:

	// declare parameter list
	DECLARE_FECORE_CLASS();
};
