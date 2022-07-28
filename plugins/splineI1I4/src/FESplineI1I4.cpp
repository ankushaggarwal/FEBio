/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "FESplineI1I4.h"

// define the material parameters
BEGIN_FECORE_CLASS(FESplineI1I4, FEUncoupledMaterial)
	ADD_PARAMETER(m_rho, "rho");
	ADD_PARAMETER(m_mu , "mu");
	ADD_PARAMETER(m_fiber, "fiber");
    ADD_PARAMETER(inp_fname,"filename");
	
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESplineI1I4::FESplineI1I4(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_rho = 1.0;
	m_mu = 0.0;
	m_fiber = vec3d(1,0,0);
    inp_fname = "spline.dat";
}

//-----------------------------------------------------------------------------
// This (optional) function is called during initialization. This can be done
// to do one-time initialization. Since for this material this is not required,
// this function could have been ommitted. 
// Make sure to always call the base class (usually first). 
bool FESplineI1I4::Init()
{
	// Don't forget the base class initialization first.
	if (FEUncoupledMaterial::Init() == false) return false;
    spline.read_file(inp_fname);
	return true;
}

//-----------------------------------------------------------------------------
mat3ds FESplineI1I4::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*m_fiber(mp); a0.unit();
	vec3d a = F * a0;
	double lam = Jm13 * a.unit();
	mat3ds m = dyad(a);

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double J1 = B.tr();
	double J4 = lam * lam;

	double W1;
	double W4;
    spline.deriv(J1,J4,W1,W4);
    W1 += m_mu;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B * W1 + m*(W4*J4);

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = m_rho*T.dev() * (2.0 / J);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FESplineI1I4::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*m_fiber(mp); a0.unit();
	vec3d a = F * a0;
	double lam = Jm13 * a.unit();
	mat3ds m = dyad(a);

	// Invariants of B (= invariants of C)
	double J1 = B.tr();
	double J4 = lam * lam;

	// put strain energy derivatives here
	// Wi = dW/dIi
	double W1, W4;
    spline.deriv(J1,J4,W1,W4);
	double W11, W44, W14;
    spline.second_deriv(J1,J4,W11,W44,W14);
    W1 += m_mu;
	// ------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B * W1 + m * (W4 * J4);

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds devs = T.dev() * (2.0 / J);

	// calculate dWdC:C
	double WC = W1 * J1 + W4 * J4;

	// calculate C:d2WdCdC:C
	double CWWC = W11 * J1 * J1 + 2.0 * W14 * J1 * J4 + W44 * J4 * J4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4 = dyad4s(B);
	tens4ds Bxm = dyad1s(B, m);
	tens4ds mxm = dyad1s(m);

	// d2W/dCdC:C
	mat3ds WCCxC = B * (W11 * J1) + B*(W14*J4) + m*(W14*J4*J1) + m*(W44*J4*J4);

	tens4ds Cw = BxB * W11 + Bxm * (W14 * J4) + mxm * (W44 * J4 * J4);

	tens4ds cw = Cw*(4.0 *Ji) - dyad1s(WCCxC, I) * (4.0 / 3.0 * Ji) + IxI * (4.0 / 9.0 * Ji * CWWC);
	tens4ds c = dyad1s(devs, I) * (-2.0 / 3.0) + (I4 - IxI / 3.0) * (4.0 / 3.0 * Ji * WC) + cw;

	return c*m_rho;
}

//-----------------------------------------------------------------------------
double FESplineI1I4::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*m_fiber(mp); a0.unit();
	vec3d a = F * a0;
	double lam = Jm13 * a.unit();

	// Invariants of B (= invariants of C)
	double J1 = B.tr();
	double J4 = lam * lam;

	// calculate sed
	double sed = spline.eval(J1,J4)+m_mu*(J1-3); 

	return sed*m_rho;
}
