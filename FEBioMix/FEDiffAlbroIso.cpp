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



#include "stdafx.h"
#include "FEDiffAlbroIso.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include <FEBioFluid/FEMultiphasicFSI.h>
#include <FEBioFluid/FEFluidSolutes.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDiffAlbroIso, FESoluteDiffusivity)
	ADD_PARAMETER(m_diff0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "free_diff");
	ADD_PARAMETER(m_cdinv , FE_RANGE_GREATER_OR_EQUAL(0.0), "cdinv"    );
	ADD_PARAMETER(m_alphad, FE_RANGE_GREATER_OR_EQUAL(0.0), "alphad"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDiffAlbroIso::FEDiffAlbroIso(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_diff0 = 1;
	m_cdinv = m_alphad = 0;
    m_lsol = -1;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEDiffAlbroIso::Init()
{
	if (FESoluteDiffusivity::Init() == false) return false;

	// get the grandparent material which must be
    // a biphasic-solute/triphasic/multiphasic material
    FESolute* pSol = dynamic_cast<FESolute*> (GetParent());
    m_lsol = pSol->GetSoluteLocalID();
    
	if (m_lsol == -1) {
		feLogError("Invalid value for sol"); 
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffAlbroIso::Free_Diffusivity(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint* spt = mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint* fspt = mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FEMultiphasicFSIMaterialPoint* mfpt = mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
	
    // solute concentration
    double ca = 0;
    if (spt)
        ca = spt->m_ca[m_lsol];
    else if (fspt)
        ca = fspt->m_ca[m_lsol];
    else if (mfpt)
        ca = mfpt->m_ca[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_cdinv*ca);
    
	return d;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffAlbroIso::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
    FESolutesMaterialPoint* spt = mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint* fspt = mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FEMultiphasicFSIMaterialPoint* mfpt = mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
	
    // solute concentration
    double ca, c, dkdc, k;
    ca = c = dkdc = k = 0;
    if (spt)
    {
        ca = spt->m_ca[m_lsol];
        c = spt->m_c[m_lsol];
        dkdc = spt->m_dkdc[m_lsol][isol];
        k = spt->m_k[m_lsol];
    }
    else if (fspt)
    {
        ca = fspt->m_ca[m_lsol];
        c = fspt->m_c[m_lsol];
        dkdc = fspt->m_dkdc[m_lsol][isol];
        k = fspt->m_k[m_lsol];
    }
    else if (mfpt)
    {
        ca = mfpt->m_ca[m_lsol];
        c = mfpt->m_c[m_lsol];
        dkdc = mfpt->m_dkdc[m_lsol][isol];
        k = mfpt->m_k[m_lsol];
    }
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_cdinv*ca);
    // derivative of d w.r.t. actual concentration
    double dc = -m_cdinv*d;
    
    
    // tangent w.r.t. concentration
    if (isol == m_lsol)
        return dc*(k+dkdc*c);
    else
        return dc*dkdc*c;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor.
mat3ds FEDiffAlbroIso::Diffusivity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint* et = mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint* ppt = mp.ExtractData<FEBiphasicMaterialPoint>();
    FEBiphasicFSIMaterialPoint* bfpt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
	FESolutesMaterialPoint* spt = mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint* fspt = mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FEMultiphasicFSIMaterialPoint* mfpt = mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
	
	// relative volume
	double J = et->m_J;
	
	// solid volume fraction in reference configuration
	double phi0 = 0;
    if (ppt)
        phi0 = ppt->m_phi0t;
    else if (bfpt)
        phi0 = bfpt->m_phi0;
    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = 0;
    if (spt)
        ca = spt->m_ca[m_lsol];
    else if (fspt)
        ca = fspt->m_ca[m_lsol];
    else if (mfpt)
        ca = mfpt->m_ca[m_lsol];
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
	
	// diffusivity tensor
    mat3dd dt(d);
	
	return dt;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4dmm FEDiffAlbroIso::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint* et = mp.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint* ppt = mp.ExtractData<FEBiphasicMaterialPoint>();
    FEBiphasicFSIMaterialPoint* bfpt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    FESolutesMaterialPoint* spt = mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint* fspt = mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FEMultiphasicFSIMaterialPoint* mfpt = mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// relative volume
	double J = et->m_J;
	
	// solid volume fraction in reference configuration
    double phi0 = 0;
    if (ppt)
        phi0 = ppt->m_phi0t;
    else if (bfpt)
        phi0 = bfpt->m_phi0;
    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca = 0;
    double c = 0;
    double dkdJ = 0;
    if (spt)
    {
        ca = spt->m_ca[m_lsol];
        c = spt->m_c[m_lsol];
        dkdJ = spt->m_dkdJ[m_lsol];
    }
    else if (fspt)
    {
        ca = fspt->m_ca[m_lsol];
        c = fspt->m_c[m_lsol];
        dkdJ = fspt->m_dkdJ[m_lsol];
    }
    else if (mfpt)
    {
        ca = mfpt->m_ca[m_lsol];
        c = mfpt->m_c[m_lsol];
        dkdJ = mfpt->m_dkdJ[m_lsol];
    }
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
    
    // derivative of (J d) w.r.t. J
    double dJ = d*(1+J*(m_alphad*phi0/(J-phi0)/(J-phi0) - m_cdinv*c*dkdJ));
		
	tens4dmm D4 = dyad1s(I)*dJ-dyad4s(I)*(2*d);
	
	return D4;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffAlbroIso::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
    FEElasticMaterialPoint* et = mp.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint* ppt = mp.ExtractData<FEBiphasicMaterialPoint>();
    FEBiphasicFSIMaterialPoint* bfpt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    FESolutesMaterialPoint* spt = mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint* fspt = mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FEMultiphasicFSIMaterialPoint* mfpt = mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
	
	// relative volume
	double J = et->m_J;
	
	// solid volume fraction in reference configuration
    double phi0 = 0;
    if (ppt)
        phi0 = ppt->m_phi0t;
    else if (bfpt)
        phi0 = bfpt->m_phi0;
    // porosity in current configuration
    double phiw = 1 - phi0/J;
    // solute concentration
    double ca, c, dkdc, k;
    ca = c = dkdc = k = 0;
    if (spt)
    {
        ca = spt->m_ca[m_lsol];
        c = spt->m_c[m_lsol];
        dkdc = spt->m_dkdc[m_lsol][isol];
        k = spt->m_k[m_lsol];
    }
    else if (fspt)
    {
        ca = fspt->m_ca[m_lsol];
        c = fspt->m_c[m_lsol];
        dkdc = fspt->m_dkdc[m_lsol][isol];
        k = fspt->m_k[m_lsol];
    }
    else if (mfpt)
    {
        ca = mfpt->m_ca[m_lsol];
        c = mfpt->m_c[m_lsol];
        dkdc = mfpt->m_dkdc[m_lsol][isol];
        k = mfpt->m_k[m_lsol];
    }
    
    // diffusivity coefficient
    double d = m_diff0*exp(-m_alphad*(1-phiw)/phiw - m_cdinv*ca);
    // derivative of d w.r.t. actual concentration
    double dc = -m_cdinv*d;
    
    // tangent w.r.t. concentration
    if (isol == m_lsol) {
        return mat3dd(dc*(k+dkdc*c));
    } else
        return mat3dd(dc*dkdc*c);
}
