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
#include "FEMassActionReversible.h"

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversible::FwdReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *pt.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // get forward reaction rate
    double k = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    int nsol = 0;
    
    // start with contribution from solutes
    if (m_pMP)
    {
        nsol = (int)spt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vR = m_vR[i];
            if (vR > 0) {
                double c = spt.m_ca[i];
                zhat *= pow(c, vR);
            }
        }
    }
    else if (m_pFS)
    {
        nsol = (int)fspt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vR = m_vR[i];
            if (vR > 0) {
                double c = fspt.m_ca[i];
                zhat *= pow(c, vR);
            }
        }
    }
    else if (m_pSM)
    {
        nsol = (int)smpt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vR = m_vR[i];
            if (vR > 0) {
                double c = smpt.m_ca[i];
                zhat *= pow(c, vR);
            }
        }
    }
    else if (m_pMF)
    {
        nsol = (int)mfpt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vR = m_vR[i];
            if (vR > 0) {
                double c = mfpt.m_ca[i];
                zhat *= pow(c, vR);
            }
        }
    }
    
    if (m_pMP)
    {
        // add contribution of solid-bound molecules
        const int nsbm = (int)spt.m_sbmr.size();
        for (int i=0; i<nsbm; ++i) {
            int vR = m_vR[nsol+i];
            if (vR > 0) {
                double c = m_pMP->SBMConcentration(pt, i);
                zhat *= pow(c, vR);
            }
        }
    }
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversible::RevReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *pt.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // get forward reaction rate
    double k = m_pRev->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    // start with contribution from solutes
    int nsol = 0;
    if (m_pMP)
    {
        nsol = (int)spt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vP = m_vP[i];
            if (vP > 0) {
                double c = spt.m_ca[i];
                zhat *= pow(c, vP);
            }
        }
    }
    else if (m_pFS)
    {
        nsol = (int)fspt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vP = m_vP[i];
            if (vP > 0) {
                double c = fspt.m_ca[i];
                zhat *= pow(c, vP);
            }
        }
    }
    else if (m_pSM)
    {
        nsol = (int)smpt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vP = m_vP[i];
            if (vP > 0) {
                double c = smpt.m_ca[i];
                zhat *= pow(c, vP);
            }
        }
    }
    else if (m_pMF)
    {
        nsol = (int)mfpt.m_ca.size();
        for (int i=0; i<nsol; ++i) {
            int vP = m_vP[i];
            if (vP > 0) {
                double c = mfpt.m_ca[i];
                zhat *= pow(c, vP);
            }
        }
    }
    
    if (m_pMP)
    {
        // add contribution of solid-bound molecules
        const int nsbm = (int)spt.m_sbmr.size();
        for (int i=0; i<nsbm; ++i) {
            int vP = m_vP[nsol+i];
            if (vP > 0) {
                double c = m_pMP->SBMConcentration(pt, i);
                zhat *= pow(c, vP);
            }
        }
    }
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversible::ReactionSupply(FEMaterialPoint& pt)
{
	double zhatF = FwdReactionSupply(pt);
	double zhatR = RevReactionSupply(pt);
	return zhatF - zhatR;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionReversible::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    FEBiphasicFSIMaterialPoint& bfpt = *pt.ExtractData<FEBiphasicFSIMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *pt.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
	
	const int nsol = m_nsol;
	const int nsbm = (int)m_v.size() - nsol;
    double J = 0.0;
    double phi0 = 0.0;
    if (m_pMF)
    {
        J = ept.m_J;
        phi0 = bfpt.m_phi0;
    }
    else
    {
        J = ept.m_J;
        phi0 = bpt.m_phi0t;
    }
	
	
	// forward reaction
	double kF = m_pFwd->ReactionRate(pt);
	mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
	double zhatF = FwdReactionSupply(pt);
	mat3ds dzhatFde = mat3dd(0);
	if (kF > 0) {
		dzhatFde += dkFde/kF;
	}
	mat3ds I = mat3dd(1);
	for (int isol=0; isol<nsol; ++isol)
    {
        if (m_pMP)
            dzhatFde += I*(m_vR[isol]*spt.m_dkdJ[isol]/spt.m_k[isol]);
        else if (m_pFS)
            dzhatFde += I*(m_vR[isol]*fspt.m_dkdJ[isol]/fspt.m_k[isol]);
        else if (m_pSM)
            dzhatFde += I*(m_vR[isol]*smpt.m_dkdJ[isol]/smpt.m_k[isol]);
        else if (m_pMF)
            dzhatFde += I*(m_vR[isol]*mfpt.m_dkdJ[isol]/mfpt.m_k[isol]);
    }
	for (int isbm = 0; isbm<nsbm; ++isbm)
		dzhatFde += I*(m_vR[nsol+isbm]/(J-phi0));
	
	dzhatFde *= zhatF;
	
	// reverse reaction
	double kR = m_pRev->ReactionRate(pt);
	mat3ds dkRde = m_pRev->Tangent_ReactionRate_Strain(pt);
	double zhatR = RevReactionSupply(pt);
	mat3ds dzhatRde = mat3dd(0);
	if (kR > 0) {
		dzhatRde += dkRde/kR;
	}
	for (int isol=0; isol<nsol; ++isol)
    {
        if (m_pMP)
            dzhatRde += I*(m_vP[isol]*spt.m_dkdJ[isol]/spt.m_k[isol]);
        else if (m_pFS)
            dzhatRde += I*(m_vP[isol]*fspt.m_dkdJ[isol]/fspt.m_k[isol]);
        else if (m_pSM)
            dzhatRde += I*(m_vP[isol]*smpt.m_dkdJ[isol]/smpt.m_k[isol]);
        else if (m_pMF)
            dzhatRde += I*(m_vP[isol]*mfpt.m_dkdJ[isol]/mfpt.m_k[isol]);
    }
	for (int isbm = 0; isbm<nsbm; ++isbm)
		dzhatRde -= I*(m_vP[nsol+isbm]/(J-phi0));
	
	dzhatRde *= zhatR;
	
	return dzhatFde - dzhatRde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionReversible::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
	// forward reaction
	double kF = m_pFwd->ReactionRate(pt);
	double dzhatFdp = 0;
	if (kF > 0) {
		double dkFdp = m_pFwd->Tangent_ReactionRate_Pressure(pt);
		double zhatF = FwdReactionSupply(pt);
		dzhatFdp = dkFdp*zhatF/kF;
	}
	
	// reverse reaction
	double kR = m_pRev->ReactionRate(pt);
	double dzhatRdp = 0;
	if (kR > 0) {
		double dkRdp = m_pRev->Tangent_ReactionRate_Pressure(pt);
		double zhatR = RevReactionSupply(pt);
		dzhatRdp = dkRdp*zhatR/kR;
	}
	
	return dzhatFdp - dzhatRdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionReversible::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    const int nsol = m_nsol;
    
    // if the derivative is taken with respect to a solid-bound molecule, return 0
    if (sol >= nsol) {
        return 0;
    }
    
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *pt.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // forward reaction
    double zhatF = FwdReactionSupply(pt);
    double dzhatFdc = 0;
    for (int isol=0; isol<nsol; ++isol) {
        if (m_pMP)
        {
            dzhatFdc += m_vR[isol]*spt.m_dkdc[isol][sol]/spt.m_k[isol];
            if ((isol == sol) && (spt.m_c[sol] > 0))
                dzhatFdc += m_vR[isol]/spt.m_c[sol];
        }
        else if (m_pFS)
        {
            dzhatFdc += m_vR[isol]*fspt.m_dkdc[isol][sol]/fspt.m_k[isol];
            if ((isol == sol) && (fspt.m_c[sol] > 0))
                dzhatFdc += m_vR[isol]/fspt.m_c[sol];
        }
        else if (m_pSM)
        {
            dzhatFdc += m_vR[isol]*smpt.m_dkdc[isol][sol]/smpt.m_k[isol];
            if ((isol == sol) && (smpt.m_c[sol] > 0))
                dzhatFdc += m_vR[isol]/smpt.m_c[sol];
        }
        else if (m_pMF)
        {
            dzhatFdc += m_vR[isol]*mfpt.m_dkdc[isol][sol]/mfpt.m_k[isol];
            if ((isol == sol) && (mfpt.m_c[sol] > 0))
                dzhatFdc += m_vR[isol]/mfpt.m_c[sol];
        }
    }
    
    dzhatFdc *= zhatF;
    
    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double dzhatRdc = 0;
    for (int isol=0; isol<nsol; ++isol) {
        if (m_pMP)
        {
            dzhatRdc += m_vP[isol]*spt.m_dkdc[isol][sol]/spt.m_k[isol];
            if ((isol == sol) && (spt.m_c[sol] > 0))
                dzhatRdc += m_vP[isol]/spt.m_c[sol];
        }
        else if (m_pFS)
        {
            dzhatRdc += m_vP[isol]*fspt.m_dkdc[isol][sol]/fspt.m_k[isol];
            if ((isol == sol) && (fspt.m_c[sol] > 0))
                dzhatRdc += m_vP[isol]/fspt.m_c[sol];
        }
        else if (m_pSM)
        {
            dzhatRdc += m_vP[isol]*smpt.m_dkdc[isol][sol]/smpt.m_k[isol];
            if ((isol == sol) && (smpt.m_c[sol] > 0))
                dzhatRdc += m_vP[isol]/smpt.m_c[sol];
        }
        else if (m_pMF)
        {
            dzhatRdc += m_vP[isol]*mfpt.m_dkdc[isol][sol]/mfpt.m_k[isol];
            if ((isol == sol) && (mfpt.m_c[sol] > 0))
                dzhatRdc += m_vP[isol]/mfpt.m_c[sol];
        }
    }
    
    dzhatRdc *= zhatR;
    
    return dzhatFdc - dzhatRdc;
}
