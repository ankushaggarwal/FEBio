/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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


#pragma once
#include "FEElasticMaterial.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include "FEReactivePlasticDamageMaterialPoint.h"

//-----------------------------------------------------------------------------
// This material models damage in any reactive plastic materials.

class FEReactivePlasticDamage : public FEElasticMaterial
{
public:
	FEReactivePlasticDamage(FEModel* pfem);
    
public:
    //! data initialization and checking
    bool Init() override;
    
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
    
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
    
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! evaluate elastic deformation gradient
    void ElasticDeformationGradient(FEMaterialPoint& pt);
    
    //! damage
    double Damage(FEMaterialPoint& pt, int k);
    
    //!< update fatigue material point at each iteration
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& pt, const FETimeInfo& tp) override;
    
	// returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override;
    
    // get the elastic material
    FEElasticMaterial* GetElasticMaterial() override { return m_pBase; }
    
    // get the yield surface normal
    mat3ds YieldSurfaceNormal(FEElasticMaterialPoint& pe);
    
    // evaluate octahedral plastic strain
    void OctahedralPlasticStrain(FEMaterialPoint& pt);
    
    // evaluate reactive heat supply
    void ReactiveHeatSupplyDensity(FEMaterialPoint& pt);
    void Serialize(DumpStream& ar);
public:
    FEElasticMaterial*  m_pBase;     // base elastic material
    FEDamageCriterion*  m_pCrit;     // yield criterion
    FEDamageCDF*        m_pYDamg;    // yield damage model
    FEDamageCriterion*  m_pYDCrit;   // yield damage criterion
    FEDamageCDF*        m_pIDamg;    // intact damage model
    FEDamageCriterion*  m_pIDCrit;   // intact damage criterion
    
private:
    vector<double>      Ky;
    vector<double>      w;
    
public:
    double      m_wmin;     // initial fraction of yielding bonds
    double      m_wmax;     // final fraction of yielding bonds
    double      m_we;       // fraction of unyielding bonds
    double      m_Ymin;     // initial yield measure
    double      m_Ymax;     // yield measure when all bonds have yielded
    int         m_n;        // number of yield levels
    bool        m_isochrc;  // flag for constraining plastic def grad to be isochoric
    double      m_rtol;     // user-defined relative tolerance
    double      m_bias;     // biasing factor for intervals in yield measures and bond fractions
    
    DECLARE_FECORE_CLASS();
};
