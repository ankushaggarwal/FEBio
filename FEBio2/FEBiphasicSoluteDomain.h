#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic-solute 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSoluteDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEBiphasicSoluteDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_BIPHASIC_SOLUTE_DOMAIN; }
	
	FEDomain* Clone()
	{
		FEBiphasicSoluteDomain* pd = new FEBiphasicSoluteDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);
	
	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);
	
	// update stresses
	void UpdateStresses(FEModel& fem);
	
protected:
	//! Calculates the internal fluid forces
	bool InternalFluidWork(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! Calculates the internal solute forces
	bool InternalSoluteWork(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! calculates the element solute-poroelastic stiffness matrix
	bool ElementBiphasicSoluteStiffness(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! calculates the solid element stiffness matrix
	void SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void BiphasicSoluteMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);
};
