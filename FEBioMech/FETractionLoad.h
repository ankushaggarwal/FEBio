#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! FETractionLoad is a surface that has a constant (deformation independant)
//! traction force on it.
//!
class FETractionLoad : public FESurfaceLoad
{
public:
	struct LOAD
	{
		vec3d	s[8];		// nodal scale factors
		int		lc;			// load curve
	};

public:
	//! constructor
	FETractionLoad(FEModel* pfem) : FESurfaceLoad(pfem) {}

	//! allocate storage
	void Create(int n) { m_TC.resize(n); }

	//! get a traction load BC
	LOAD& TractionLoad(int n) { return m_TC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver) {}

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

public:
	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

protected:
	vector<LOAD>	m_TC;		//!< traction boundary cards
};
