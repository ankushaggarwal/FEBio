#pragma once

#include "mat3d.h"
#include "DumpFile.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Material point class

//! This class implements the concept of a material point. This point carries
//! with it not only information about its location, both in the reference and  
//! current configuration but also about the local deformation. In addition
//! it contains the state information that is associated with the current
//! point.

class FEMaterialPoint
{
public:
	FEMaterialPoint(FEMaterialPoint* ppt = 0) : m_pt(ppt){}
	virtual ~FEMaterialPoint() { if (m_pt) { delete m_pt; m_pt = 0; } }

	//! The init function is used to intialize data
	virtual void Init(bool bflag) = 0;

	virtual FEMaterialPoint* Copy() = 0;

	virtual void Serialize(DumpFile& ar) = 0;

	template <class T>
	T* ExtractData()
	{
		T* p = dynamic_cast<T*>(this);
		if (p) return p; 
		else
		{
			if (m_pt) return m_pt->ExtractData<T>();
			else return 0;
		}
	}

protected:
	FEMaterialPoint*	m_pt;	//<! nested point data

public:
	static double time;	// time value
	static double dt; // time increment
};

//-----------------------------------------------------------------------------

class FEElasticMaterialPoint : public FEMaterialPoint
{
public:
	FEElasticMaterialPoint()
	{
		F.zero();
		Q.unit();
		J = 1;
		s.zero();
		s0.zero();
	}

	FEMaterialPoint* Copy()
	{
		FEElasticMaterialPoint* pt = new FEElasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << F << J << Q << s << s0;
		}
		else
		{
			ar >> F >> J >> Q >> s >> s0;
		}

		if (m_pt) m_pt->Serialize(ar);
	}

	mat3ds Strain();

	mat3ds RightCauchyGreen();
	mat3ds LeftCauchyGreen ();

	mat3ds DevRightCauchyGreen();
	mat3ds DevLeftCauchyGreen ();

	mat3ds pull_back(const mat3ds& A);
	mat3ds push_forward(const mat3ds& A);

public:
	void Init(bool bflag)
	{
		if (bflag)
		{
			F.unit();

			J = 1;

			s.zero();
			s0.zero();

//			Q.unit();
		}

		if (m_pt) m_pt->Init(bflag);
	}

public:
	// position 
	vec3d	r0;	//!< material position
	vec3d	rt;	//!< spatial position

	// deformation data
	mat3d	F;	//!< deformation gradient
	double	J;			//!< determinant8 of F
	mat3d	Q;			//!< local material orientation

	// solid material data
	mat3ds		s;			//!< Cauchy stress
	mat3ds		s0;			//!< Initial stress (only used by linear solid solver)
};
