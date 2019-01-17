#include "stdafx.h"
#include "ILU0_Preconditioner.h"
#include <FECore/CompactUnSymmMatrix.h>

// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

//=================================================================================================

ILU0_Preconditioner::ILU0_Preconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;

	m_K = 0;
}

bool ILU0_Preconditioner::Create()
{
	m_K = dynamic_cast<CRSSparseMatrix*>(GetSparseMatrix());
	if (m_K == 0) return false;
	assert(m_K->Offset() == 1);

	int N = m_K->Rows();
	int NNZ = m_K->NonZeroes();

	double* pa = m_K->Values();
	int* ia = m_K->Pointers();
	int* ja = m_K->Indices();

	MKL_INT ipar[128] = { 0 };
	double dpar[128] = { 0.0 };

	// parameters affecting the pre-conditioner
	if (m_checkZeroDiagonal)
	{
		ipar[30] = 1;
		dpar[30] = m_zeroThreshold;
		dpar[31] = m_zeroReplace;
	}

	m_tmp.resize(N, 0.0);

	m_bilu0.resize(NNZ);
	int ierr = 0;
	dcsrilu0(&N, pa, ia, ja, &m_bilu0[0], ipar, dpar, &ierr);
	if (ierr != 0) return false;

	return true;
}

bool ILU0_Preconditioner::mult_vector(double* x, double* y)
{
	int ivar = m_K->Rows();
	int* ia = m_K->Pointers();
	int* ja = m_K->Indices();

	char cvar1 = 'L';
	char cvar = 'N';
	char cvar2 = 'U';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilu0[0], ia, ja, &x[0], &m_tmp[0]);
	cvar1 = 'U';
	cvar = 'N';
	cvar2 = 'N';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilu0[0], ia, ja, &m_tmp[0], &y[0]);

	return true;
}
