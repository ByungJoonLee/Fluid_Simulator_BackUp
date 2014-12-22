#pragma once

#include "LINEAR_SOLVER.h"

class CG_METHOD : public LINEAR_SOLVER
{
public: // Typedef
	typedef LINEAR_SOLVER BASE;

public: // Using keyword
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
	using BASE::multithreading;
	
public: // Essential Data
	VECTOR_ND<T> res, p, Ap;		
	
public: // For Bi-CG
	VECTOR_ND<T> res_0, s, As;
	
public: // Constructor and Destructor
	CG_METHOD(void)
	{}

	~CG_METHOD(void)
	{}

public: // Solver
	void Solve(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_3D<int>& bc)
	{
		CGMethod(thread_id, A, x, b);
	}

	void CGMethod(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void biCGSTABMethod(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
};




					 