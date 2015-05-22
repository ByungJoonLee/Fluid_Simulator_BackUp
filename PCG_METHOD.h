#pragma once

#include "LINEAR_SOLVER.h"
#include "CG_METHOD.h"

class PCG_METHOD : public LINEAR_SOLVER
{
public: // Typedef
	typedef LINEAR_SOLVER BASE;

public: // Using Keyword
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
    using BASE::one_point_fixed;
    using BASE::projecting_space_without_null;

public: // Essential Data
	VECTOR_ND<T> res, p, Ap, s;

public: // Subdata for Multithreading
	CSR_MATRIX<T> M;
	VECTOR_ND<T>  y;

public: // Constructor and Destructor
	PCG_METHOD(void)
	{}

	~PCG_METHOD(void)
	{}

public: // Initialization Function
	void Initialize(MULTITHREADING* multithreading_input, const T& tolerance_input, const int& max_iteration_input)
	{
		BASE::Initialize(multithreading_input, tolerance_input, max_iteration_input);
	}

public: // Solver
	void Solve(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_3D<int>& bc)
	{
		if (projecting_space_without_null)
		{
            PCGMethodProjected(thread_id, bc.grid.i_res, bc.grid.j_res, bc.grid.k_res, A, x, b);
		}
		else
		{
            PCGMethod(thread_id, bc.grid.i_res, bc.grid.j_res, bc.grid.k_res, A, x, b);
		}
	}

    void PCGMethod(const int& thread_id, const int& i_res_input, const int& j_res_input, const int& k_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void PCGMethodProjected(const int & thread_id, const int& i_res_input, const int& j_res_input, const int& k_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void MultiplicationByMinverseAsDiagonal(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void MultiplicationByMinverseAsDiagonalProjected(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const T& constant_value);
	void MultiplicationByMinverse(const int& thread_id, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
};