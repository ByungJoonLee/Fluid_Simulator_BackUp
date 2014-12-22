#pragma once

#include "LINEAR_SOLVER.h"
#include "LOG.h"

class GAUSS_SEIDEL_METHOD : public LINEAR_SOLVER
{
public:
	typedef LINEAR_SOLVER BASE;

public:
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
	using BASE::multithreading;
	
public: // Constructor and Destructor
	GAUSS_SEIDEL_METHOD(void)
	{}

	~GAUSS_SEIDEL_METHOD(void)
	{}

public: // Solver
	void Solve(const int& thread_id, const CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, const VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<int>& bc)
	{
		GaussSeidelMethod(thread_id, A_matrix, x_vector, b_vector);
	}

	void GaussSeidelMethod(const int& thread_id, const CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, const VECTOR_ND<T>& b_vector)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			multithreading->SplitDomainIndex1D(0, A_matrix.N);
		}
		END_HEAD_THREAD_WORK
		
		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration = 0;
		}
		END_HEAD_THREAD_WORK

		residual = (T)0;

		for (int i = 0; i < max_iteration; i++)
		{
			residual = GaussSeidelStep(thread_id, A_matrix, x_vector, b_vector);

			BEGIN_HEAD_THREAD_WORK
			{
				num_iteration++;
			}
			END_HEAD_THREAD_WORK
			
			if(residual < tolerance) break;
		}
		
		multithreading->Sync(thread_id);

		BEGIN_HEAD_THREAD_WORK
		{
			LOG::cout << "GS method iteration = " << num_iteration << " residual = " << residual << endl;
		}
		END_HEAD_THREAD_WORK
	}

	T GaussSeidelStep(const int& thread_id, const CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, const VECTOR_ND<T>& div)
	{
		assert(A_matrix.N == x_vector.num_dimension);
		assert(x_vector.num_dimension == div.num_dimension);

		int k_start(multithreading->start_ix_1D[thread_id]), k_end(multithreading->end_ix_1D[thread_id]);

		T *divval(div.values), *xval(x_vector.values);

		int v_start, v_end, vix;
		T v, A_ii, residual, residual_sum(0);

		for (int row = k_start; row <= k_end; row++)
		{
			v_start = A_matrix.row_ptr[row];
			v_end = A_matrix.row_ptr[row + 1] - 1;

			v = (T)0;
			A_ii = (A_matrix.values[A_matrix.row_ptr[row + 1] - 1]);		// Note: the last value of each row is the diagonal term in our POISSON_SOLVER_3D::BuildLinearSystem function
			
			for (vix = v_start; vix <= v_end; vix++)
			{
				v += A_matrix.values[vix]*xval[A_matrix.column_index[vix]];
			}

			residual = divval[row] - v;

			// Gauss-Seidel Procedure : x = (D-L)^(-1)(Ux + b) = x + (D-L)^(-1)(b-Ax) = x + (D-L)^(-1)(residual)
			if (A_ii != 0)
			{
				xval[row] += residual/A_ii;
			}

			residual_sum += POW2(residual);
		}

		multithreading->SyncSum(thread_id, residual_sum);

		residual_sum = sqrt(residual_sum);

		return residual_sum;
	}
};