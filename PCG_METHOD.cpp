#include "stdafx.h"
#include "PCG_METHOD.h"

void PCG_METHOD::PCGMethod(const int& thread_id, const int& i_res_input, const int& j_res_input, const int& k_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	LOG::Begin(thread_id, "PCGMethod");

	BEGIN_HEAD_THREAD_WORK
	{
		multithreading->SplitDomainIndex1D(0, A.N);
	}
	END_HEAD_THREAD_WORK
	
	const int N(x.num_dimension), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

	BEGIN_HEAD_THREAD_WORK
	{
		res.Initialize(N);
	}
	END_HEAD_THREAD_WORK

	BEGIN_HEAD_THREAD_WORK
	{
		s.Initialize(N);
	}
	END_HEAD_THREAD_WORK

	BEGIN_HEAD_THREAD_WORK
	{
		p.Initialize(N);
	}
	END_HEAD_THREAD_WORK
	
	BEGIN_HEAD_THREAD_WORK
	{
		Ap.Initialize(N);
	}
	END_HEAD_THREAD_WORK

	BEGIN_HEAD_THREAD_WORK
	{
		num_iteration = 0;
	}
	END_HEAD_THREAD_WORK

	T *rval(res.values), *pval(p.values), *Apval(Ap.values), *xval(x.values), *sval(s.values);

	T alpha, beta, res_old, res_new, dot_result;

	A.ComputeResidual(thread_id, x, b, res);

	//IncompleteCholeskyDecomposition(multithreading, i_res_input, j_res_input, k_res_input, A, M);

	//MultiplicationByMinverse(thread_id, M, p, res);

	MultiplicationByMinverseAsDiagonal(thread_id, A, p, res);

	DotProduct(multithreading, thread_id, res, p, res_new);
		
	while (num_iteration < max_iteration)
	{
		A.Multiply(thread_id, p, Ap);

		DotProduct(multithreading, thread_id, p, Ap, dot_result);

		alpha = res_new/dot_result;

		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] += alpha*pval[i];
		}
		multithreading->Sync(thread_id);
		
		for (int i = start_ix; i <= end_ix; i++)
		{
			rval[i] -= alpha*Apval[i];
		}
		multithreading->Sync(thread_id);

		//MultiplicationByMinverse(thread_id, M, p, res);
		MultiplicationByMinverseAsDiagonal(thread_id, A, s, res);

		res_old = res_new;

		DotProduct(multithreading, thread_id, res, s, res_new);

		if (res_new < sqr_tolerance)
		{
			break;
		}

		beta = res_new/res_old;

		for (int i = start_ix; i <= end_ix; i++)
		{
			pval[i] *= beta;
			pval[i] += sval[i];
		}
		multithreading->Sync(thread_id);

		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration++;
		}
		END_HEAD_THREAD_WORK;
	}
	multithreading->Sync(thread_id);

	BEGIN_HEAD_THREAD_WORK
	{
		residual = sqrt(res_new);
		

		if (use_detailed_log)
		{
			LOG::cout << "[PCG] Iteration = " << num_iteration << ", Residual = " << residual << endl;
		}
		else
		{
			LOG::cout << "Iteration = " << num_iteration << ", Residual = " << residual << endl;
		}
	}
	END_HEAD_THREAD_WORK

	LOG::End(thread_id);
}

void PCG_METHOD::MultiplicationByMinverseAsDiagonal(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	const int N = A.N;

	BEGIN_HEAD_THREAD_WORK
	{
		M.Initialize(multithreading, N, N);
		multithreading->SplitDomainIndex1D(0, N);
		
		M.start_ix[0] = 0;
		M.end_ix[0] = multithreading->sync_value_int[0] - 1;
		M.prev_row_array[0] = -1;
		M.values_ix_array[0] = 0;
		for (int id = 1; id < multithreading->num_threads; id++)
		{
			M.start_ix[id] = M.end_ix[id - 1] + 1;
			M.end_ix[id] = M.end_ix[id - 1] + multithreading->sync_value_int[id];
			M.prev_row_array[id] = -1;
			M.values_ix_array[id] = M.start_ix[id];
		}
	}
	END_HEAD_THREAD_WORK

	DiagonalPreconditioner(multithreading, thread_id, A, M);
		
	const int start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

	for (int i = start_ix; i <= end_ix; i++)
	{
		x[i] = b[i]/M(i, i);
	}
	multithreading->Sync(thread_id);
}

void PCG_METHOD::MultiplicationByMinverse(const int& thread_id, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	BEGIN_HEAD_THREAD_WORK
	{
		y.Initialize(x.num_dimension, true);
	}
	END_HEAD_THREAD_WORK
	 
	T one_over_M_start((T)1/M(0,0));
	
	int number(0), num_2(0), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);;
	
	for (int i = start_ix; i <= end_ix; i++)
	{
		T sum(0);

		for (int k = M.row_ptr[i]; k < (M.row_ptr[i + 1] - 1); k++)
		{
			sum += M.values[k]*y[M.column_index[k]];
		}
			
		T one_over_Mii = 1/M(i, i);
		y[i] = (b[i] - sum)*one_over_Mii;
	}
	multithreading->Sync(thread_id);

	T* summation = new T[x.num_dimension];

	for (int i = start_ix; i <= end_ix; i++)
	{
		summation[i] = 0;
	}
	multithreading->Sync(thread_id);

	// Matrix-Transpose-Vector Multiplication
	T one_over_M_end = 1/M(x.num_dimension - 1, x.num_dimension - 1);
	
	BEGIN_HEAD_THREAD_WORK
	{
		x[x.num_dimension - 1] = y[x.num_dimension - 1]*one_over_M_end;
	}
	END_HEAD_THREAD_WORK

	for (int i = end_ix; i >= start_ix; i--)
	{
		for (int k = M.row_ptr[i + 1] - 1; k >= M.row_ptr[i]; k--)
		{
			if (k != M.row_ptr[i + 1] - 1)
			{
				summation[M.column_index[k]] += M.values[k]*x[i];
			}
		}
			
		T one_over_M = 1/M(i - 1, i - 1);
		x[i - 1] = (y[i - 1] - summation[i - 1])*one_over_M;
	}
	multithreading->Sync(thread_id);

	delete summation;
}

