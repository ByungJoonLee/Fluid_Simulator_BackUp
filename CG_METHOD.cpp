#include "stdafx.h"
#include "CG_METHOD.h"

void CG_METHOD::CGMethod(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	LOG::Begin(thread_id, "CGMethod");

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

	T *rval(res.values), *pval(p.values), *Apval(Ap.values), *xval(x.values);

	T alpha, res_old, res_new;

	A.ComputeResidual(thread_id, x, b, res);

	for (int i = start_ix; i <= end_ix; i++)
	{
		p.values[i] = res.values[i];
	}
	multithreading->Sync(thread_id);

	res_old = DotProduct(multithreading, thread_id, res, p);
	
	while(num_iteration < max_iteration)
	{
		A.Multiply(thread_id, p, Ap);

		alpha = res_old/ DotProduct(multithreading, thread_id, p, Ap);
		
		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] += alpha*pval[i];
			rval[i] -= alpha*Apval[i];
		}
		multithreading->Sync(thread_id);

		res_new = DotProduct(multithreading, thread_id, res, res);

		if(res_new < sqr_tolerance) break;					// In L2 Norm

		for (int i = start_ix; i <= end_ix; i++)
		{
			const T k = res_new/res_old;

			pval[i] = res.values[i] + k*p.values[i];
		}
		multithreading->Sync(thread_id);

		res_old = res_new;

        BEGIN_HEAD_THREAD_WORK
		{
			num_iteration++;
		}
		END_HEAD_THREAD_WORK
	}

	multithreading->Sync(thread_id);

	BEGIN_HEAD_THREAD_WORK
	{
		residual = sqrt(res_new);

		if(use_detailed_log)
		{
			LOG::cout << "[CG] Iteration = " << num_iteration << ", Residual = " << residual << endl;
		}
		else
		{
			LOG::cout << "Iteration = " << num_iteration << ", Residual = " << residual << endl;
		}
	}
	END_HEAD_THREAD_WORK

	LOG::End(thread_id);
}

void CG_METHOD::biCGSTABMethod(const int& thread_id, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	LOG::Begin(thread_id, "bi-CGMethod");

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
		res_0.Initialize(N);
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
		As.Initialize(N);
	}
	END_HEAD_THREAD_WORK
	
	BEGIN_HEAD_THREAD_WORK
	{
		num_iteration = 0;
	}
	END_HEAD_THREAD_WORK

	T *rval(res.values), *pval(p.values), *Apval(Ap.values), *xval(x.values), *sval(s.values), *Asval(As.values);

	T alpha, res_old, res_new, omega, beta, s_norm, residual_check;

	// Subvariables
	T deno_omega, nu_omega;

	A.ComputeResidual(thread_id, x, b, res);

	// The value of r_0^*
	for (int i = start_ix; i <= end_ix; i++)
	{
		res_0.values[i] = res.values[i];
	}
	multithreading->Sync(thread_id);

	for (int i = start_ix; i <= end_ix; i++)
	{
		p.values[i] = res.values[i];
	}
	multithreading->Sync(thread_id);

	res_old = DotProduct(multithreading, thread_id, res, res_0);

	while(num_iteration < max_iteration)
	{
		A.Multiply(thread_id, p, Ap);
	
		alpha = res_old/DotProduct(multithreading, thread_id, Ap, res_0);
		
		for (int i = start_ix; i <= end_ix; i++)
		{
			sval[i] = rval[i] - alpha*Apval[i];
		}
		multithreading->Sync(thread_id);

		s_norm = sqrt(DotProduct(multithreading, thread_id, s, s));

		if (s_norm < tolerance)
		{
			for (int i = start_ix; i <= end_ix; i++)
			{
				xval[i] = xval[i] + alpha*pval[i];
			}
			multithreading->Sync(thread_id);
			
			break;
		}

		A.Multiply(thread_id, s, As);

		DotProduct(multithreading, thread_id, As, s, nu_omega);
		DotProduct(multithreading, thread_id, As, As, deno_omega);
		
		omega = nu_omega/deno_omega;

		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] += alpha*pval[i] + omega*sval[i];
			rval[i] = sval[i] - omega*Asval[i];
		}
		multithreading->Sync(thread_id);

		res_new = DotProduct(multithreading, thread_id, res, res_0);

		residual_check = DotProduct(multithreading, thread_id, res, res);

		if(residual_check < sqr_tolerance) break;					// In L2 Norm

		beta = (res_new/res_old)*(alpha/omega);
		
		for (int i = start_ix; i <= end_ix; i++)
		{
			pval[i] = rval[i] + beta*(pval[i] - omega*Apval[i]);
		}
		multithreading->Sync(thread_id);

		res_old = res_new;

		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration++;
		}
		END_HEAD_THREAD_WORK
	}

	multithreading->Sync(thread_id);

	BEGIN_HEAD_THREAD_WORK
	{
		residual = sqrt(residual_check);

		if(use_detailed_log)
		{
			LOG::cout << "[bi-CG] Iteration = " << num_iteration << ", Residual = " << residual << endl;
		}
		else
		{
			LOG::cout << "Iteration = " << num_iteration << ", Residual = " << residual << endl;
		}
	}
	END_HEAD_THREAD_WORK

	LOG::End(thread_id);
}