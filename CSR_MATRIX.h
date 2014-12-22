#pragma once

#include "COMMON_DEFINITIONS.h"
#include "VECTOR_ND.h"
#include "MULTITHREADING.h"

// compressed sparse row (CSR or CRS) http://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR_or_CRS.29
// see http://developer.download.nvidia.com/compute/DevZone/docs/html/CUDALibraries/doc/CUSPARSE_Library.pdf for CUDA compatable descriptions.

template <class T>
class CSR_MATRIX
{
public:
	int						N;					// the number of row
	int						nz;					// the number of nonzero elements
	T*						values;				// the values of the nonzero elements
	int*					row_ptr;			// the locations in the val vector that start a row
	int*					column_index;		// the column indexes of the elements in the val vector
	
	int						values_ix;
	int						row_count;			// row_count is not used.
	int						prev_row;

	MULTITHREADING*			multithreading;

	int						*start_ix, *end_ix;
	int*					prev_row_array;
	int*					values_ix_array;

public:
	CSR_MATRIX(void)
		: values(0), column_index(0), row_ptr(0), multithreading(0), start_ix(0), end_ix(0), prev_row_array(0), values_ix_array(0)
	{}

	~CSR_MATRIX(void)
	{
		DeleteMemory();
	}

	void DeleteMemory()
	{
		DELETE_ARRAY(values);
		DELETE_ARRAY(row_ptr);
		DELETE_ARRAY(column_index);
		DELETE_ARRAY(start_ix);
		DELETE_ARRAY(end_ix);
		DELETE_ARRAY(prev_row_array);
		DELETE_ARRAY(values_ix_array);
	}

public: // Initialization Functions
	void Initialize(MULTITHREADING* multithreading_input, const int &N_input, const int &nz_input)		// N by N square matrix. nz: number of non zero elements
	{
		DeleteMemory();

		multithreading = multithreading_input;
		
		start_ix = new int [multithreading->num_threads];
		end_ix = new int [multithreading->num_threads];
		prev_row_array = new int [multithreading->num_threads];
		values_ix_array = new int [multithreading->num_threads];

		N = N_input;
		nz = nz_input;

		values = new T [nz];
		row_ptr = new int [N+1];				// +1 is for final row iteration
		column_index = new int [nz];

		values_ix = 0;
		row_count = 0;
		prev_row = -1;

		row_ptr[N_input] = nz_input;
	}

	void Initialize(MULTITHREADING* multithreading_input, const CSR_MATRIX<T>& matrix_input)
	{
		Initialize(multithreading, matrix_input.N, matrix_input.nz);

		for (int i = 0; i < nz; i++)
		{
			values[i] = matrix_input.values[i];
		}

		for (int i = 0; i < N + 1; i++)
		{
			row_ptr[i] = matrix_input.row_ptr[i];
		}

		for (int i = 0; i < nz; i++)
		{
			column_index[i] = matrix_input.column_index[i];
		}
		
		for (int i = 0; i < multithreading->num_threads; i++)
		{
			start_ix[i] = matrix_input.start_ix[i];
			end_ix[i] = matrix_input.end_ix[i];
			prev_row_array[i] = matrix_input.prev_row_array[i];
			values_ix_array[i] = matrix_input.values_ix_array[i];
		}

		value_ix = matrix_input.value_ix;
		prev_row = matrix_input.prev_row;
	}

public: // Operator Overloading
	void operator*=(const T& s)
	{
		for(int i = 0; i < nz; i++)
		{
			values[i] *= s;
		}
	}

	inline T& operator()(const int& row_input, const int& column_input) const
	{
		bool is_nonzero(false);
		
		int vix;
		vix = row_ptr[row_input];

		while (true)
		{
			if (column_index[vix] == column_input)
			{
				return values[vix];
				is_nonzero = true;
			}
			else
			{
				vix += 1;
			}
		}

		/*for (int vix = row_ptr[row_input]; vix < row_ptr[row_input + 1]; vix++)
		{
			if (column_index[vix] == column_input)
			{
				return values[vix];
				is_nonzero = true;
			}
		}*/

		if (is_nonzero == false)
		{
			T zero(0);
			T& temp = zero;
			return temp;
		}
	}
	
	CSR_MATRIX<T>& operator=(const CSR_MATRIX<T>& matrix_input)
	{
		Initialize(multithreading, matrix_input);

		return (*this);
	}

public: // Member Functions
	void AssignValue(const int& row_input, const int& column_input, const T& values_input)
	{
		values[values_ix] = values_input;

		if(row_input != prev_row)
		{
//			assert(row_input == prev_row_+1);

			row_ptr[row_input] = values_ix;
			prev_row = row_input;
		}

		column_index[values_ix] = column_input;

		values_ix ++;
	}

	void AssignValue(const int& thread_id, const int& row_input, const int& column_input, const T& values_input)
	{
		values[values_ix_array[thread_id]] = values_input;

		if(row_input != prev_row_array[thread_id])
		{
			//check whether the matrix is well-made or not.
			//it can be omitted in simulation stage.
//			assert(row_input == prev_row_array_[thread_id]+1);

			row_ptr[row_input] = values_ix_array[thread_id];
			prev_row_array[thread_id] = row_input;
		}

		column_index[values_ix_array[thread_id]] = column_input;

		values_ix_array[thread_id] ++;
	}
	
	inline void Multiply(const VECTOR_ND<T>& x, VECTOR_ND<T>& b) const // this_matrix * x -> b
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);

		T *bval(b.values), *xval(x.values);

		for(int row = 0; row < N; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			T v=0;
			for(int vix = row_ptr[row]; vix < row_ptr[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values[vix]*xval[column_index[vix]];
			}

			bval[row] = v;
		}
	}

	inline void Multiply(const VECTOR_ND<T>& x, VECTOR_ND<T>& b, const int& k_start, const int& k_end) const // this_matrix * x -> b
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);
		T *bval(b.values), *xval(x.values);

		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			T v=0;
			for(int vix = row_ptr[row]; vix < row_ptr[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values[vix]*xval[column_index_[vix]];
			}

			bval[row] = v;
		}
	}

	inline void Multiply(const int& thread_id, const VECTOR_ND<T>& x, VECTOR_ND<T>& b) const // this_matrix * x -> b
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);

		T *bval(b.values), *xval(x.values);

		const int k_start(multithreading->start_ix_1D[thread_id]), k_end(multithreading->end_ix_1D[thread_id]);
		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			T v=0;
			const int vix_start = row_ptr[row];	assert(vix_start >= 0);
			const int vix_end = row_ptr[row+1];	assert(vix_start < nz);
			for(int vix = vix_start; vix < vix_end; vix ++) // iterate all components of 'row'th row of this matrix
			{
				assert(column_index[vix] < N);

				v += values[vix]*xval[column_index[vix]];
			}

			bval[row] = v;
		}

		multithreading->Sync(thread_id);
	}

	inline void ComputeResidual(const VECTOR_ND<T>& x, const VECTOR_ND<T>& b, VECTOR_ND<T>& residual) const // residual = b - this_matrix*x
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);
		assert(residual.num_dimension == N);

		T *bval(b.values), *xval(x.values), *rval(residual.values);

		for(int row = 0; row < N; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			// compute A*x
			T v=0;
			for(int vix = row_ptr[row]; vix < row_ptr[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values[vix]*xval[column_index[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}
	}

	inline void ComputeResidual(const VECTOR_ND<T>& x, const VECTOR_ND<T>& b, VECTOR_ND<T>& residual, const int& k_start, const int& k_end) const // residual = b - this_matrix*x
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);
		assert(residual.num_dimension == N);

		// speed-up pointers
		T *bval(b.values), *xval(x.values), *rval(residual.values);

		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			// compute A*x
			// TODO: we may optimize row_ptr_[row] and row_ptr_[row+1] access
			T v=0;
			for(int vix = row_ptr[row]; vix < row_ptr[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values[vix]*xval[column_index[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}
	}

	inline void ComputeResidual(const int& thread_id, const VECTOR_ND<T>& x, const VECTOR_ND<T>& b, VECTOR_ND<T>& residual) const // residual = b - this_matrix*x
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);
		assert(residual.num_dimension == N);

		// speed-up pointers
		T *bval(b.values), *xval(x.values), *rval(residual.values);

		const int k_start(multithreading->start_ix_1D[thread_id]), k_end(multithreading->end_ix_1D[thread_id]);
		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			// compute A*x
			// TODO: we may optimize row_ptr_[row] and row_ptr_[row+1] access
			T v=0;
			for(int vix = row_ptr[row]; vix < row_ptr[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values[vix]*xval[column_index[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}

		multithreading->Sync(thread_id);
	}
	
	T* GetValue(const int& row, const int& column)
	{
		for(int vix = row_ptr[row]; vix < row_ptr[row+1]; vix ++)
		{
			if(column_index[vix] == column) 
			{
				return &values[vix];
			}
		}
	}
};

inline std::ostream& operator<<(std::ostream& output, const CSR_MATRIX<T>& A)
{
	output << "--Matrix Information in CSR Form--" << endl;
	
	output << "nz = [ ";
	for (int i = 0; i < A.nz; i++)
	{
		output << A.values[i] << " ";
	}
	output << "]" << endl;
	
	output << "ci = [ ";
	for (int i = 0; i < A.nz; i++)
	{
		output << A.column_index[i] << " ";
	}
	output << "]" << endl;

	output << "rp = [ ";
	for (int i = 0; i <= A.N; i++)
	{
		output << A.row_ptr[i] << " ";
	}
	output << "]" << endl;

	return output;
}

template<class TT>
static void IncompleteCholeskyDecomposition(MULTITHREADING* multithreading_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const CSR_MATRIX<TT>& A, CSR_MATRIX<TT>& L)
{
	const int N = A.N;
	const int nz = A.nz;

	FIELD_STRUCTURE_3D<int> index_field;
	index_field.Initialize(multithreading_input, i_res_input, j_res_input, k_res_input, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3);

	const int i_start(index_field.i_start), i_end(index_field.i_end), j_start(index_field.j_start), j_end(index_field.j_end), k_start(index_field.k_start), k_end(index_field.k_end);
	const int i_res(index_field.grid.i_res), j_res(index_field.grid.j_res), k_res(index_field.grid.k_res);

	int start_ix(0);

	GRID_ITERATION_3D(index_field.grid)
	{
		index_field(i, j, k) = -1;
	}

	GRID_ITERATION_3D(index_field.grid_ghost)
	{
		if (i < i_start || i > i_end || j < j_start || j > j_end || k < k_start || k > k_end)
		{
			index_field(i, j, k) = -1;
		}
	}

	GRID_ITERATION_3D(index_field.grid)
	{
		index_field(i, j, k) = start_ix++;
	}
	
	int number(0);

	int i, j, k;
	LOOPS_3D(i, j, k, i_start, j_start, k_start, i_end, j_end, k_end)
	{
		T sum(0), coef(0);

		if (index_field(i, j, k - 1) > -1)
		{
			if (index_field(i, j, k) == i_res*j_res)
			{
				coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res*j_res]]*(A(index_field(i, j, k - 1), index_field(i, j, k))); 
			}
			else if ((index_field(i, j, k) > i_res*j_res) && (index_field(i, j, k) <= 2*i_res*j_res))
			{
				if ((index_field(i, j, k) <= i_res*j_res + i_res) || (index_field(i, j, k) % i_res == 0))
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res*j_res] + 1]*(A(index_field(i, j, k - 1), index_field(i, j, k)));
				}
				else
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res*j_res] + 2]*(A(index_field(i, j, k - 1), index_field(i, j, k)));
				}
				
			}
			else 
			{
				if ((index_field(i, j, k) % i_res*j_res == 0))
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res*j_res] + 2]*(A(index_field(i, j, k - 1), index_field(i, j, k)));
				}
				else if ((index_field(i, j, k) % i_res*j_res > 0) && (index_field(i, j, k) % i_res*j_res < i_res))
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res*j_res] + 2]*(A(index_field(i, j, k - 1), index_field(i, j, k)));
				}
				else
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res*j_res] + 3]*(A(index_field(i, j, k - 1), index_field(i, j, k)));
				}
			}

			L.AssignValue(index_field(i, j, k), index_field(i, j, k - 1), coef);
			number += 1;
		}

		if (index_field(i, j - 1, k) > -1)
		{
			if (index_field(i, j, k) < i_res*j_res)
			{
				if (index_field(i, j, k) == i_res)
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res]]*(A(index_field(i, j - 1, k), index_field(i, j, k)));
				}
				else if (((index_field(i, j, k) > i_res) && (index_field(i, j, k) < 2*i_res)) || (index_field(i, j, k) % i_res == 0))
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res] + 1]*(A(index_field(i, j - 1, k), index_field(i, j, k)));
				}
				else 
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res] + 2]*(A(index_field(i, j - 1, k), index_field(i, j, k)));
				}
			}
			else
			{
				if (index_field(i, j, k) % i_res*j_res == i_res)
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res] + 1]*(A(index_field(i, j - 1, k), index_field(i, j, k)));
				}
				else if ((index_field(i, j, k) % i_res*j_res > i_res) && (index_field(i, j, k) % i_res*j_res <= i_res + j_res))
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res] + 2]*(A(index_field(i, j - 1, k), index_field(i, j, k)));
				}
				else
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - i_res] + 3]*(A(index_field(i, j - 1, k), index_field(i, j, k)));
				}
			}

			L.AssignValue(index_field(i, j, k), index_field(i, j - 1, k), coef);
			number += 1;
		}
		
		if (index_field(i - 1, j, k) > -1)
		{
			if (index_field(i, j, k) < i_res*j_res)
			{
				if (index_field(i, j, k) == 1)
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - 1]]*(A(index_field(i - 1, j, k), index_field(i, j, k)));
				}
				else if (index_field(i, j, k) < index_field.grid.i_res)
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k) - 1] + 1]*(A(index_field(i - 1, j, k), index_field(i, j, k)));
				}
				else
				{
					coef = (T)1/L.values[L.row_ptr[index_field(i, j, k)] - 1]*(A(index_field(i - 1, j, k), index_field(i, j, k)));
				}
			}
			else
			{
				coef = (T)1/L.values[L.row_ptr[index_field(i, j, k)] - 1]*(A(index_field(i - 1, j, k), index_field(i, j, k)));
			}

			L.AssignValue(index_field(i, j, k), index_field(i - 1, j, k), coef);
			number += 1;
		}

		if (index_field(i, j, k) == 0)
		{
			coef = sqrt(A(index_field(i, j, k), index_field(i, j, k)));
		}
		else
		{
			sum = 0;

			for (int t = L.row_ptr[index_field(i, j, k)] ; t  < number; t++)
			{
				sum += POW2(L.values[t]);
			}
			coef = sqrt(A(index_field(i, j, k), index_field(i, j, k)) - sum);
		}
		L.AssignValue(index_field(i, j, k), index_field(i, j, k), coef);
		
		number += 1;
	}
}

template<class TT>
static void DiagonalPreconditioner(MULTITHREADING* multithreading, const int& thread_id, const CSR_MATRIX<TT>& A, CSR_MATRIX<TT>& D)
{
	const int N = A.N;

	BEGIN_HEAD_THREAD_WORK
	{
		D.Initialize(multithreading, N, N);
		
		multithreading->SplitDomainIndex1D(0, D.N);
		
		D.start_ix[0] = 0;
		D.end_ix[0] = multithreading->sync_value_int[0] - 1;
		D.prev_row_array[0] = -1;
		D.values_ix_array[0] = 0;
		for (int id = 1; id < multithreading->num_threads; id++)
		{
			D.start_ix[id] = D.end_ix[id - 1] + 1;
			D.end_ix[id] = D.end_ix[id - 1] + multithreading->sync_value_int[id];
			D.prev_row_array[id] = -1;
			D.values_ix_array[id] = D.start_ix[id];
		}
	}
	END_HEAD_THREAD_WORK;

	const int start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

	for (int i = start_ix; i <= end_ix; i++)
	{
		D.row_ptr[i] = i;
	}
	multithreading->Sync(thread_id);

	for (int i = start_ix; i <= end_ix; i++)
	{
		for (int vix = D.row_ptr[i]; vix < D.row_ptr[i + 1]; vix++)
		{
			D.column_index[vix] = vix;
		}
	}
	multithreading->Sync(thread_id);

	for (int i = start_ix; i <= end_ix; i++)
	{
		D(i, i) = A(i, i);
	}
	multithreading->Sync(thread_id);
}



