#pragma once

#include "VECTOR_ND.h"
#include "MATRIX_NXN.h"

template<class T>
class MATRIX_MXN
{
public: // Essential Data
	int m, n;			// Size of matrix
	T* x;				// Entries in matrix
	MATRIX_MXN<T> *Q, *V;
	MATRIX_NXN<T> *R;

public: // Constructors and Destructor
	MATRIX_MXN(void)
		: m(0), n(0), x(0), Q(0), V(0), R(0)
	{}

	MATRIX_MXN(const int& m_input, const int& n_input)
		: m(0), n(0), x(0), Q(0), V(0), R(0)
	{
		Initialize(m_input, n_input);
	}

	MATRIX_MXN(const MATRIX_MXN<T>& A)
		: m(0), n(0), x(0), Q(0), V(0), R(0)
	{
		Initialize(A);
	}

	~MATRIX_MXN(void)
	{
		if(x != 0) delete [] x;
		if(Q != 0) delete Q;
		if(V != 0) delete V;
		if(R != 0) delete R;
	}

public: // Initialization Functions
	void Initialize(const int& m_input, const int& n_input)
	{
		m = m_input; 
		n = n_input;
		
		if(x != 0) delete [] x;
		x = new T[m*n];
		
		for(int i = 0; i < m*n; i++) x[i] = 0;
	}

	void Initialize(const MATRIX_MXN<T>& A)
	{
		m = A.m;
		n = A.n;

		if(x != 0) delete [] x;
		x = new T[m*n];

		for(int i = 0; i < m*n; i++) x[i] = A.x[i];
	}

public: // Operator Overloading
	T& operator()(const int& i, const int& j)
	{
		assert(i >= 0 && i < m);
		assert(j >= 0 && j < m);

		return *(x + j*m + i);
	}

	const T& operator()(const int& i, const int& j) const
	{
		assert(i >= 0 && i < m);
		assert(j >= 0 && j < m);

		return *(x + j*m + i);
	}

	MATRIX_MXN<T>& operator=(const T& a)
	{
		int size = m*n;

		for(int i = 0; i < size; i++) x[i] = a;
		
		return *this;
	}
	
	MATRIX_MXN<T>& operator=(const MATRIX_MXN<T>& A)
	{
		delete Q;
		Q = 0;	
		
		delete V;
		V = 0;
		
		delete R;
		R = 0;

		if(x != 0 || m*n != A.m*A.n) 
		{
			delete [] x;
			T* x = new T[A.m*A.n];
		}

		int size = A.m*A.n;
		
		for(int i = 0; i < size; i++) x[i] = A[i];

		return *this;
	}

	MATRIX_MXN<T>& operator*=(const T& a)
	{
		int size = m*n;

		for(int i = 0; i < size; i++) x[i] *= a;
		
		return *this;
	}

	MATRIX_MXN<T>& operator+=(const MATRIX_MXN<T>& A)
	{
		assert(m == A.m && n == A.n);

		int size = m*n;

		for(int i = 0; i < size; i++) x[i] += A.x[i];

		return *this;
	}

	MATRIX_MXN<T> operator-(const MATRIX_MXN<T>& A) const
	{
		assert(m == A.m && n == A.n);

		MATRIX_MXN<T> matrix(m, n);
		
		int size = m*n;

		for(int i = 0; i < size; i++) matrix.x[i] = x[i] - A.x[i];

		return matrix;
	}

	MATRIX_MXN<T> operator-() const
	{
		MATRIX_MXN<T> matrix(m, n);

		int size = m*n;

		for(int i = 0; i < size; i++) matrix.x[i] = -x[i];

		return matrix;
	}

	VECTOR_ND<T> operator*(const VECTOR_ND<T>& y) const
	{
		assert(y.num_dimension == n);

		VECTOR_ND<T> result(m);

		for(int j = 0; j < m; j++)
		{
			for(int i = 0; i < n; i++)
			{
				result[i] += *(x + j*m + i)*y[j];
			}
		}

		return result;
	}

	MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
	{
		assert(n == A.m);

		MATRIX_MXN<T> matrix(m, A.n);
		
		for(int j = 0; j < A.n; j++)
		{
			for(int i = 0; i < m; i++)
			{
				for(int k = 0; k < n; k++)
				{
					matrix(i,j) += *(x + k*m + i)*A(k,j);
				}
			}
		}

		return matrix;
	}

	MATRIX_MXN<T> operator*(const MATRIX_NXN<T>& A) const
	{
		assert(n == A.n);

		MATRIX_MXN<T> matrix(m, A.n);

		for(int j = 0; j < A.n; j++)
		{
			for(int i = 0; i < m; i++)
			{
				for(int k = 0; k < n; k++)
				{
					matrix(i,j) += *(x + k*m + i)*A(k,j);
				}
			}
		}

		return matrix;
	}

public: // Member Functions
	void Resize(const int& m_new, const int& n_new, const bool copy_exciting_elements = true)
	{
		if(m_new == m && n == n_new) return;

		T* x_new = new T[m_new*n_new];

		int new_size = m_new*n_new;

		for(int i = 0; i < new_size; i++) x_new[i] = (T)0;

		int m1 = min(m, m_new), n1 = min(n, n_new);

		for(int i = 0; i < m1; i++)
		{
			for(int j = 0; j < n1; j++)
			{
				*(x_new + j*m_new + i) = (*this)(i,j);
			}
		}

		delete [] x;
		x = x_new; m = m_new; n = n_new;
				
		delete Q; Q = 0;
		delete V; V = 0;
		delete R; R = 0;
	}

	inline void Transposed_Multiply(const VECTOR_ND<T>& y, VECTOR_ND<T>& result)
	{
		assert(y.num_dimension == m && result.num_dimension == n);

		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				result[j] += *(x + j*m + i)*y[i];
			}
		}
	}

	inline VECTOR_ND<T> Transposed_Multiply(const VECTOR_ND<T>& y) const
	{
		assert(y.num_dimension == m);
		
		VECTOR_ND<T> result;
		result.Initialize(n, true);

		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				result.values[j] += *(x + j*m + i)*y.values[i];
			}
		}

		return result;
	}

	// Need to check this! 
	inline void Multiply(const VECTOR_ND<T>& y, VECTOR_ND<T>& result)
	{
		assert(y.num_dimension == n && result.num_dimension == m);

		for(j = 0; j < n; j++)
		{
			for(i = 0; i < m; i++)
			{
				result[i] += (x + j*n + i)*y[j];
			}
		}
	}

	MATRIX_MXN<T> Transpose()
	{
		MATRIX_MXN<T> matrix(n, m);

		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				matrix(j,i) = *(x + j*m + i);
			}
		}

		return matrix;
	}

	VECTOR_ND<T> Transpose_Times(const VECTOR_ND<T>& y) const
	{
		assert(y.num_dimension == m);

		VECTOR_ND<T> result(n, true);

		for(int j = 0; j < n; j++)
		{
			for(int i = 0; i < m; i++)
			{
				result[j] += *(x + j*m + i)*y[i];
			}
		}

		return result;
	}

	MATRIX_MXN<T> Transpose_Times(const MATRIX_MXN<T>& A) const
	{
		assert(m == A.m);
		
		MATRIX_MXN<T> matrix(n, A.n);

		for(int j = 0; j < A.n; j++)
		{
			for(int i = 0; i < n; i++)
			{
				for(int k = 0; k < m; k++)
				{
					matrix(i,j) += *(x + i*m + k)*A(k,j);
				}
			}
		}

		return matrix;
	}

	MATRIX_MXN<T> Permute_Columns(const VECTOR_ND<int>& p) const
	{
		assert(n == p.num_dimension);

		MATRIX_MXN<T> matrix(m, n);

		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				matrix(i,j) = (*this)(i,p(j));
			}
		}

		return matrix;
	}

	MATRIX_MXN<T> Unpermute_Columns(const VECTOR_ND<int>& p) const
	{
		assert(n == p.num_dimension);

		MATRIX_MXN<T> matrix(m, n);

		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++)
			{
				matrix(i,p(j)) = (*this)(i,j);
			}
		}

		return matrix;
	}

	void Set_Zero_Matrix()
	{
		int size = m*n;

		for(int i = 0; i < size; i++) x[i] = 0;
	}

	T Frobenius_Norm_Squared() const
	{
		T sum(0);

		int size = m*n;

		for(int i = 0; i < size; i++)
		{
			sum += sqr(x[i]);
		}

		return sum;
	}

	MATRIX_NXN<T> Normal_Equations_Matrix() const
	{
		MATRIX_NXN<T> result(n, true);

		for(int j = 0; j < n; j++)
		{
			for(int i = 0; i < n; i++)
			{
				for(int k = 0; k < m; k++)
				{
					result(i,j) += *(x + i*m + k)**(x + j*m + k);
				}
			}
		}

		return result;
	}

	MATRIX_NXN<T> Weighted_Normal_Equations_Matrix(const VECTOR_ND<T>& w) const
	{
		MATRIX_NXN<T> result(n);
		for(int j = 0; j < n; j++)
		{
			for(int i = 0; i < n; i++)
			{
				for(int k = 0; k < m; k++)
				{
					result(i,j) += *(x + i*m + k)**(x + j*m + k)*w.values[k]*w.values[k];
				}
			}
		}

		return result;
	}

	VECTOR_ND<T> Normal_Equations_Solve(const VECTOR_ND<T>& b) const
	{
		MATRIX<T> A_transpose_A(Normal_Equations_Matrix(A));
		VECTOR_ND<T> A_transpose_b(Transpose_Times(b));

		return A_transpose_A.Cholesky_Solve(A_transpose_b);
	}

	void Set_Column(const int& j, const VECTOR_ND<T>& a)
	{
		assert(a.num_dimension == m);
		for(int i = 0; i < m; i++) (*this)(i,j) = a[i];
	}

	void Set_Row(const int& i, const VECTOR_ND<T>& a)
	{
		assert(a.num_dimension == n);
		for(int j = 0; j < n; j++) (*this)(i,j) = a[j];
	}

	void Get_Column(const int& j, VECTOR_ND<T>& a) const
	{
		assert(a.num_dimension == m);
		for(int i = 0; i < m; i++) a[i] = (*this)(i,j);
	}

	VECTOR_ND<T> Get_Column(const int& j) const
	{
		VECTOR_ND<T> a(m);
		for(int i = 0; i < m; i++) a[i] = (*this)(i,j);

		return a;
	}

	void Add_To_Submatrix(const int& i_start, const int& j_start, const VECTOR_ND<T>& a)
	{
		for(int i = 0; i < a.num_dimension; i++)
		{
			(*this)(i_start + i, j_start) += a[i];
		}
	}

	void Add_To_Submatrix(const int& i_start, const int& j_start, const MATRIX_MXN<T>& A)
	{
		for(int i = 0; i < A.m; i++)
		{
			for(int j = 0; j < A.n; j++)
			{
				(*this)(i_start + i, j_start + j) += A(i,j);
			}
		}
	}

	void Add_To_Submatrix(const int& i_start, const int& j_start, const MATRIX_NXN<T>& A)
	{
		for(int i = 0; i < A.n; i++)
		{
			for(int j = 0; j < A.n; j++)
			{
				(*this)(i_start + i, j_start + j) += A(i,j);
			}
		}
	}
};

template<class T> inline ostream& operator<<(ostream& output_stream, const MATRIX_MXN<T>& A)
{
	for(int i = 0; i < A.m; i++)
	{
		for(int j = 0; j < A.n; j++)
		{
			output_stream << A(i,j) << " " ;
			output_stream << endl;
		}
	}

	return output_stream;
}

template<class T> inline MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A, const MATRIX_MXN<T>& B)
{
	assert(A.n == B.m);
	
	MATRIX_MXN<T> matrix(A.m, B.n);

	for(int j = 0; j < B.n; j++)
	{
		for(int i = 0; i < A.m; i++)
		{
			for(int k = 0; k < A.n; k++)
			{
				matrix(i,j) += A(i,k)*B(k,j);
			}
		}
	}

	return matrix;
}

// Need to check this header one more time!











		






