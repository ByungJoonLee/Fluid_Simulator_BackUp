#pragma once

#include "VECTOR_ND.h"
#include "MATRIX_MXN.h"

#ifndef __exchange__
#define __exchange__
template<class T> inline void exchange(T& a, T& b)
{
	T c = a; a = b; b = c;
}
#endif

template<class T>
class MATRIX_NXN
{
public: // Essential Data
	int n;			// Size of matrix
	T* x;				// Entries in matrix
	T small_number;
	MATRIX_NXN<T> *L, *U, *inverse;
	VECTOR_ND<int>* p;
	
public: // Constructors and Destructor
	MATRIX_NXN(void)
		: n(0), x(0), small_number((T)1e-8), L(0), U(0), inverse(0), p(0)
	{}

	MATRIX_NXN(const int& n_input)
		: n(0), small_number((T)1e-8), L(0), U(0), inverse(0), p(0)
	{
		Initialize(n_input);
	}

	MATRIX_NXN(const MATRIX_NXN<T>& A)
		: n(0), small_number((T)1e-8), L(0), U(0), inverse(0), p(0)
	{
		Initialize(A);
	}

	~MATRIX_NXN(void)
	{
		if(x != 0) delete [] x;
		if(L != 0) delete L;
		if(U != 0) delete U;
		if(p != 0) delete p;
		if(inverse != 0) delete inverse;
	}

public: // Initialization Functions
	void Initialize(const int& n_input)
	{
		n = n_input;
		
		if(x != 0) delete [] x;
		x = new T[n*n];
		
		for(int i = 0; i < n*n; i++) x[i] = 0;
	}

	void Initialize(const MATRIX_NXN<T>& A)
	{
		n = A.n;

		if(x != 0) delete [] x;
		x = new T[n*n];

		for(int i = 0; i < n*n; i++) x[i] = A.x[i];
	}

public: // Operator Overloading
	T& operator()(const int& i, const int& j) 
	{
		assert(i >= 0 && i < n);
		assert(j >= 0 && j < n);

		return *(x + i*n + j);
	}
	
	const T& operator()(const int& i, const int& j) const
	{
		assert(i >= 0 && i < n);
		assert(j >= 0 && j < n);

		return *(x + i*n + j);
	}

	bool operator==(const MATRIX_NXN<T>& A) const
	{
		assert(n == A.n);
		
		for(int i = 0; i < n*n; i++)
		{
			if(x[i] != A.x[i])
				return false;
		}

		return true;
	}

	bool operator!=(const MATRIX_NXN<T>& A) const
	{
		assert(n == A.n);

		for(int i = 0; i < n*n; i++)
		{
			if(x[i] != A.x[i])
				return true;
		}

		return false;
	}

	MATRIX_NXN<T>& operator=(const MATRIX_NXN<T>& A)
	{
		assert(n == A.n);

		delete L;
		L = 0;

		delete U;
		U = 0;

		delete inverse;
		inverse = 0;

		delete p;
		p = 0;
		
		if(x != 0 || n != A.n)
		{
			delete [] x;
			x = new T[A.n*A.n];
		}

		int size = A.n*A.n;
		
		for(int i = 0; i < size; i++) x[i] = A.x[i];

		return *this;
	}

	MATRIX_NXN<T>& operator+=(const MATRIX_NXN<T>& A)
	{
		assert(n == A.n);

		for(int i = 0; i < n*n; i++) x[i] += A.x[i];
		return *this;
	}

	MATRIX_NXN<T>& operator-=(const MATRIX_NXN<T>& A)
	{
		assert(n == A.n);

		for(int i = 0; i < n*n; i++) x[i] -= A.x[i];
		return *this;
	}

	MATRIX_NXN<T>& operator=(const T& a)
	{
		for(int i = 0; i < n*n; i++) x[i] = a;
		return *this;
	}

	MATRIX_NXN<T>& operator*=(const T& a)
	{
		for(int i = 0; i < n*n; i++) x[i] *= a;
		return *this;
	}

	MATRIX_NXN<T> operator+(const MATRIX_NXN<T>& A) const
	{
		assert(A.n == n);
		
		MATRIX_NXN<T> result(n, true);
		for(int i = 0; i < n*n; i++) result.x[i] = x[i] + A.x[i];

		return result;
	}

	MATRIX_NXN<T> operator-(const MATRIX_NXN<T>& A) const
	{
		assert(A.n == n);

		MATRIX_NXN<T> result(n, true);
		for(int i = 0; i < n*n; i++) result.x[i] = x[i] - A.x[i];
		
		return result;
	}

	MATRIX_NXN<T> operator*(const T& a) const
	{
		MATRIX_NXN<T> result(n, true);
		for(int i = 0; i < n*n; i++) result.x[i] = a*x[i];

		return result;
	}

	VECTOR_ND<T> operator*(const VECTOR_ND<T>& y) const
	{
		VECTOR_ND<T> result(n, true);
		Multiply(y, result);
		
		return result;
	}

	MATRIX_NXN<T> operator*(const MATRIX_NXN<T>& A) const
	{
		assert(A.n == n);

		MATRIX_NXN<T> result(n);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				T total = (T)0;
				for(int k = 0; k < n; k++)
				{
					total += (*(x + i*n + k))*(*(A.x + k*n + j));
					(*(result.x + i*n + j)) = total;
				}
			}
		}

		return result;
	}

public: // Member Functions
	void Multiply(const VECTOR_ND<T>& y, VECTOR_ND<T>& result) const
	{
		assert(y.num_dimension == n && result.num_dimension ==n);

		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				(*(result.x + i)) += (*(x + i*n + j))*(*(y.x + j));
			}
		}
	}

	VECTOR_ND<T> Transpose_Times(const VECTOR_ND<T>& y) const
	{
		assert(y.num_dimension == n);

		VECTOR_ND<T> result(n);
		for(int j = 0; j < n; j++)
		{
			for(int i = 0; i < n; i++)
			{
				(*(result.x + i)) += (*(x + j*n + i))*(*(y.x + j));
			}
		}

		return result;
	}

	MATRIX_NXN<T> Transpose_Times(const MATRIX_NXN<T>& A) const
	{
		assert(A.n == n);

		MATRIX_NXN<T> matrix(n);
		for(int j = 0; j < A.n; j++)
		{
			for(int i = 0; i < n; i++)
			{
				for(int k = 0; k < n; k++)
				{
					matrix(i,j) += (*this)(k,i)*A(k,j);
				}
			}
		}

		return matrix;
	}

	MATRIX_NXN<T> Permute_Columns(const VECTOR_ND<int>& p) const
	{
		assert(n == p.num_dimension);

		MATRIX_NXN<T> matrix(n);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; k < n; j++)
			{
				matrix(i,j) = (*this)(i,p(j));
			}
		}

		return matrix;
	}

	MATRIX_NXN<T> Unpermute_Columns(const VECTOR_ND<int>& p) const
	{
		assert(n == p.num_dimension);

		MATRIX_NXN<T> matrix(n);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				matrix(i,p(j)) = (*this)(i,j);
			}
		}

		return matrix;
	}

	VECTOR_ND<T> Lower_Triangular_Solve(const VECTOR_ND<T>& b) const
	{
		assert(n == b.num_dimension);

		VECTOR_ND<T> v(n);
		for(int i = 0; i < n; i++)
		{
			v.values[i] = b.values[i];

			for(int j = 0; j <= i-1; j++)
			{
				v.values[i] -= (*this)(i,j)*v.values[j];
			}

			v.values[i]/=(*this)(i,i);
		}

		return v;
	}
	
	VECTOR_ND<T> Transpose_Lower_Triangular_Solve(const VECTOR_ND<T>& b) const
	{
		assert(n == b.num_dimension);

		VECTOR_ND<T> v(b);

		for(int i = 0; i < n; i++) 
		{
			v.values[i] /= (*this)(i,i);
			
			for(int j = i+1; j < n; j++)
			{
				v.values[j] -=(*this)(i,j)*v.values[i];
			}
		}

		return v;
	}

	VECTOR_ND<T> Upper_Triangular_Solve(const VECTOR_ND<T>& b) const
	{
		assert(n == b.num_dimension);
		
		VECTOR_ND<T> v(n);
		for(int i = 0; i < n; i++)
		{
			v.values[i] = b[i];

			for(int j = n-1; j >= i+1; j--)
			{
				v.values[i] -= (*this)(i,j)*v.values[j];
			}
		}

		return v;
	}

	VECTOR_ND<T> Transpose_Upper_Triangular_Solve(const VECTOR_ND<T>& b) const
	{
		assert(n == b.num_dimension);

		VECTOR_ND<T> v(b);
		for(int i = n; i >= 1; i--)
		{
			v.values[i] /= (*this)(i,i);

			for(int j = i-1; j >= 1; j--)
			{
				v.values[j] -= (*this)(i,j)*v.values[i];
			}
		}

		return v;
	}

	void LU_Factorization()
	{
		int i, j, k;
		
		if(L != 0) delete L;
		MATRIX_NXN<T> L = new MATRIX_NXN<T>(n);

		if(U != 0) delete U;
		MATRIX_NXN<T> U = new MATRIX_NXN<T>(*this);

		for(j = 0; j < n; j++)
		{
			for(i = j; i < n; i++)
			{
				(*L)(i,j) = (*U)(i,j)/(*U)(j,j);
			}
			for(i = j+1; i < n; i++)
			{
				for(k = j; k < n; k++)
				{
					(*U)(i,k) -= (*L)(i,j)*(*U)(j,k);
				}
			}
		}
	}

	VECTOR_ND<T> LU_Solve(const VECTOR_ND<T>& b)
	{
		LU_Factorization();

		return U->Upper_Triangular_Solve(L->Lower_Triangular_Solve(b));
	}

	void LU_Inverse()
	{
		if(inverse != 0) delete inverse;
		MATRIX_NXN<T> inverse = new MATRIX_NXN<T>(n);
		
		LU_Factorization();
	
		for(int j = 0; j < n; j++)
		{
			VECTOR_ND<T> b(n);
			b(j) = 1;
			
			VECTOR_ND<T> x = U->Upper_Triangular_Solve(L->Lower_Triangular_Solve(b));
			for(int i = 0; i < n; i++)
			{
				(*inverse)(i,j) = x.values[i];
			}
		}
	}

	void PLU_Factorization()
	{
		int i, j, k;
		if(L != 0) delete L;
		MATRIX_NXN<T> L = new MATRIX_NXN<T>(n);

		if(U != 0) delete U;
		MATRIX_NXN<T> U = new MATRIX_NXN<T>(*this);
		
		if(p != 0) delete p;
		VECTOR_ND<T> p = new VECTOR_ND<T>(n, false);

		for(int i = 0; i < n; i++) (*p)(i) = i;

		for(j = 0; j < n; j++)
		{
			int row = j;
			T value = abs((*U)(i,j));
			for(i = j+1; i < n; i++)
			{
				if(abs((*U)(i,j)) > value)
				{
					row = i;
					value = abs((*U)(i,j));
				}
			}
			if(row != j)
			{
				exchange((*p)(j), (*p)(row));
				for(k = 0; k < j-1; k++)
				{
					exchange((*L)(j,k), (*L)(row,k));
				}
				for(k = j; k < n; k++)
				{
					exchange((*U)(j,k), (*U)(row, k));
				}
			}
			for(i = j; i < n; i++) 
			{
				(*L)(i,j) = (*U)(i,j)/(*U)(j,i);
			}
			for(i = j+1; i < n; i++)
			{
				for(k = j; k < n; k++)
				{
					(*U)(i,k) -= (*L)(i,j)*(*U)(j,k);
				}
			}
		}
	}

	// Need to fix this!!
	VECTOR_ND<T> PLU_Solve(const VECTOR_ND<T>& b)
	{
		PLU_Factorization();
		VECTOR_ND<T> x = b;
		return U->Upper_Triangular_Solve(L->Lower_Triangular_Solve(x.Permute(*p)));
	}

	void Cholesky_Factorization()
	{
		int i, j, k;

		if(L != 0) delete L;
		L = new MATRIX_NXN<T>(n);

		if(U != 0) delete U;
		U = new MATRIX_NXN<T>(*this);

		for(j = 0; j < n; j++)
		{
			for(k = 0; k <= j-1; k++)
			{
				for(i = j; i < n; i++)
				{
					(*U)(i,j) -= (*L)(j,k)*(*L)(i,k);
				}
			}
			
			(*L)(j,i) = sqrt((*U)(j,i));

			for(i = j+1; i < n; i++)
			{
				(*L)(i,j) = (*U)(i,j)/(*L)(j,j);
			}
		}
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
			{
				(*U)(i,j) = (*L)(i,j);
			}
		}
	}

	VECTOR_ND<T> Cholesky_Solve(const VECTOR_ND<T>& b)
	{
		Cholesky_Factorization();
		return U->Upper_Triangular_Solve(L->Lower_Triangular_Solve(b));
	}

	void Cholesky_Inverse()
	{
		if(inverse != 0) delete inverse;
		MATRIX_NXN<T> inverse = new MATRIX_NXN<T>(n);

		Cholesky_Factorization();
		
		for(int j = 0; j < n; j++)
		{
			VECTOR_ND<T> b(n);
			b(j) = 1;
			
			VECTOR_ND<T> x = U->Upper_Triangular_Solve(L->Lower_Triangular_Solve(b));
			for(int i = 0; i < n; i++) 
			{
				(*inverse)(i,j) = x.values[i];
			}
		}
	}

	void Set_Column(const int& j, const VECTOR_ND<T>& a)
	{
		assert(a.num_dimension == n);
		
		for(int i = 0; i < n; i++)
		{
			(*this)(i,j) = a[i];
		}
	}

	void Get_Column(const int& j, const VECTOR_ND<T>& a) const
	{
		assert(a.num_dimension == n);

		for(int i = 0; i < n; i++)
		{
			a[i] = (*this)(i,j);
		}
	}

	VECTOR_ND<T> Get_Column(const int& j) const
	{
		VECTOR_ND<T> a(n);

		for(int i = 0; i < n; i++)
		{
			a[i] = (*this)(i,j);
		}

		return a;
	}

	void Add_To_Submatrix(const int& i_start, const int& j_start, const MATRIX_NXN<T>& a)
	{
		for(int i = 0; i < a.m; i++)
		{
			for(int j = 0; j < a.n; j++)
			{
				(*this)(i_start + i, j_start + j) += a(i,j);
			}
		}
	}
};

template<class T> inline MATRIX_NXN<T> operator*(const T& a, const MATRIX_NXN<T>& A)
{
	MATRIX_NXN<T> result(A.n);
	
	int size = A.n*A.n;

	for(int i = 0; i < size; i++)
	{
		result.x[i] = a*A[i];
	}
}

template<class T> inline MATRIX_NXN<T> operator*(const int& a, const MATRIX_NXN<T>& A)
{
	MATRIX_NXN<T> result(A.n);
	
	int size = A.n*A.n;

	for(int i = 0; i < size; i++)
	{
		result.x[i] = a*A[i];
	}
}

template<class T> inline istream& operator>>(istream& input_stream, MATRIX_NXN<T>& A)
{
	for(int i = 0; i < A.n; i++)
	{
		for(int j = 0; j < A.n; j++)
		{
			input_stream >> A(i,j);
		}
	}
	
	return input_stream;
}

template<class T> inline ostream& operator<<(ostream& output_stream, const MATRIX_NXN<T>& A)
{
	for(int i = 0; i < A.n; i++)
	{
		for(int j = 0; j < A.n; j++)
		{
			output_stream << A(i,j) << " ";
			output_stream << endl;
		}
	}

	return output_stream;
}





		






					

		

		




		

				

