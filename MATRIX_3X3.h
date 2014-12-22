#pragma once

#include "COMMON_DEFINITIONS.h"

class MATRIX_3X3
{
public: // Essential Data
	T x[9];

public: // Constructors and Destructor
	MATRIX_3X3(void)
	{}

	MATRIX_3X3(const T x11, const T x21, const T x31, const T x12, const T x22, const T x32, const T x13, const T x23, const T x33)
	{
		Initialize(x11, x21, x31, x12, x22, x32, x13, x23, x33);
	}

	~MATRIX_3X3(void)
	{}
	
public: // Initialization Function
	void Initialize(const T x11, const T x21, const T x31, const T x12, const T x22, const T x32, const T x13, const T x23, const T x33)
	{
		x[0] = x11;
		x[1] = x21;
		x[2] = x31;
		x[3] = x12;
		x[4] = x22;
		x[5] = x32;
		x[6] = x13;
		x[7] = x23;
		x[8] = x33;
	}

public: // Member Functions
	T& operator()(const int& i, const int& j)
	{
		assert(i >= 1 && i <= 3);
		assert(j >= 1 && j <= 3);
		return x[i - 1 + 3*(j - 1)];
	}

	MATRIX_3X3 operator*(const MATRIX_3X3& A) const
	{
		return MATRIX_3X3(
			x[0]*A.x[0] + x[3]*A.x[1] + x[6]*A.x[2], x[1]*A.x[0] + x[4]*A.x[1] + x[7]*A.x[2], x[2]*A.x[0] + x[5]*A.x[1] + x[7]*A.x[2],
			x[0]*A.x[3] + x[3]*A.x[4] + x[6]*A.x[5], x[1]*A.x[3] + x[4]*A.x[4] + x[7]*A.x[5], x[2]*A.x[3] + x[5]*A.x[4] + x[7]*A.x[5],
			x[0]*A.x[6] + x[3]*A.x[7] + x[6]*A.x[8], x[1]*A.x[6] + x[4]*A.x[7] + x[7]*A.x[8], x[2]*A.x[6] + x[5]*A.x[7] + x[7]*A.x[8]
		);
	}

	MATRIX_3X3 Transposed() const
	{
		return MATRIX_3X3(x[0], x[3], x[6], x[1], x[4], x[7], x[2], x[5], x[8]);
	}

	MATRIX_3X3 Inversed()
	{
		T cofactor11 = x[4]*x[8] - x[7]*x[5], cofactor12 = x[7]*x[2] - x[1]*x[8], cofactor13 = x[1]*x[5] - x[4]*x[2];
		T determinant = x[0]*cofactor11 + x[3]*cofactor12 + x[6]*cofactor13;

		assert(determinant != 0);

		T s = (T)1/determinant;

		return MATRIX_3X3(s*cofactor11             , s*cofactor12             ,	s*cofactor13             ,
						  s*(x[6]*x[5] - x[3]*x[8]), s*(x[0]*x[8] - x[6]*x[2]), s*(x[3]*x[2] - x[0]*x[5]),
						  s*(x[3]*x[7] - x[6]*x[4]), s*(x[6]*x[1] - x[0]*x[7]), s*(x[0]*x[4] - x[3]*x[1])
						  );
	}

	MATRIX_3X3 operator*(const T a) const
	{
		return MATRIX_3X3(a*x[0], a*x[1], a*x[2], a*x[3], a*x[4], a*x[5], a*x[6], a*x[7], a*x[8]);
	}

	static MATRIX_3X3 Identity()
	{
		return MATRIX_3X3(1, 0 , 0, 0, 1, 0, 0, 0, 1);
	}

	void operator+=(const MATRIX_3X3& A)
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] += A.x[i];
		}
	}

	void operator-=(const MATRIX_3X3& A)
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] -= A.x[i];
		}
	}
	
	void operator*=(const T& a) 	
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] *= a;
		}
	}

	void operator=(const MATRIX_3X3& A)
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] = A.x[i];
		}
	}
};

inline VT operator*(const MATRIX_3X3& A, const VT& v)
{
	return VT(v.x*A.x[0] + v.y*A.x[3] + v.z*A.x[6], v.x*A.x[1] + v.y*A.x[4] + v.z*A.x[7], v.x*A.x[2] + v.y*A.x[5] + v.z*A.x[8]);
}

inline MATRIX_3X3 operator*(const T& s, const MATRIX_3X3& A)
{
	return MATRIX_3X3(s*A.x[0], s*A.x[1], s*A.x[2], s*A.x[3], s*A.x[4], s*A.x[5], s*A.x[6], s*A.x[7], s*A.x[8]);
}

