#pragma once

#include "COMMON_DEFINITIONS.h"

// Coefficients for Gauss-Seidel iteration
class GS_COEFFICIENTS
{
public: // Essential Data
	T plus[3], minus[3];
	T center, inv_center;

public: // Constructor and Destructor
	GS_COEFFICIENTS(void)
	{
		Initialize();
	}

public: // Initialization Functions
	inline void Initialize(void)
	{
		plus[0] = (T)0;
		plus[1] = (T)0;
		plus[2] = (T)0;
		minus[0] = (T)0;
		minus[1] = (T)0;
		minus[2] = (T)0;
		center = (T)0;
		inv_center = (T)1;
	}

public: // Operator Overloading
	inline void operator=(const GS_COEFFICIENTS& v)
	{
		for (int i = 0; i < 3; i++)
		{
			plus[i] = v.plus[i];
			minus[i] = v.minus[i];
		}
		center = (T)0;
		inv_center = (T)1;
	}

	void Cout()
	{
		cout << center << " " << inv_center << " " << plus[0] << " " << plus[1] << " " << plus[2] << " " << minus[0] << " " << minus[1] << " " << minus[2] << endl;
	}
};
