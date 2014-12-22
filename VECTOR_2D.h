#pragma once

#include <iostream>

template<class TT>
class VECTOR_2D
{
public: // Essential member variables
	union
	{
		struct{TT x, y;};
		struct{TT i, j;};
		TT values[2];
	};	// Zero based indexing

public: // Constructors and Destructor
	inline VECTOR_2D(void)
		: x(TT()), y(TT())
	{}

	inline VECTOR_2D(const TT& x_input, const TT& y_input)
		: x(TT()), y(TT())
	{
		Initialize(x_input, y_input);
	}

	inline VECTOR_2D(const VECTOR_2D<TT>& vector_2d_input)
		: x(TT()), y(TT())
	{
		Initialize(vector_2d_input);
	}

	inline VECTOR_2D(const TT values_input[2])
		: x(TT()), y(TT())
	{
		Initialize(values_input);
	}

	~VECTOR_2D(void)
	{}
	
public: // Initialization Function
	void Initialize(const TT& x_input, const TT& y_input)
	{
		x = x_input;
		y = y_input;
	}

	void Initialize(const VECTOR_2D<TT>& vector_2d_input)
	{
		x = vector_2d_input.x;
		y = vector_2d_input.y;
	}

	void Initialize(const TT values_input[2])
	{
		x = values_input[0];
		y = values_input[1];
	}

public: // Operation overloading
	inline void operator = (const VECTOR_2D& v)
	{
		x = v.x;
		y = v.y;
	}

	inline void operator += (const VECTOR_2D& v)
	{
		x += v.x;
		y += v.y;
	}

	inline void operator -= (const VECTOR_2D& v)
	{
		x -= v.x;
		y -= v.y;
	}

	inline void operator *= (const TT& s)
	{
		x *= s;
		y *= s;
	}

	inline void operator /= (const TT& s)
	{
		TT one_over_s = (TT)1/s;
		x *= one_over_s;
		y *= one_over_s;
	}

	inline VECTOR_2D operator + (const VECTOR_2D& v) const
	{
		return VECTOR_2D(x + v.x, y + v.y);
	}

	inline VECTOR_2D operator - (const VECTOR_2D& v) const
	{
		return VECTOR_2D(x - v.x, y - v.y);
	}

	inline VECTOR_2D operator * (const TT& s) const
	{
		return VECTOR_2D(x*s, y*s);
	}

	inline VECTOR_2D operator / (const TT& s) const
	{
		TT one_over_s = (TT)1/s;
	
		return VECTOR_2D(x*one_over_s, y*one_over_s);
	}

public: // Member Functions
	inline TT Magnitude() const
	{
		return sqrt(x*x + y*y);
	}

	inline TT SqrMagnitude() const
	{
		return x*x + y*y;
	}

	inline void Normalize()
	{
		TT M = Magnitude();
		
		if(M != 0)
		{
			TT one_over_M = (TT)1/M;
			
			x *= one_over_M;
			y *= one_over_M;
		}
	}

	inline bool IsSqrMagnitudeSmallerThan(const TT& sqrmagnitude) const
	{
		if(sqrmagnitude > (x*x + y*y)) 
			return true;
		else
			return false;
	}

	// Be caution let normal vector be a unit vector
	inline void ScalingComponents(const VECTOR_2D<TT>& normal, const TT& normal_coef, const TT& tangential_coef)
	{
		const TT alpha = DotProduct(*this, normal);
		(*this) -= alpha*normal;
		(*this) = (normal_coef*alpha)*normal + tangential_coef*(*this);
	}
	
	inline void Assign(const TT& x_input, const TT& y_input)
	{
		x = x_input;
		y = y_input;
	}

	inline void AssignZeroVector()
	{
		x = (TT)0;
		y = (TT)0;
	}

	inline void AssignDifference(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
	{
		x = v1.x - v2.x;
		y = v1.y - v2.y;
	}

	inline void AssignDifferencePlusScaledDifference(const VECTOR_2D<TT>& p1, const VECTOR_2D<TT>& p2, const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2, const TT& dt)
	{
		x = (p1.x - p2.x) + (v1.x - v2.x)*dt;
		y = (p1.y - p2.y) + (v1.y - v2.y)*dt;
	}

	inline void AssignScaledDifference(const TT& scalar, const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
	{
		x = scalar*(v1.x - v2.x);
		y = scalar*(v1.y - v2.y);
	}

	inline void AssignScaledVector(const TT& scalar, const VECTOR_2D<TT>& v1)
	{
		x = scalar*v1.x;
		y = scalar*v1.y;
	}

	inline void AddSum(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
	{
		x += (v1.x + v2.x);
		y += (v1.y + v2.y);
	}

	inline void SubtractSum(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
	{
		x -= (v1.x + v2.x);
		y -= (v1.y + v2.y);
	}
};

// Miscellaneous Free Operators and Functions
template<class TT> 
inline static VECTOR_2D<TT> operator*(const TT& a, const VECTOR_2D<TT>& v)
{
	return VECTOR_2D<TT>(a*v.x, a*v.y);
}

template<class TT>
inline static TT DotProduct(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
{
	return v1.x*v2.x + v1.y*v2.y;
}

template<class TT>
inline static TT CrossProduct(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
{
	return v1.x*v2.y - v2.x*v1.y;
}

template<class TT>
inline static bool IsSqrDistanceSmallerThan(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2, const TT& sqrmagnitude)
{
	const TT diff_x(v1.x - v2.x), diff_y(v1.y - v2.y);
	if(sqrmagnitude > (diff_x*diff_x + diff_y*diff_y))
		return true;
	else 
		return false;
}

template<class TT>
inline std::ostream& operator << (std::ostream& output, const VECTOR_2D<TT>& v)
{
	return output << v.x << " " << v.y;
}







