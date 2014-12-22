#pragma once

#include "COMMON_DEFINITIONS.h"
#include "MATRIX_3X3.h"

class QUATERNION
{
public: // Essential Data
	T  s;
	VT v;

public: // Constructor and Destructor
	QUATERNION(void)
		: s(1), v(0, 0, 0)
	{}

	QUATERNION(const T& angle_input, const VT& direction_input)
		: s(cos((T)0.5*angle_input)), v(direction_input)
	{
		v.Normalize();
		v *= sin((T)0.5*angle_input);
	}

	QUATERNION(const T s_input, const T x_input, const T y_input, const T z_input)
		: s(s_input), v(x_input, y_input, z_input)
	{}

	QUATERNION(const VT& v_input)		// Make a vector into a quaternion (Note difference from From_Rotation_Vector)
		: s(0), v(v_input)
	{}

	QUATERNION(MATRIX_3X3& A)
	{
		T trace = 1 + A(1,1) + A(2,2) + A(3,3);
		if (trace > 1)
		{
			s = (T)0.5*sqrt(trace);
			v.x = A(3,2) - A(2,3);
			v.y = A(1,3) - A(3,1);
			v.z = A(2,1) - A(1,2);
			v *= (T)0.25/s;
		}
		else
		{
			int i = (A(1,1) > A(2,2)) ? 1: 2;
			i = (A(i,i) > A(3,3)) ? i: 3;					// Set i to be the index of the dominating diagonal term
			switch(i)
			{
			case 1: v.x = (T)0.5*sqrt(1 + A(1,1) - A(2,2) - A(3,3));
					v.y = (T)0.25*(A(2,1) + A(1,2))/v.x;
					v.z = (T)0.25*(A(1,3) + A(3,1))/v.x;
					s	= (T)0.25*(A(3,2) - A(2,3))/v.x;
					break;
			case 2: v.y = (T)0.5*sqrt(1 - A(1,1) + A(2,2) - A(3,3));
					v.x = (T)0.25*(A(2,1) + A(1,2))/v.y;
					v.z = (T)0.25*(A(3,2) + A(2,3))/v.y;
					s	= (T)0.25*(A(1,3) + A(3,1))/v.y; 
					break;
			case 3: v.z = (T)0.5*sqrt(1 - A(1,1) - A(2,2) + A(3,3));
					v.x = (T)0.25*(A(1,3) + A(3,1))/v.z;
					v.y = (T)0.25*(A(3,2) + A(2,3))/v.z;
					s	= (T)0.25*(A(2,1) + A(1,2))/v.z;
					break;
			}
		}
	}
	
	~QUATERNION(void)
	{}
	
public: // Operator Overloading
	QUATERNION operator*(const QUATERNION& q) const
	{
		VT r = s*q.v + q.s*v + CrossProduct(v, q.v);
		return QUATERNION(s*q.s - DotProduct(v, q.v), r.x, r.y, r.z);
	}

public: // Member Functions
	MATRIX_3X3 Matrix3X3() const
	{
		T vx2 = POW2(v.x), vy2 = POW2(v.y), vz2 = POW2(v.z), vxvy = v.x*v.y, vxvz = v.x*v.z, vyvz = v.y*v.z, svx = s*v.x, svy = s*v.y, svz = s*v.z;

		return MATRIX_3X3(1 - 2*(vy2 + vz2), 2*(vxvy + svz)   , 2*(vxvz - svy)   ,
						  2*(vxvy - svz)   , 1 - 2*(vx2 + vz2), 2*(vyvz + svx)   ,
						  2*(vxvz + svy)   , 2*(vyvz - svx)   , 1 - 2*(vx2 + vy2)
						  );
	}

	T Magnitude() const
	{
		return sqrt(POW2(s) + POW2(v.x) + POW2(v.y) + POW2(v.z));
	}
	
	void Normalize()
	{
		T magnitude = Magnitude();
		assert(magnitude != 0);

		T r = (T)1/magnitude;
		
		s *= r;
		v.x *= r;
		v.y *= r;
		v.z *= r;
	}

	QUATERNION Inverse() const
	{
		return QUATERNION(s, -v.x, -v.y, -v.z);
	}
	
	VT Rotate(const VT& v) const
	{
		QUATERNION q = *this*QUATERNION(v)*Inverse();
		return q.v;
	}
	
	VT Inverse_Rotate(const VT& v) const
	{
		QUATERNION q = Inverse()*QUATERNION(v)*(*this);
		return q.v;
	}

	static QUATERNION FromRotationVector(const VT& v)
	{
		const T magnitude = v.Magnitude();
		if(magnitude <= (T)1e-8)
		{
			return QUATERNION();
		}

		QUATERNION q;
		q.s = cos((T)0.5*magnitude);
		q.v = ((T)sin((T)0.5*magnitude)/magnitude)*v;
		return q;
	}

	static QUATERNION Slerp(const QUATERNION& a_input, const QUATERNION& b_input, T t)
	{
		QUATERNION a(a_input);
		QUATERNION b(b_input);

		// T sina, sinat, sinaomt, angle
		T fAlpha, fTheta, fBeta, fCosTheta, oosinTheta;

		fCosTheta = a.v.x*b.v.x + a.v.y*b.v.y + a.v.z*b.v.z + a.s*b.s;

		if (fCosTheta < (T)0)
		{
			b.v *= (T)-1;
			b.s *= (T)-1;
			fCosTheta *= (T)-1;
		}

		fAlpha = t;
		fBeta = (T)1 - fAlpha;

		if ((T)1 - fCosTheta > (T)0.001)
		{
			fTheta = acos(fCosTheta);
			oosinTheta = (T)1/sin(fTheta);
			fBeta = sin(fTheta*fBeta)*oosinTheta;
			fAlpha = sin(fTheta*fAlpha)*oosinTheta;
		}

		return QUATERNION(fBeta*a.s + fBeta*b.s, fBeta*a.v.x + fBeta*b.v.x, fBeta*a.v.y + fBeta*b.v.y, fBeta*a.v.z + fBeta*b.v.z);
	}
};



