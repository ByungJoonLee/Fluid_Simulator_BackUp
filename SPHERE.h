#pragma once

#include "COMMON_DEFINITIONS.h"
#include "LEVELSET_OBJECT.h"

class SPHERE : public LEVELSET_OBJECT
{
public: // Essential Data
	VT center;
	T radius;

public: // Constructors and Desturctor
	SPHERE(const VT& center_input, const T& radius_input)
	{
		Initialize(center_input, radius_input);
	}

	~SPHERE(void)
	{}

public: // Initialization Functions
	void Initialize(const VT& center_input, const T& radius_input)
	{
		center = center_input;
		radius = radius_input;
	}

public: // Member Functions
	inline const T SignedDistance(const VT& position) const
	{
		return (position - center).Magnitude() - radius;
	}

	inline void Normal(const VT& position, VT& normal) const
	{
		normal = position - center;
	}

	inline const VT Normal(const VT& position) const
	{
		return position - center;
	}

	inline void UnitNormal(const VT& position, VT& normal) const
	{
		normal = position - center;
		normal.Normalize();
	}

	inline const VT UnitNormal(const VT& position) const
	{
		VT unitnormal = position - center;
		unitnormal.Normalize();

		return unitnormal;
	}
};




