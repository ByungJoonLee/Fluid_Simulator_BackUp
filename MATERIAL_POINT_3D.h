#pragma once

#include "COMMON_DEFINITIONS.h"

class MATERIAL_POINT_3D
{
public: // Essential Data
	VT position;
	VT velocity;
	T  T00, T01, T02, T11, T12, T22;

	MATERIAL_POINT_3D* next;

public: // Constructor and Destructor
	MATERIAL_POINT_3D(void)
	{
		next = 0;
		T00 = T01 = T02 = T11 = T12 = T22 = (T)0;
	}

	~MATERIAL_POINT_3D(void)
	{}
};

