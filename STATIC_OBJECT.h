#pragma once

#include "SIMULATION_OBJECT.h"
#include "PLANE.h"

class STATIC_OBJECT : public SIMULATION_OBJECT
{
public: // Essential Data
	using SIMULATION_OBJECT::aabb;
	using SIMULATION_OBJECT::velocity;
	using SIMULATION_OBJECT::levelset_object;

public: // Constructors and Destructor
	STATIC_OBJECT(void)
	{}

	STATIC_OBJECT(PLANE* plane, const GRID_STRUCTURE_3D& domain_grid)
	{
		levelset_object = plane;

		// Build AABB
		aabb.Initialize(domain_grid.min, domain_grid.max);
		aabb.Enlarge((T)1*domain_grid.dx);
		if (aabb.Intersection(*plane) == true)
		{
			aabb.Enlarge((T)3*domain_grid.dx);		// Enlarge aabb if plane-domain intersection volume exists
		}

		// Initialization
		fluid_control_options.source_velocity = true;
		velocity = VT((T)0, (T)0, (T)0);
		render = false;
		object_type = STATIC;
	}

	~STATIC_OBJECT(void)
	{}

public: // Member Functions
	void UpdateBoundingBoxes()
	{}

	bool InsideAABB(const VT& position) const
	{
		return aabb.Inside(position);
	}
};

