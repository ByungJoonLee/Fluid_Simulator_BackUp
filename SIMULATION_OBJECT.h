#pragma once

#include "GRID_STRUCTURE_3D.h"
#include "LEVELSET_OBJECT.h"
#include "BOX.h"
#include "FLUID_CONTROL_OPTIONS.h"
#include "SCRIPT_READER.h"

enum SIMULATION_OBJECT_TYPE
{
	STATIC = 0,
	RIGID,
	FLUID,
	DEFORMABLE, 
	PARTICLE,
	SEQUENCE,
	KEYFRAME,
};

class SIMULATION_OBJECT
{
public: // Object Properties
	SIMULATION_OBJECT_TYPE object_type;
	string name; // name is required to identify simulation object in property controls in plug-ins
	bool discretize;

	VT velocity;
	bool render;
	VT color;

public: // Dynamics solver options
	FLUID_CONTROL_OPTIONS fluid_control_options;

public: // Geometry
	LEVELSET_OBJECT* levelset_object;

public: // Spatial partitioning
	BOX aabb;							// AABB (axis-aligned bounding box) can make spatial partitioning faster
	BOX obb;							// OBB (oriented bounding box) can make spatial partitioning faster
	GRID_STRUCTURE_3D aabb_grid;		// This must include ghost indices
	ARRAY<GRID_STRUCTURE_3D> partial_aabb_grids;

public: // Constructors and Destructor
	SIMULATION_OBJECT(void)
		: levelset_object(0), velocity(VT((T)0, (T)0, (T)0)), render(false), object_type(STATIC), discretize(true)
	{}

	SIMULATION_OBJECT(LEVELSET_OBJECT* levelset_object_input)
	{
		Initialize(levelset_object_input);
	}

	~SIMULATION_OBJECT(void)
	{
		DELETE_POINTER(levelset_object);
	}

public: // Initialization Function
	void Initialize(LEVELSET_OBJECT* levelset_object_input)
	{
		levelset_object = levelset_object_input;

		velocity = VT((T)0, (T)0, (T)0);
		render = false;
		object_type = STATIC;
		discretize = true;
	}

public: // Virtual Functions
	virtual const T		SignedDistance(const VT& position) const
	{
		assert(levelset_object != 0);

		return levelset_object->SignedDistance(position);
	}

	virtual const bool	Inside(const VT& position) const
	{
		assert(levelset_object != 0);

		if(levelset_object->SignedDistance(position) <= 0) return true;
		else return false;
	}

	virtual const VT	ClosestPoint(const VT& position) const
	{
		assert(levelset_object != 0);

		return position - levelset_object->UnitNormal(position)*levelset_object->SignedDistance(position);
	}

	virtual const VT	Velocity(const VT& position) const
	{
		assert(levelset_object != 0);

		return velocity;
	}

	virtual const void	Normal(const VT& position, VT& normal) const
	{
		assert(levelset_object != 0);

		return levelset_object->UnitNormal(position, normal);
	}

	virtual const VT	Normal(const VT& position) const
	{
		assert(levelset_object != 0);

		return levelset_object->UnitNormal(position);
	}

	virtual const void	UnitNormal(const VT& position, VT& normal) const
	{
		assert(levelset_object != 0);

		return levelset_object->UnitNormal(position, normal);
	}

	virtual const VT	UnitNormal(const VT& position) const
	{
		assert(levelset_object != 0);

		return levelset_object->UnitNormal(position);
	}

	virtual void UpdateBoundingBoxes()
	{
		PRINT_AND_EXIT("virtual void SIMULATION_OBJECT::UpdateBoundingBoxes()");
	}

	virtual bool InsideAABB(const VT& position) const
	{
		return true;
	}

	virtual bool InsideOBB(const VT& position) const
	{
		return true;
	}

public: // Member Functions
	bool CompareName(const char* object_name)
	{
		if(name.compare(object_name) == 0) return true;
		else return false;
	}

	void BuildAABBGrid(const int& number_of_threads, const GRID_STRUCTURE_3D& world_grid_ghost)
	{
		const VI index_min(world_grid_ghost.ClampedCell(aabb.min)), index_max(world_grid_ghost.ClampedCell(aabb.max));
		const VI res(index_max - index_min + VI(1, 1, 1));
		const VT aabb_max(aabb.min.x + (T)res.x*world_grid_ghost.dx, aabb.min.y + (T)res.y*world_grid_ghost.dy, aabb.min.z + (T)res.z*world_grid_ghost.dz);

		// Build a partial grid of world_grid_ghost thst is included by AABB
		aabb_grid.Initialize(res, index_min, aabb.min, aabb_max);

		// Split aabb_grid for multithreading
		aabb_grid.SplitInMaxDirection(number_of_threads, partial_aabb_grids);
	}
};



		
