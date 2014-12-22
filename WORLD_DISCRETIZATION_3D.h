#pragma once

#include "COMMON_DEFINITIONS.h"
#include "SIMULATION_OBJECT.h"
#include "ARRAY_3D.h"
#include "DYNAMIC_ARRAY.h"
#include "MULTITHREADING.h"
#include "SCRIPT_READER.h"

class WORLD_DISCRETIZATION_3D
{
public: // Essential Data 
	GRID_STRUCTURE_3D						world_grid;
	GRID_STRUCTURE_3D						world_grid_ghost;

	DYNAMIC_ARRAY<SIMULATION_OBJECT*>		object_list;

	FIELD_STRUCTURE_3D<SIMULATION_OBJECT*>	object_field;

	MULTITHREADING*							multithreading;

    // Simulation Options
    bool                                    air_water_simulation;
    bool                                    oil_water_simulation;

	int										ghost_width;

	T										object_padding_width;

	bool									use_grid_uniform;

	// Option For Air Water Simulation
	bool									large_bubble, small_bubble;

	// Option For Oil Water Simulation
	bool									is_vertical, is_parallel;

public: // Speedup constants
	T dx, dy, dz;
	T dx_over_two;
	T dy_over_two;
	T dz_over_two;

public: // Constructors and Destructor
	WORLD_DISCRETIZATION_3D(void)
		: object_list(0), multithreading(0), use_grid_uniform(true), large_bubble(false), small_bubble(false), air_water_simulation(false), oil_water_simulation(false), is_vertical(false), is_parallel(false)
	{}

	~WORLD_DISCRETIZATION_3D(void)
	{
		for (int i = 0; i < object_list.num_of_elements; ++i)
		{
			DELETE_POINTER(object_list.values[i]);
		}
	}

public: // Initialization Functions
	// Need to update this after you update the script block
	void Initialize(MULTITHREADING* multithreading_input, const SCRIPT_BLOCK& script_block)
	{
		for (int i = 0; i < object_list.num_of_elements; ++i)
		{
			DELETE_POINTER(object_list.values[i]);
		}
		
		multithreading = multithreading_input;

		// Display
		cout << "--------------DISCRETIZATION VARIABLES--------------" << endl;
		
		ghost_width = script_block.GetInteger("ghost_width", 1);
			
		if (air_water_simulation)
		{
            large_bubble = script_block.GetBoolean("large_bubble", (bool)false);
			small_bubble = script_block.GetBoolean("small_bubble", (bool)false);
			
			if (large_bubble)
		    {
		    	world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_3D_LARGE"));
		    	world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
		    	use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
		    }
		    if (small_bubble)
		    {
		    	world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_3D_SMALL"));
		    	world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
		    	use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
		    }
		}
		if (oil_water_simulation)
		{
			is_vertical = script_block.FindBlock("PIPE_OPTIONS").GetBoolean("is_vertical", (bool)false);
			is_parallel = script_block.FindBlock("PIPE_OPTIONS").GetBoolean("is_parallel", (bool)false);
			
			if (is_vertical)
			{
				world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_3D_VERTICAL"));
				cout << "Vertical Pipe Simulation is activated!" << endl;
				world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
				use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
			}
			if (is_parallel)
			{
				world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_3D_PARALLEL"));
				cout << "Parallel Pipe Simulation is activated!" << endl;
				world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
				use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
			}
		}

		object_field.Initialize(multithreading, world_grid, ghost_width);

		object_padding_width = (T)script_block.GetInteger("object_padding_width", 2)*world_grid.dx;
		use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);

		InitializeSpeedupConstants();
	}

	void InitializeSpeedupConstants()
	{
		dx = world_grid.dx;
		dy = world_grid.dy;
		dz = world_grid.dz;

		dx_over_two = dx*(T)0.5;
		dy_over_two = dy*(T)0.5;
		dz_over_two = dz*(T)0.5;
	}

	void Update(const int& thread_id)
	{
		if(use_grid_uniform)
		{
			HEAD_THREAD_WORK(memset(object_field.values, 0, sizeof(SIMULATION_OBJECT*) * object_field.ij_res_g););
						
			UpdateAABBGrids(thread_id, world_grid_ghost);

			// Registering object to cells
			for(int x = 0; x < object_list.num_of_elements; ++x)
			{
				SIMULATION_OBJECT* object = object_list.values[x];

				BEGIN_GRID_ITERATION_3D(object->partial_aabb_grids[thread_id])
				{
					const VT cell_center = world_grid_ghost.CellCenter(i, j, k);
					const T obj_phi = object->SignedDistance(cell_center);

					if (obj_phi <= object_padding_width)
					{
						SIMULATION_OBJECT* old_object = object_field(i, j, k);

						if (old_object)
						{
							T old_phi = old_object->SignedDistance(cell_center);

							if (ABS(obj_phi) < ABS(old_phi))
							{
								object_field(i, j, k) = object;
							}
						}
						else
						{
							object_field(i, j, k) = object;
						}
					}
				}
				END_GRID_ITERATION_3D;
			}
		}
	}

	void UpdateAABBGrids(const int& thread_id, const GRID_STRUCTURE_3D& grid)
	{
		PREPARE_FOR_1D_ITERATION(object_list.num_of_elements);
		
		// Build AABB grids of all objects
		BEGIN_PARTICLE_ITERATION
		{
			assert(object_list->values[i] != 0);

			object_list.values[i]->BuildAABBGrid(multithreading->num_threads, grid);
		}
		END_PARTICLE_ITERATION
	}

	SIMULATION_OBJECT* GetObjectList(const VT& position)
	{
		return object_field(world_grid_ghost.ClampedCell(position));
	}
};
			



			



				

			






