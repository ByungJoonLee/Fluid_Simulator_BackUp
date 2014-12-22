#include "stdafx.h"
#include "SIMULATION_WORLD.h"

void SIMULATION_WORLD::InitializeObjectListFromScript(SCRIPT_READER& script_reader)
{
	GRID_STRUCTURE_3D& base_grid = world_discretization.world_grid;

	SCRIPT_BLOCK script_block = script_reader.FindBlock("SIMULATION_WORLD");

	const bool use_field_boundary = script_block.GetBoolean("use_field_boundary", false);
	const bool top = script_block.GetBoolean("field_boundary_top", false);
	const bool bottom = script_block.GetBoolean("field_boundary_bottom", false);
	const bool left = script_block.GetBoolean("field_boundary_left", false);
	const bool right = script_block.GetBoolean("field_boundary_right", false);
	const bool back = script_block.GetBoolean("field_boundary_back", false);
	const bool front = script_block.GetBoolean("field_boundary_front", false);

	const T wall_width = (T)0;
	const VT left_bottom(base_grid.x_min + base_grid.dx*wall_width, base_grid.y_min + base_grid.dy*wall_width, base_grid.z_min + base_grid.dz*wall_width);
	const VT right_up(base_grid.x_max - base_grid.dx*wall_width, base_grid.y_max - base_grid.dy*wall_width, base_grid.z_max - base_grid.dz*wall_width);

	// Boundary Check
	if (use_field_boundary)
	{
		if (left == true)
		{
			SIMULATION_OBJECT* wall_object = new STATIC_OBJECT(new PLANE(left_bottom, VT(1,0,0)), world_discretization.world_grid);
			object_list.Push(wall_object);
		}
		if (right == true)
		{
			SIMULATION_OBJECT* wall_object = new STATIC_OBJECT(new PLANE(right_up, VT(-1,0,0)), world_discretization.world_grid);
			object_list.Push(wall_object);
		}
		if (bottom == true)
		{
			SIMULATION_OBJECT* wall_object = new STATIC_OBJECT(new PLANE(left_bottom, VT(0,1,0)), world_discretization.world_grid);
			object_list.Push(wall_object);
		}
		if (top == true)
		{
			SIMULATION_OBJECT* wall_object = new STATIC_OBJECT(new PLANE(right_up, VT(0,-1,0)), world_discretization.world_grid);
			object_list.Push(wall_object);
		}
		if (front == true)
		{
			SIMULATION_OBJECT* wall_object = new STATIC_OBJECT(new PLANE(left_bottom, VT(0,0,1)), world_discretization.world_grid);
			object_list.Push(wall_object);
		}
		if (back == true)
		{
			SIMULATION_OBJECT* wall_object = new STATIC_OBJECT(new PLANE(right_up, VT(0,0,-1)), world_discretization.world_grid);
			object_list.Push(wall_object);
		}
	}

	// Object list
	script_block = script_reader.FindBlock("OBJECT_LIST");
	// Water block
	READ_SCRIPT_BLOCKS(script_block, water_block)
	{
		if (!strcmp("WATER", (*water_block)->GetString("shape_type", "NONE")))
		{
			if ((*water_block)->GetBoolean("use_sea_surface", false))
			{
				eulerian_solver.water_height = (*water_block)->GetFloat("height", (T)0);
				PLANE water_surface(VT((T)0, eulerian_solver.water_height, (T)0), VT((T)0, (T)1, (T)0));

				LEVELSET_3D& water_levelset = *eulerian_solver.water_levelset;
				GRID_STRUCTURE_3D& water_grid = water_levelset.grid;

				int max_i, max_j, max_k, i, j, k;
				water_grid.ClampedCell(VT((T)0, eulerian_solver.water_height, (T)0), max_i, max_j, max_k);

				// Assigning the levelset of given surface
				LOOPS_3D(i, j, k, water_grid.i_start, water_grid.j_start, water_grid.k_start, water_grid.i_end, max_j + 3, water_grid.k_end)
				{
					water_levelset(i, j, k) = water_surface.SignedDistance(water_grid.CellCenter(i, j, k));
				}
				water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
			}
			else if ((*water_block)->GetBoolean("use_sphere_bubble", false))
			{
				const VT position = (*water_block)->GetVector3("position", VT(0, 0, 0));
				const T size = (*water_block)->GetFloat("size", T(0));
				
				SPHERE sphere_surface(position, size);

				LEVELSET_3D& water_levelset = *eulerian_solver.water_levelset;
				GRID_STRUCTURE_3D& water_grid = water_levelset.grid;

				int i, j, k;
									
				// Assigning the levelset of given surface
				LOOPS_3D(i, j, k, water_grid.i_start, water_grid.j_start, water_grid.k_start, water_grid.i_end, water_grid.j_end, water_grid.k_end)
				{
					water_levelset(i, j, k) = sphere_surface.SignedDistance(water_grid.CellCenter(i, j, k));
				}
                
                water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
            }
			else if ((*water_block)->GetBoolean("user_defined", false))
			{
                LEVELSET_3D& water_levelset = *eulerian_solver.water_levelset;
                GRID_STRUCTURE_3D& water_grid = water_levelset.grid;

                int i, j, k;

                LOOPS_3D(i, j, k, water_grid.i_start, water_grid.j_start, water_grid.k_start, water_grid.i_end, water_grid.j_end, water_grid.k_end)
				{
				    T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy, z_coor = water_grid.z_min + k*water_grid.dz;
                    water_levelset(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - ((T)0.01*cos((T)2*x_coor) + (T)1/3);
				}
                water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
			}
			if (air_water_simulation)
			{
                LEVELSET_3D& water_levelset = *eulerian_solver.water_levelset;
                GRID_STRUCTURE_3D& water_grid = water_levelset.grid;

                int i, j, k;

                LOOPS_3D(i, j, k, water_grid.i_start, water_grid.j_start, water_grid.k_start, water_grid.i_end, water_grid.j_end, water_grid.k_end)
				{
				    T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy, z_coor = water_grid.z_min + k*water_grid.dz;
                    water_levelset(i, j, k) = sqrt(POW2(x_coor) + POW2(y_coor) + POW2(z_coor)) - (T)1/3;
				}
                water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
            }
			if (oil_water_simulation)
			{
                LEVELSET_3D& water_levelset = *eulerian_solver.water_levelset;
                LEVELSET_3D& boundary_levelset = *eulerian_solver.boundary_levelset;
				GRID_STRUCTURE_3D& water_grid = water_levelset.grid;
				GRID_STRUCTURE_3D& boundary_grid = boundary_levelset.grid;

				// Interface
				/*GRID_ITERATION_3D(water_grid)
				{
					T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy, z_coor = water_grid.z_min + k*water_grid.dz;
					
					if (is_vertical)
					{
						T magnitude = sqrt(POW2(x_coor) + POW2(z_coor));
						
						if (magnitude < (A_0*cos(alpha*y_coor) + (T)1))
						{
							water_levelset(i, j, k) = (T)-100;
						}
						else if (magnitude > (A_0*cos(alpha*y_coor) + (T)1))
						{
							water_levelset(i, j, k) = (T)100;
						}
						else
						{
							water_levelset(i, j, k) = (T)0;
						}
					}
					if (is_parallel)
					{
						T magnitude = sqrt(POW2(y_coor) + POW2(z_coor));
						
						if (magnitude < (A_0*cos(alpha*x_coor) + (T)1))
						{
							water_levelset(i, j, k) = (T)-100;
						}
						else if (magnitude > (A_0*cos(alpha*x_coor) + (T)1))
						{
							water_levelset(i, j, k) = (T)100;
						}
						else
						{
							water_levelset(i, j, k) = (T)0;
						}
					}
				}*/

				/*GRID_ITERATION_3D(water_grid)
				{
					T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy, z_coor = water_grid.z_min + k*water_grid.dz;
					
					if (is_vertical)
					{
						T magnitude = sqrt(POW2(x_coor) + POW2(z_coor));
						
						if (magnitude < (A_0*cos(alpha*y_coor) + 1))
						{
							water_levelset(i, j, k) = (T)-100;
						}
						else if (magnitude > (A_0*cos(alpha*y_coor) + 1))
						{
							water_levelset(i, j, k) = (T)100;
						}
						else
						{
							water_levelset(i, j, k) = (T)0;
						}
					}
					if (is_parallel)
					{
						T magnitude = sqrt(POW2(y_coor) + POW2(z_coor));
						
						if (magnitude < (A_0*cos(alpha*x_coor) + 1))
						{
							water_levelset(i, j, k) = (T)-100;
						}
						else if (magnitude > (A_0*cos(alpha*x_coor) + 1))
						{
							water_levelset(i, j, k) = (T)100;
						}
						else
						{
							water_levelset(i, j, k) = (T)0;
						}
					}
				}

				int number_for_interface(101), number_for_axis(101);
				
				ARRAY<VT> point_for_interface;
				point_for_interface.Initialize(number_for_interface*number_for_axis);

				T step_for_interface = 2*PI/(number_for_interface - 1);
				
				T step_for_axis;
				
				if (is_vertical)
				{
					step_for_axis = (water_grid.y_max - water_grid.y_min)/(number_for_axis - 1);
				}
				if (is_parallel)
				{
					step_for_axis = (water_grid.x_max - water_grid.x_min)/(number_for_axis - 1);
				}
				
				int index_for_interface(0), index_for_axis(0);
				
				if (is_vertical)
				{
					for (index_for_axis = 0; index_for_axis < number_for_axis; index_for_axis++)
					{
						for (index_for_interface = 0; index_for_interface < number_for_interface; index_for_interface++)
						{
							T radius_for_pipe = A_0*cos(alpha*step_for_axis*index_for_axis) + 1;
							point_for_interface[index_for_interface + index_for_axis*number_for_axis] = VT(radius_for_pipe*cos(step_for_interface*index_for_interface), step_for_axis*index_for_axis, radius_for_pipe*sin(step_for_interface*index_for_interface));
						}
					}
				}
				
				if (is_parallel)
				{
					for (index_for_axis = 0; index_for_axis < number_for_axis; index_for_axis++)
					{
						for (index_for_interface = 0; index_for_interface < number_for_interface; index_for_interface++)
						{
							T radius_for_pipe = A_0*cos(alpha*step_for_axis*index_for_axis) + 1;
							point_for_interface[index_for_interface + index_for_axis*number_for_axis] = VT(step_for_axis*index_for_axis, radius_for_pipe*cos(step_for_interface*index_for_interface), radius_for_pipe*sin(step_for_interface*index_for_interface));
						}
					}
				}

				GRID_ITERATION_3D(water_grid)
				{
					T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy, z_coor = water_grid.z_min + k*water_grid.dz;
					
					T min_value = sqrt(POW2(x_coor - point_for_interface[0].x) + POW2(y_coor - point_for_interface[0].y) + POW2(z_coor - point_for_interface[0].z));
				
					for (int l = 0; l < point_for_interface.num_elements; l++)
					{
						min_value = min(min_value, sqrt(POW2(x_coor - point_for_interface[l].x) + POW2(y_coor - point_for_interface[l].y) + POW2(z_coor - point_for_interface[l].z)));
					}

					if (water_levelset(i, j, k) > (T)0)
					{
						water_levelset(i, j, k) = min_value;
					}
					else if (water_levelset(i, j, k) < (T)0)
					{
						water_levelset(i, j, k) = -min_value;
					}
					else
					{
						water_levelset(i, j, k) = 0;
					}
				}*/
				
				int i, j, k;

                LOOPS_3D(i, j, k, water_grid.i_start, water_grid.j_start, water_grid.k_start, water_grid.i_end, water_grid.j_end, water_grid.k_end)
				{
				    T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy, z_coor = water_grid.z_min + k*water_grid.dz;
					if (is_vertical)
					{
						water_levelset(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - (A_0*cos(alpha*y_coor) + (T)1);
					}
					if (is_parallel)
					{
						water_levelset(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - (A_0*cos(alpha*x_coor) + (T)1);
					}
				}

                //water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
				water_levelset.FillGhostCellsContinuousDerivatesFromThreaded(&(boundary_levelset.phi), false);

				LOOPS_3D(i, j, k, boundary_grid.i_start, boundary_grid.j_start, boundary_grid.k_start, boundary_grid.i_end, boundary_grid.j_end, boundary_grid.k_end)
				{
					T x_coor = boundary_grid.x_min + i*boundary_grid.dx, y_coor = boundary_grid.y_min + j*boundary_grid.dy, z_coor = boundary_grid.z_min + k*boundary_grid.dz;
					//boundary_levelset(i, j, k) = -1000;
					if (is_vertical)
					{
						boundary_levelset(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - eulerian_solver.a;
					}
					if (is_parallel)
					{
						boundary_levelset(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - eulerian_solver.a;
					}
				}
				//boundary_levelset.FillGhostCellsFromThreaded(&(boundary_levelset.phi), false);
				boundary_levelset.FillGhostCellsContinuousDerivatesFromThreaded(&(boundary_levelset.phi), false);
			}
		}
	}

	READ_SCRIPT_BLOCKS(script_block, itr_block)
	{
		SIMULATION_OBJECT* simulation_object = 0;

		if (!strcmp("RIGID", (*itr_block)->GetString("object_type", "NONE")))
		{
			RIGID_OBJECT* rigid_object = 0;

			if (!strcmp("OBJ", (*itr_block)->GetString("shape_type", "NONE")))
			{
				rigid_object = new RIGID_OBJECT;
				rigid_object->InitializeFromScriptBlock(multithreading, *(*itr_block), world_discretization.world_grid, 3);
			}

			assert(rigid_object);

			simulation_object = rigid_object;
			simulation_object->object_type = RIGID;

			rigid_object->position = (*itr_block)->GetVector3("position");
			rigid_object->rotation = QUATERNION((*itr_block)->GetFloat("rotation_angle"), (*itr_block)->GetVector3("rotation_direction"));
			rigid_object->linear_momentum = (*itr_block)->GetVector3("linear_momentum");
			rigid_object->angular_momentum = (*itr_block)->GetVector3("angular_momentum");
			rigid_object->force = (*itr_block)->GetVector3("force");
			rigid_object->torque = (*itr_block)->GetVector3("torque");
			rigid_object->static_object = (*itr_block)->GetBoolean("static_object");
			rigid_object->no_gravity = (*itr_block)->GetBoolean("no_gravity");	
		}
		// Need to be updated later
		else
		{
			if (!strcmp("SPHERE", (*itr_block)->GetString("shape_type", "NONE")))
			{
				const VT position = (*itr_block)->GetVector3("position");
				const T size = (*itr_block)->GetFloat("size");
				const T bb_half_width = size + (T)3*world_discretization.dx;

				simulation_object = new SIMULATION_OBJECT(new SPHERE(position, size));
				simulation_object->aabb.Initialize(position - VT(bb_half_width, bb_half_width, bb_half_width), position + VT(bb_half_width, bb_half_width, bb_half_width));
				simulation_object->obb.Initialize(position - VT(bb_half_width, bb_half_width, bb_half_width), position + VT(bb_half_width, bb_half_width, bb_half_width));
			}
			else if (!strcmp("CYLINDER", (*itr_block)->GetString("shape_type", "NONE")))
			{
				CYLINDER* cylinder = new CYLINDER();

				const VT start_position = (*itr_block)->GetVector3("start_position");
				const VT end_position = (*itr_block)->GetVector3("end_position");
				const T radius = (*itr_block)->GetFloat("radius");

				const T bb_half_width = radius + (T)3*world_discretization.dx;

				VT p0, p1;
				if ((start_position - VT()).Magnitude() < ((end_position - VT()).Magnitude()))
				{
					p0 = start_position;
					p1 = end_position;
				}
				else
				{
					p0 = end_position;
					p1 = start_position;
				}
				cylinder->Initialize(start_position, end_position, radius, (*itr_block)->GetInteger("slices"), multithreading, world_discretization.world_grid, 3);

				simulation_object = new SIMULATION_OBJECT(cylinder);
				simulation_object->object_type = STATIC;

				simulation_object->aabb.Initialize(p0 - VT(bb_half_width, bb_half_width, bb_half_width), p1 + VT(bb_half_width, bb_half_width, bb_half_width));
				simulation_object->obb.Initialize(p0 - VT(bb_half_width, bb_half_width, bb_half_width), p1 + VT(bb_half_width, bb_half_width, bb_half_width));
			}
			else if (!strcmp("WALL", (*itr_block)->GetString("shape_type", "NONE")))
			{
				// simulation_object = wall_object;
			}
			else if (!strcmp("OBJ", (*itr_block)->GetString("shape_type", "NONE")))
			{
				TRIANGULAR_SURFACE triangular_surface;
				triangular_surface.ReadOBJ((*itr_block)->GetString("file", 0));
				triangular_surface.Scale((*itr_block)->GetFloat("obj_scale", (T)1));

				GRID_STRUCTURE_3D& grid = world_discretization.world_grid;

				BOX domain(triangular_surface.bounding_box);
				domain.Enlarge(grid.dx*(T)5);

				VT center = (domain.min + domain.max)/(T)2;
				VT position = (*itr_block)->GetVector3("position", center);

				GRID_STRUCTURE_3D levelset_grid;
				levelset_grid.Initialize(domain, grid.dx*(*itr_block)->GetFloat("dx_scale", (T)1));

				// Generate levelset uniform
				LEVELSET_3D* levelset_3d = new LEVELSET_3D;
				levelset_3d->Initialize(multithreading, levelset_grid, 3);

				T old_value = levelset_grid.dx*(T)-3;
				T new_value = levelset_grid.dx*(T)3;

				levelset_3d->AssignAllValuesLevelsetThreaded(old_value);
				levelset_3d->AssignTriangularSurfaceLevelSet(triangular_surface);

				SCAN_LINE_ALGORITHM<T> scan_line_algorithm;
				scan_line_algorithm.Initialize(multithreading, levelset_3d->grid.ijk_res);

				// Register first seed cells
				for (int t = 0; t < multithreading->num_threads; t++)
				{
					GRID_STRUCTURE_3D& bc_grid = levelset_3d->partial_grids[t];
					int i, j, k;
					LOOPS_3D(i, j, k, bc_grid.i_start, bc_grid.j_start, bc_grid.k_start, bc_grid.i_end, bc_grid.j_end, bc_grid.k_end)
					{
						scan_line_algorithm.stack_from_thread[t].Push(VI(i, j, k));
					}
				}

				scan_line_algorithm.ScanFieldNoBlockThreaded(&(levelset_3d->signed_distance_field), old_value, new_value);

				levelset_3d->FullFastSweepingMethodThreaded();
				levelset_3d->UpdateThreaded();

				domain.min = domain.min - center + position;
				domain.max = domain.max - center + position;

				levelset_3d->Translate(-center + position);

				simulation_object = new SIMULATION_OBJECT(levelset_3d);

				simulation_object->aabb.Initialize(domain.min, domain.max);
				simulation_object->obb.Initialize(domain.min, domain.max);
			}

			if (simulation_object != 0)
			{
				if (!strcmp("FLUID", (*itr_block)->GetString("object_type", "NONE")))
				{
					simulation_object->object_type = FLUID;
				}
				else
				{
					simulation_object->object_type = STATIC;
				}
			}
		}

		if (simulation_object != 0)
		{
			// After updating the particle simulation system, you can add the particle related one
			simulation_object->velocity = (*itr_block)->GetVector3("velocity");
			simulation_object->render = (*itr_block)->GetBoolean("render");
			simulation_object->color = (*itr_block)->GetVector3("color");
			
			simulation_object->discretize = (*itr_block)->GetBoolean("discretize", true);

			// Block name is used as simulation object's name
			simulation_object->name = (*itr_block)->block_name;

			// Initialize control option
			SCRIPT_BLOCK option_block;

			option_block = (*itr_block)->FindBlock("FLUID_CONTROL_OPTIONS");
			simulation_object->fluid_control_options.InitializeFromScript(option_block);

			{
				object_list.Push(simulation_object);
			}
		}
	}
}

