#pragma once

#include "LEVELSET_3D.h"
#include "ADVECTION_METHOD_3D.h"
//#include "FAST_ITERATIVE_METHOD.h"

class ADVECTION_3D
{
public:
	MULTITHREADING&			multithreading;
	
public:
	bool					use_water_solver;
	
	// Water Levelset
	LEVELSET_3D&			water_levelset;
	FIELD_STRUCTURE_3D<T>&	water_signed_distance_field;
	
	// Boundary Levelset
	LEVELSET_3D&			boundary_levelset;

	// Density
	FIELD_STRUCTURE_3D<T>*	pressure_field;

	GRID_STRUCTURE_3D&		base_grid;
	
	// For MAC Grid
	FIELD_STRUCTURE_3D<T>&  water_velocity_field_mac_x;
	FIELD_STRUCTURE_3D<T>&  water_velocity_field_mac_y;
	FIELD_STRUCTURE_3D<T>&  water_velocity_field_mac_z;

	FIELD_STRUCTURE_3D<T>&	scalar_field_ghost;
	FIELD_STRUCTURE_3D<T>&	velocity_field_mac_ghost_x;
	FIELD_STRUCTURE_3D<T>&  velocity_field_mac_ghost_y;
	FIELD_STRUCTURE_3D<T>&  velocity_field_mac_ghost_z;

	// For boundary condition
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_x_x;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_x_y;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_x_z;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_y_x;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_y_y;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_y_z;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_z_x;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_z_y;
	FIELD_STRUCTURE_3D<T>	velocity_field_mac_ghost_z_z;

public: // For Boundary Condition of Cylindrical Pipe
	FIELD_STRUCTURE_3D<T>	boundary_phi_x;
	FIELD_STRUCTURE_3D<T>	boundary_phi_y;
	FIELD_STRUCTURE_3D<T>	boundary_phi_z;

public: // Options
	bool					air_water_simulation;
	bool					oil_water_simulation;
	
	// For Levelset Advection
	bool					use_csl, use_5th_weno, use_3rd_eno;

	// For Velocity Advection
	bool					use_5th_weno_v, use_3rd_eno_v;
		
	// For Pipe
	bool					is_vertical, is_parallel;
	
	// Options for Boundary Condition
	bool					levelset_advection, velocity_advection;

public: // Condition variable
	T						epsilon, epsilon_v;

public: // Pipe Radius
	T						a;

public: // Constructor and Desructor
	ADVECTION_3D(MULTITHREADING &multithreading_input, LEVELSET_3D& water_levelset_input, FIELD_STRUCTURE_3D<T>& scalar_field_ghost_input, LEVELSET_3D& boundary_levelset_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_x_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_y_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_z_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_ghost_x_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_ghost_y_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_ghost_z_input) 
				 :multithreading(multithreading_input), water_levelset(water_levelset_input), scalar_field_ghost(scalar_field_ghost_input), boundary_levelset(boundary_levelset_input),
				 water_velocity_field_mac_x(velocity_field_mac_x_input), water_velocity_field_mac_y(velocity_field_mac_y_input), water_velocity_field_mac_z(velocity_field_mac_z_input), velocity_field_mac_ghost_x(velocity_field_mac_ghost_x_input), velocity_field_mac_ghost_y(velocity_field_mac_ghost_y_input), velocity_field_mac_ghost_z(velocity_field_mac_ghost_z_input),
				 base_grid(water_levelset_input.grid), 
				 water_signed_distance_field(water_levelset_input.signed_distance_field), pressure_field(0),
				 air_water_simulation(false), oil_water_simulation(false),
				 is_vertical(false), is_parallel(false)
	{
		water_signed_distance_field.is_scalar = true;
		use_5th_weno = false;
		use_3rd_eno = false;
		use_5th_weno_v = false;
		use_3rd_eno_v = false;
	}

	~ADVECTION_3D(void)
	{}

public: // Initialization Function
	void InitializeFromBlock(const SCRIPT_BLOCK& outer_block)
	{
		// script_block is "FLUID_SOLVER_UNIFORM"
		const SCRIPT_BLOCK& advection_block(outer_block.FindBlock("LEVELSET_ADVECTION"));
		use_csl = advection_block.GetBoolean("use_csl", false);
		use_5th_weno = advection_block.GetBoolean("use_5th_weno", false);
		use_3rd_eno = advection_block.GetBoolean("use_3rd_eno", false);
		
		epsilon = advection_block.GetFloat("epsilon", (T)10e-6);

		const SCRIPT_BLOCK& advection_block_v(outer_block.FindBlock("VELOCITY_ADVECTION"));
		use_5th_weno_v = advection_block_v.GetBoolean("use_5th_weno_v", false);
		use_3rd_eno_v = advection_block_v.GetBoolean("use_3rd_eno_v", false);
		epsilon_v = advection_block_v.GetFloat("epsilon", (T)10e-6);

		// Display
		cout << "--------------LEVELSET_ADVECTION--------------" << endl;
		if (use_csl)
		{
			cout << "use csl: " << "true" << endl;
		}
		else
		{
			cout << "use csl: " << "false" << endl;
		}
		
		if (use_5th_weno)
		{
			cout << "use 5th weno: " << "true" << endl;
		}
		else
		{
			cout << "use 5th weno: " << "false" << endl;
		}
		
		if (use_3rd_eno)
		{
			cout << "use 3rd eno: " << "true" << endl;
		}
		else
		{
			cout << "use 3rd eno: " << "false" << endl;
		}
		
		cout << "epsilon: " << epsilon << endl;

		cout << "--------------VELOCITY_ADVECTION--------------" << endl;
		if (use_5th_weno_v)
		{
			cout << "use 5th weno: " << "true" << endl;
		}
		else
		{
			cout << "use 5th weno: " << "false" << endl;
		}
		
		if (use_3rd_eno_v)
		{
			cout << "use 3rd eno: " << "true" << endl;
		}
		else
		{
			cout << "use 3rd eno: " << "false" << endl;
		}
		
		cout << "epsilon: " << epsilon_v << endl;

		// For boundary condition of cylindrical pipe
		boundary_phi_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		boundary_phi_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		boundary_phi_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		
		// For boundary condition
		velocity_field_mac_ghost_x_x.Initialize(&multithreading, velocity_field_mac_ghost_x.grid, 3);
		velocity_field_mac_ghost_x_y.Initialize(&multithreading, velocity_field_mac_ghost_x.grid, 3);
		velocity_field_mac_ghost_x_z.Initialize(&multithreading, velocity_field_mac_ghost_x.grid, 3);
		velocity_field_mac_ghost_y_x.Initialize(&multithreading, velocity_field_mac_ghost_y.grid, 3);
		velocity_field_mac_ghost_y_y.Initialize(&multithreading, velocity_field_mac_ghost_y.grid, 3);
		velocity_field_mac_ghost_y_z.Initialize(&multithreading, velocity_field_mac_ghost_y.grid, 3);
		velocity_field_mac_ghost_z_x.Initialize(&multithreading, velocity_field_mac_ghost_z.grid, 3);
		velocity_field_mac_ghost_z_y.Initialize(&multithreading, velocity_field_mac_ghost_z.grid, 3);
		velocity_field_mac_ghost_z_z.Initialize(&multithreading, velocity_field_mac_ghost_z.grid, 3);

		velocity_field_mac_ghost_x_x.is_x_component = true;
		velocity_field_mac_ghost_x_y.is_x_component = true;
		velocity_field_mac_ghost_x_z.is_x_component = true;
		velocity_field_mac_ghost_y_x.is_y_component = true;
		velocity_field_mac_ghost_y_y.is_y_component = true;
		velocity_field_mac_ghost_y_z.is_y_component = true;
		velocity_field_mac_ghost_z_x.is_z_component = true;
		velocity_field_mac_ghost_z_y.is_z_component = true;
		velocity_field_mac_ghost_z_z.is_z_component = true;
	}

public: // Member Functions
	void Solve_Levelset(const int& thread_id, const T& dt)
	{
		levelset_advection = true;

		SetupBoundaryConditionForVelocity(thread_id);

		if (use_5th_weno)
		{
			if (is_vertical)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInYDirection(thread_id, water_signed_distance_field.array_for_this, true);
				multithreading.Sync(thread_id);
			}
			
			if (is_parallel)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInXDirection(thread_id, water_signed_distance_field.array_for_this, true);
				multithreading.Sync(thread_id);
			}

			// Solve for levelset
			if (air_water_simulation)
			{
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon);
			}
			if (oil_water_simulation)
			{
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_signed_distance_field, scalar_field_ghost, boundary_levelset.signed_distance_field, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon);
			}

			levelset_advection = false;

			
		}
		else if (use_3rd_eno)
		{

		}
	}

	void Solve_Velocity(const int& thread_id, const T& dt)
	{
		velocity_advection = true;
		
		if (is_vertical)
		{
			// Set up the centerline velocity to be constant
			GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
			{
				if ((i == (int)water_velocity_field_mac_x.i_end/2 || i == (int)water_velocity_field_mac_x.i_end/2 + 1) && k == water_velocity_field_mac_x.k_end/2)
				{
					water_velocity_field_mac_x.fixed(i, j, k) = true;
				}
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
			{
				if (i == water_velocity_field_mac_y.i_end/2 && k == water_velocity_field_mac_y.k_end/2)
				{
					water_velocity_field_mac_y.fixed(i, j, k) = true;
				}
				
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
			{
				if (i == water_velocity_field_mac_z.i_end/2 && (k == (int)water_velocity_field_mac_z.k_end/2 || k == (int)water_velocity_field_mac_z.k_end/2 + 1))
				{
					water_velocity_field_mac_z.fixed(i, j, k) = true;
				}
			}
			multithreading.Sync(thread_id);
		}
		if (is_parallel)
		{
			// Set up the centerline velocity to be constant
			GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
			{
				if (j == water_velocity_field_mac_x.j_end/2 && k == water_velocity_field_mac_x.k_end/2)
				{
					water_velocity_field_mac_x.fixed(i, j, k) = true;
				}
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
			{
				if ((j == (int)water_velocity_field_mac_y.j_end/2 || j == (int)water_velocity_field_mac_y.j_end/2 + 1) && k == water_velocity_field_mac_y.k_end/2)
				{
					water_velocity_field_mac_y.fixed(i, j, k) = true;
				}
				
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
			{
				if (j == water_velocity_field_mac_z.j_end/2 && (k == (int)water_velocity_field_mac_z.k_end/2 || k == (int)water_velocity_field_mac_z.k_end/2 + 1))
				{
					water_velocity_field_mac_z.fixed(i, j, k) = true;
				}
			}
			multithreading.Sync(thread_id);
		}

		SetupBoundaryConditionForVelocity(thread_id);

		if (air_water_simulation)
		{
			if (use_5th_weno_v)
			{
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon_v);
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon_v);
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon_v);
			}
			if (use_3rd_eno_v)
			{
				ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
				ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
				ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
			} 
		}

		if (oil_water_simulation)
		{
			if (use_5th_weno_v)
			{
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon_v);
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon_v);
				ADVECTION_METHOD_3D<T>::WENO5th(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading, epsilon_v);
			}
			if (use_3rd_eno_v)
			{
				if (is_vertical)
				{
					ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_z, boundary_phi_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
					ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_z, boundary_phi_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
					ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_x, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_z, boundary_phi_z, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading); 
				}
				if (is_parallel)
				{
					ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_y, velocity_field_mac_ghost_x_z, boundary_phi_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
					ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_y, velocity_field_mac_ghost_y_z, boundary_phi_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
					ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_y, velocity_field_mac_ghost_z_z, boundary_phi_z, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading); 
				}
				/*ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, boundary_phi_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
				ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, boundary_phi_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);
				ADVECTION_METHOD_3D<T>::ENO3rd(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, boundary_phi_z, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_z, dt, multithreading);*/
			} 
		}

		velocity_advection = false;
	}

	void ReinitializationBySussman(const int& thread_id, const T& dt, const FIELD_STRUCTURE_3D<T>& sign_function)
	{
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInYDirection(thread_id, water_signed_distance_field.array_for_this, true);
			}
			if (is_parallel)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInXDirection(thread_id, water_signed_distance_field.array_for_this, true);
			}
		}
		if (air_water_simulation)
		{
			scalar_field_ghost.FillGhostCellsContinuousDerivativesFrom(thread_id, water_signed_distance_field.array_for_this, true);
		}
		ADVECTION_METHOD_3D<T>::WENO5thReinitialization(thread_id, water_signed_distance_field, scalar_field_ghost, dt, multithreading, epsilon, sign_function);
	}

	void SetupBoundaryConditionForVelocity(const int& thread_id)
	{
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_x.x_min + i*boundary_phi_x.dx, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
					boundary_phi_x(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
					//boundary_phi_x(i, j, k) = (T)0.5*(boundary_levelset(i, j, k) + boundary_levelset(i - 1, j, k));	
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_y.x_min + i*boundary_phi_y.dx, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
					boundary_phi_y(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
					//boundary_phi_y(i, j, k) = (T)0.5*(boundary_levelset(i, j, k) + boundary_levelset(i, j - 1, k));	
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_z.x_min + i*boundary_phi_z.dx, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
					boundary_phi_z(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
					//boundary_phi_z(i, j, k) = (T)0.5*(boundary_levelset(i, j, k) + boundary_levelset(i, j, k - 1));	
				}
				multithreading.Sync(thread_id); 
			}
			if (is_parallel)
			{
				GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_x.y_min + j*boundary_phi_x.dy, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
					boundary_phi_x(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
					//boundary_phi_x(i, j, k) = (T)0.5*(boundary_levelset(i, j, k) + boundary_levelset(i - 1, j, k));	
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_y.y_min + j*boundary_phi_y.dy, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
					boundary_phi_y(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
					//boundary_phi_y(i, j, k) = (T)0.5*(boundary_levelset(i, j, k) + boundary_levelset(i, j - 1, k));	
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_z.y_min + j*boundary_phi_z.dy, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
					boundary_phi_z(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
					//boundary_phi_z(i, j, k) = (T)0.5*(boundary_levelset(i, j, k) + boundary_levelset(i, j, k - 1));	
				}
				multithreading.Sync(thread_id); 
			}

			if (levelset_advection)
			{
				if (is_vertical)
				{
					// Periodic Boundary Condition
					velocity_field_mac_ghost_x.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_y.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_z.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);

					GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
					{
						if (boundary_phi_x(i, j, k) > 0)
						{	
							water_velocity_field_mac_x(i, j, k) = (T)0;

							if ((boundary_phi_x(i - 1, j, k) > 0) && (boundary_phi_x(i + 1, j, k) > 0) && (boundary_phi_x(i, j, k - 1) > 0) && (boundary_phi_x(i, j, k + 1) > 0))
							{
								velocity_field_mac_ghost_x(i, j, k) = (T)0;
							}
							if (boundary_phi_x(i - 1, j, k) < 0)
							{
								T theta = abs(boundary_phi_x(i - 1, j, k))/(abs(boundary_phi_x(i - 1, j, k)) + abs(boundary_phi_x(i, j, k)));
								velocity_field_mac_ghost_x(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i - 1, j, k);
							}
							if (boundary_phi_x(i + 1, j, k) < 0)
							{
								T theta = abs(boundary_phi_x(i + 1, j, k))/(abs(boundary_phi_x(i + 1, j, k)) + abs(boundary_phi_x(i, j, k)));
								velocity_field_mac_ghost_x(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i + 1, j, k);
							}
						}
					}
					multithreading.Sync(thread_id);

					GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
					{
						if (boundary_phi_z(i, j, k) > 0)
						{
							water_velocity_field_mac_z(i, j, k) = (T)0;

							if ((boundary_phi_z(i - 1, j, k) > 0) && (boundary_phi_z(i + 1, j, k) > 0) && (boundary_phi_z(i, j, k - 1) > 0) && (boundary_phi_z(i, j, k + 1) > 0))
							{
								velocity_field_mac_ghost_z(i, j, k) = (T)0;
							}
							if (boundary_phi_z(i, j, k - 1) < 0)
							{
								T theta = abs(boundary_phi_z(i, j, k - 1))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
								velocity_field_mac_ghost_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k - 1);
							}
							if (boundary_phi_z(i, j, k + 1) < 0)
							{
								T theta = abs(boundary_phi_z(i, j, k + 1))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
								velocity_field_mac_ghost_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k + 1);
							}
						}	
					}
					multithreading.Sync(thread_id); 
				}
				if (is_parallel)
				{
					// Periodic Boundary Condition
					velocity_field_mac_ghost_x.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_y.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_z.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);

					GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
					{
						if (boundary_phi_y(i, j, k) > 0)
						{	
							water_velocity_field_mac_y(i, j, k) = (T)0;

							if ((boundary_phi_y(i, j - 1, k) > 0) && (boundary_phi_y(i, j + 1, k) > 0) && (boundary_phi_y(i, j, k - 1) > 0) && (boundary_phi_y(i, j, k + 1) > 0))
							{
								velocity_field_mac_ghost_y(i, j, k) = (T)0;
							}
							if (boundary_phi_y(i, j - 1, k) < 0)
							{
								T theta = abs(boundary_phi_y(i, j - 1, k))/(abs(boundary_phi_y(i, j - 1, k)) + abs(boundary_phi_y(i, j, k)));
								velocity_field_mac_ghost_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j - 1, k);
							}
							if (boundary_phi_y(i, j - 1, k) < 0)
							{
								T theta = abs(boundary_phi_y(i, j + 1, k))/(abs(boundary_phi_y(i, j + 1, k)) + abs(boundary_phi_y(i, j, k)));
								velocity_field_mac_ghost_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j + 1, k);
							}
						}
					}
					multithreading.Sync(thread_id);

					GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
					{
						if (boundary_phi_z(i, j, k) > 0)
						{
							water_velocity_field_mac_z(i, j, k) = (T)0;

							if ((boundary_phi_z(i, j - 1, k) > 0) && (boundary_phi_z(i, j + 1, k) > 0) && (boundary_phi_z(i, j, k - 1) > 0) && (boundary_phi_z(i, j, k + 1) > 0))
							{
								velocity_field_mac_ghost_z(i, j, k) = (T)0;
							}
							if (boundary_phi_z(i, j, k - 1) < 0)
							{
								T theta = abs(boundary_phi_z(i, j, k - 1))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
								velocity_field_mac_ghost_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k - 1);
							}
							if (boundary_phi_z(i, j, k + 1) < 0)
							{
								T theta = abs(boundary_phi_z(i, j, k + 1))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
								velocity_field_mac_ghost_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k + 1);
							}
						}	
					}
					multithreading.Sync(thread_id); 
				}
			}
			
			if (velocity_advection)
			{
				if (is_vertical)
				{
					// Periodic Boundary Condition
					velocity_field_mac_ghost_x.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_y.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_z.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);

					velocity_field_mac_ghost_x_x.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_x_z.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_y_x.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_y_z.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_z_x.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);
					velocity_field_mac_ghost_z_z.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);

					GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
					{
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i + 1, j, k) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i + 1, j, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_x(i + 1, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_x(i + t, j, k) = pow(10*t,t); 
							}
						}
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i - 1, j, k) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i - 1, j, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_x(i - 1, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_x(i - t, j, k) = pow(10*t,t);
							}
						}
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i, j, k + 1) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i, j, k + 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k + 1) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_z(i, j, k + t) = pow(10*t,t); 
							}
						}
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i, j, k - 1) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i, j, k - 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k - 1) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_z(i, j, k - t) = pow(10*t,t); 
							}
						}
					}
					multithreading.Sync(thread_id);

					GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
					{
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i + 1, j, k) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i + 1, j, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_x(i + 1, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_x(i + t, j, k) = pow(10*t,t); 
							}
						}
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i - 1, j, k) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i - 1, j, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_x(i - 1, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_x(i - t, j, k) = pow(10*t,t);
							}
						}
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i, j, k + 1) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k + 1) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_z(i, j, k + t) = pow(10*t,t); 
							}
						}
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i, j, k - 1) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k - 1) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_z(i, j, k - t) = pow(10*t,t); 
							}
						}
					}
					multithreading.Sync(thread_id);

					GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
					{
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i + 1, j, k) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i + 1, j, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_x(i + 1, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_x(i + t, j, k) = pow(10*t,t); 
							}
						}
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i - 1, j, k) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i - 1, j, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_x(i - 1, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_x(i - t, j, k) = pow(10*t,t);
							}
						}
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i, j, k + 1) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i, j, k + 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k + 1) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_z(i, j, k + t) = pow(10*t,t); 
							}
						}
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i, j, k - 1) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i, j, k - 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k - 1) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_z(i, j, k - t) = pow(10*t,t); 
							}
						}
					}
					multithreading.Sync(thread_id);
				}
				if (is_parallel)
				{
					// Periodic Boundary Condition
					velocity_field_mac_ghost_x.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_y.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_z.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);

					velocity_field_mac_ghost_x_y.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_x_z.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
					velocity_field_mac_ghost_y_y.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_y_z.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
					velocity_field_mac_ghost_z_y.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);
					velocity_field_mac_ghost_z_z.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);

					GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
					{
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i, j + 1, k) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i, j + 1, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_y(i, j + 1, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_y(i, j + t, k) = pow(10*t,t); 
							}
						}
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i, j - 1, k) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i, j - 1, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_y(i, j - 1, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_y(i, j - t, k) = pow(10*t,t);
							}
						}
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i, j, k + 1) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i, j, k + 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k + 1) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_z(i, j, k + t) = pow(10*t,t); 
							}
						}
						if (boundary_phi_x(i, j, k) < 0 && boundary_phi_x(i, j, k - 1) > 0)
						{
							T theta = abs(boundary_phi_x(i, j, k))/(abs(boundary_phi_x(i, j, k - 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k - 1) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_x_z(i, j, k - t) = pow(10*t,t); 
							}
						}
					}
					multithreading.Sync(thread_id);

					GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
					{
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i, j + 1, k) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i, j + 1, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_y(i, j + 1, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_y(i, j + t, k) = pow(10*t,t); 
							}
						}
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i, j - 1, k) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i, j - 1, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_y(i, j - 1, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_y(i, j - t, k) = pow(10*t,t);
							}
						}
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i, j, k + 1) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k + 1) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_z(i, j, k + t) = pow(10*t,t); 
							}
						}
						if (boundary_phi_z(i, j, k) < 0 && boundary_phi_z(i, j, k - 1) > 0)
						{
							T theta = abs(boundary_phi_z(i, j, k))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k - 1) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_z_z(i, j, k - t) = pow(10*t,t); 
							}
						}
					}
					multithreading.Sync(thread_id);

					GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
					{
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i, j + 1, k) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i, j + 1, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_y(i, j + 1, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_y(i, j + t, k) = pow(10*t,t); 
							}
						}
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i, j - 1, k) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i, j - 1, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_y(i, j - 1, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_y(i, j - t, k) = pow(10*t,t);
							}
						}
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i, j, k + 1) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i, j, k + 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k + 1) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_z(i, j, k + t) = pow(10*t,t); 
							}
						}
						if (boundary_phi_y(i, j, k) < 0 && boundary_phi_y(i, j, k - 1) > 0)
						{
							T theta = abs(boundary_phi_y(i, j, k))/(abs(boundary_phi_y(i, j, k - 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k - 1) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k);

							for (int t = 2; t <= 4 ; t++)
							{
								velocity_field_mac_ghost_y_z(i, j, k - t) = pow(10*t,t); 
							}
						}
					}
					multithreading.Sync(thread_id);
				}
			}
		}
	}
};