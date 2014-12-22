#pragma once

#include "LEVELSET_3D.h"
#include "POISSON_SOLVER_3D.h"
#include "SCAN_LINE_ALGORITHM.h"

enum PROJECTION_TYPE				{FREE_SURFACE_WATER, AIR, MULTIPHASE};

class PROJECTION_3D
{
public: // References of the variables which are defined in EULERIAN_FLUID_SOLVER_3D class
	MULTITHREADING*					multithreading;

	bool							use_water_solver;
	
	// For MAC grid
	FIELD_STRUCTURE_3D<T>&			water_velocity_field_mac_x;
	FIELD_STRUCTURE_3D<T>&			water_velocity_field_mac_y;
	FIELD_STRUCTURE_3D<T>&			water_velocity_field_mac_z;
	
	// For nth step - When you update the given value, nth step must not be changed
	FIELD_STRUCTURE_3D<T>&			vector_field_mac_ghost_x;
	FIELD_STRUCTURE_3D<T>&			vector_field_mac_ghost_y;
	FIELD_STRUCTURE_3D<T>&			vector_field_mac_ghost_z;

	LEVELSET_3D&					boundary_levelset;
	LEVELSET_3D&					water_levelset;

	FIELD_STRUCTURE_3D<T>			scalar_field_ghost;

public: // Variables defined and used in this class. EULERIAN_FLUID_SOLVER_3D doesn't need to have these
	FIELD_STRUCTURE_3D<T>			pressure_field;
	FIELD_STRUCTURE_3D<T>			projection_density_field;
	FIELD_STRUCTURE_3D<int>			boundary_condition_field;
	FIELD_STRUCTURE_3D<T>			divergence_field;

public: // Jump Condition Field
	FIELD_STRUCTURE_3D<T>			jc_on_solution;
	FIELD_STRUCTURE_3D<T>			jc_on_derivative;

public: // Face value for levelset
	FIELD_STRUCTURE_3D<T>			levelset_x_c;
	FIELD_STRUCTURE_3D<T>			levelset_y_c;
	FIELD_STRUCTURE_3D<T>			levelset_z_c;

public: // Face value for density
	FIELD_STRUCTURE_3D<T>			density_half_x;
	FIELD_STRUCTURE_3D<T>			density_half_y;
	FIELD_STRUCTURE_3D<T>			density_half_z;

public: // Control Options
	enum POISSON_SOLVER_TYPE        poisson_solver_type;
	enum SMOOTHER_TYPE				smoother_type;			// for non-coarsest level of multigrid
	enum POISSON_SOLVER_TYPE		coarsest_level_poisson_solver_type;
	enum PROJECTION_TYPE			projection_type;

public: // Subsolver
	POISSON_SOLVER_3D				poisson_solver;
	
	SCAN_LINE_ALGORITHM<int>		scan_line_algorithm;	// For boundary condition field

public: // Convenient variables and references
	GRID_STRUCTURE_3D&				base_grid;
	const T							dx, dy, dz;
	const T							dxdx;
	const T							one_over_dx, one_over_dy, one_over_dz;
	const T							one_over_2dx, one_over_2dy, one_over_2dz;
	const T							half_dx;
	const T							inv_2dx;
	const T							inv_2dy;
	const T							inv_2dz;

	ARRAY_3D<T>&					pressure_array;
	ARRAY_3D<T>&					divergence_array;
	ARRAY_3D<T>&					boundary_levelset_array;
	ARRAY_3D<int>&					boundary_condition_array;
	ARRAY_3D<T>&					pressure_density_array;

	T								tolerance;
	int								max_iteration;
	
	T								max_velocity_magnitude;
	T								max_velocity_x;
	T								max_velocity_y;
	T								max_velocity_z;

	T								water_density, air_density, oil_density;
	T								surface_tension;

	bool							use_previous_pressure;
	bool							use_2nd_order;
	bool							remove_neumann_pocket;
	bool							place_dirichlet_bc_at_face;
	bool							use_variable_density;

	// Option For Bubble
	bool							air_bubble_rising, water_drop;

	// Option For Cylinder
	bool							regular_boundary, cylindrical_boundary;

	// Options For Simulation
	bool							air_water_simulation;
	bool							oil_water_simulation;
	bool							dimensionless_form;
	bool							CSF_model;

	// Option For Pipe
	bool							is_vertical, is_parallel;
	
	// Option For Viscosity treatment
	bool							use_delta_function_formulation;

	// Option For Boundary Condition
	bool							Dirichlet_Boundary_Condition, Neumann_Boundary_Condition;

	int&							frame;

	// Helpful Field
	FIELD_STRUCTURE_3D<T>			one_over_density;
		 
public: // Dimensionless Variable
	T								We;
	T								R1;

	T								m, eta;

public: // Pipe Radius
	T								a;

public: // For Cylinder Pipe
	FIELD_STRUCTURE_3D<T>			boundary_phi_x;
	FIELD_STRUCTURE_3D<T>			boundary_phi_y;
	FIELD_STRUCTURE_3D<T>			boundary_phi_z;

public: // Constructors and Destructor
	PROJECTION_3D(MULTITHREADING* multithreading_input, LEVELSET_3D& water_levelset_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_x_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_y_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_z_input, 
				  FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_ghost_x_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_ghost_y_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_ghost_z_input, LEVELSET_3D& boundary_levelset_input, int& frame_input, 
				  bool use_water_solver_input)
				  : // Main variables
					multithreading(multithreading_input),
				    water_velocity_field_mac_x(water_velocity_field_mac_x_input), 
					water_velocity_field_mac_y(water_velocity_field_mac_y_input),
					water_velocity_field_mac_z(water_velocity_field_mac_z_input),
					vector_field_mac_ghost_x(water_velocity_field_mac_ghost_x_input),
					vector_field_mac_ghost_y(water_velocity_field_mac_ghost_y_input),
					vector_field_mac_ghost_z(water_velocity_field_mac_ghost_z_input),
					water_levelset(water_levelset_input),
					boundary_levelset(boundary_levelset_input),
					// Speedup variables
					base_grid(water_levelset_input.grid),
					dx(base_grid.dx), dy(base_grid.dy), dz(base_grid.dz),
					dxdx(base_grid.dx*base_grid.dx), 
					one_over_dx((T)1/dx), one_over_dy((T)1/dy), one_over_dz((T)1/dz),
					one_over_2dx((T)0.5*one_over_dx), one_over_2dy((T)0.5*one_over_dy), one_over_2dz((T)0.5*one_over_dz),
					inv_2dx((T)0.5/dx), inv_2dy((T)0.5/dy), inv_2dz((T)0.5/dz),
					half_dx((T)0.5*dx),
					pressure_array(pressure_field.array_for_this), divergence_array(divergence_field.array_for_this), boundary_levelset_array(boundary_levelset.arr),
					boundary_condition_array(boundary_condition_field.array_for_this), pressure_density_array(projection_density_field.array_for_this),
					max_velocity_magnitude((T)0), use_previous_pressure(false),
					use_water_solver(use_water_solver_input), 
					frame(frame_input), air_bubble_rising(false), water_drop(false),
					max_velocity_x((T)0), max_velocity_y((T)0), max_velocity_z((T)0),
					air_density((T)0), water_density((T)0), oil_density((T)0),
					air_water_simulation(false), oil_water_simulation(false), dimensionless_form(false), Dirichlet_Boundary_Condition(false), Neumann_Boundary_Condition(false), CSF_model(false), regular_boundary(false), cylindrical_boundary(false), is_vertical(false), is_parallel(false)
	{}

	~PROJECTION_3D(void)
	{
	}

public: // Initialization Functions
	void InitializeFromBlock(const SCRIPT_BLOCK& projection_block)
	{
		tolerance = projection_block.GetFloat("tolerance", (T)1e-4);
		max_iteration = projection_block.GetInteger("max_iteration", 30);
		use_previous_pressure = projection_block.GetBoolean("use_previous_pressure", false);
		use_2nd_order = projection_block.GetBoolean("use_2nd_order", false);
		use_variable_density = projection_block.GetBoolean("use_variable_density", false);
		remove_neumann_pocket = projection_block.GetBoolean("remove_neumann_pocket", false);
		place_dirichlet_bc_at_face = projection_block.GetBoolean("place_dirichlet_bc_at_face", true);
		Dirichlet_Boundary_Condition = projection_block.GetBoolean("Dirichlet_Boundary_Condition", (bool)false);
		Neumann_Boundary_Condition = projection_block.GetBoolean("Neumann_Boundary_Condition", (bool)false);
		
		if (air_water_simulation)
		{
			air_bubble_rising = projection_block.GetBoolean("air_bubble_rising", false);
			water_drop = projection_block.GetBoolean("water_drop", false);
		}
		
		if (oil_water_simulation)
		{
			regular_boundary = projection_block.GetBoolean("regular_boundary", (bool)false);
			cylindrical_boundary = projection_block.GetBoolean("cylindrical_boundary", (bool)false);
		}

		// Poisson solver type from scipt
		const char* poisson_solver_type_input = projection_block.GetString("poisson_solver_type", "NULL");
		if (!strcmp(poisson_solver_type_input, "CG"))
		{
			poisson_solver_type = CG;
		}
		else if (!strcmp(poisson_solver_type_input, "PCG"))
		{
			poisson_solver_type = PCG;
		}
		else if (!strcmp(poisson_solver_type_input, "HYBRID"))
		{
			poisson_solver_type = HYBRID;
		}
		else if (!strcmp(poisson_solver_type_input, "GS"))
		{
			poisson_solver_type = GS;
		}
		else
		{
			poisson_solver_type = CG;
		}

		const char* smoother_name = projection_block.GetString("smoother_type", "NULL");
		if (!strcmp(smoother_name, "GS_SMOOTHER"))
		{
			smoother_type = GS_SMOOTHER;
		}
		else if (!strcmp(smoother_name, "NB_GS_SMOOTHER"))
		{
			smoother_type = NB_GS_SMOOTHER;
		}
		else if (!strcmp(smoother_name, "WEIGHTED_JACOBI_SMOOTHER"))
		{
			smoother_type = WEIGHTED_JACOBI_SMOOTHER;
		}
		else                      // Default
		{
			smoother_type = WEIGHTED_JACOBI_SMOOTHER;
		}

		const char* coarsest_level_poisson_solver_type_input = projection_block.GetString("coarsest_levelset_poisson_solver_type", "NULL");
		if (!strcmp(coarsest_level_poisson_solver_type_input, "CG"))
		{
			coarsest_level_poisson_solver_type = CG;
		}
		else if (!strcmp(coarsest_level_poisson_solver_type_input, "PCG"))
		{
			coarsest_level_poisson_solver_type = PCG;
		}
		else if (!strcmp(coarsest_level_poisson_solver_type_input, "HYBRID"))
		{
			coarsest_level_poisson_solver_type = HYBRID;
		}
		else if (!strcmp(coarsest_level_poisson_solver_type_input, "GS"))
		{
			coarsest_level_poisson_solver_type = GS;
		}
		else							// Default
		{
			coarsest_level_poisson_solver_type = CG;
		}

		cout << "--------------PROJECTION--------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max iteration: " << max_iteration << endl;
				
		switch (poisson_solver_type)
		{
		case NO_SOLVER:
			cout << "poisson solver type: " << "NO SOLVER" << endl;
			break;
		
		case MULTIGRID:
			cout << "poisson solver type: " << "MULTIGRID" << endl;
			break;

		case CG:
			cout << "poisson solver type: " << "CG" << endl;
			break;
		
		case PCG:
			cout << "poisson solver type: " << "PCG" << endl;
			break;
		
		case HYBRID:
			cout << "poisson solver type: " << "HYBRID" << endl;
			break;
		
		case GS:
			cout << "poisson solver type: " << "GS" << endl;
			break;
		default:
			break;
		}
		
		if (air_water_simulation)
		{
			if (air_bubble_rising)
			{
				cout << "Air Bubble Rising is activated!" << endl;
			}
			if (water_drop)
			{
				cout << "Water Drop is activated!" << endl;
			}
		}
		
		if (oil_water_simulation)
		{
			cout << "Oil Water Simulation is activated!" << endl;
		}

		if (use_variable_density)
		{
            cout << "use variable density: " << "true" << endl;
		}
		else
		{
            cout << "use variable density: " << "false" << endl;
		}

		if (place_dirichlet_bc_at_face)
		{
            cout << "place dirichlet bc at face: " << "true" << endl;
		}
		else
		{
            cout << "place dirichlet bc at face: " << "false" << endl;
		}

		if (Dirichlet_Boundary_Condition)
		{
			cout << "Dirichlet Boundary Condition: " << "true" << endl;
		}
		else
		{
			cout << "Dirichlet Boundary Condition: " << "false" << endl;
		}

		if (Neumann_Boundary_Condition)
		{
			cout << "Neumann Boundary Condition: " << "true" << endl;
		}
		else
		{
			cout << "Neumann Boundary Condition: " << "false" << endl;
		}

		// Jump Condition Field
		jc_on_solution.Initialize(multithreading, base_grid, 1);
		jc_on_derivative.Initialize(multithreading, base_grid, 1);
		jc_on_solution.AssignAllValues((T)0);
		jc_on_derivative.AssignAllValues((T)0);

		// Initialize fields
		pressure_field.Initialize(multithreading, base_grid, 1);					// place 1 padding for distributed 
		boundary_condition_field.Initialize(multithreading, base_grid, 1);
		divergence_field.Initialize(multithreading, base_grid, 1);					// div field also requires ghost cells because restriction is interpolation
		projection_density_field.Initialize(multithreading, base_grid, 1);

		// Initialize Pressure field as 0
		pressure_field.AssignAllValues((T)0);

		// Assigning initial field values
		projection_density_field.array_for_this.AssignAllValues((T)1);

		// Initialize scan line algorithm
		scan_line_algorithm.Initialize(multithreading, boundary_condition_field.grid.ijk_res);

		// Initialize Poisson solvers - After you add the other projection solver, add this one more time
		poisson_solver.Initialize(multithreading, tolerance, max_iteration);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
		poisson_solver.Dirichlet_Boundary_Condition = Dirichlet_Boundary_Condition;
		poisson_solver.Neumann_Boundary_Condition = Neumann_Boundary_Condition;

		if (use_variable_density)
		{
			poisson_solver.use_variable_density = true;
		}
		poisson_solver.place_dirichlet_bc_at_face = place_dirichlet_bc_at_face;
		poisson_solver.air_water_simulation = air_water_simulation;
		poisson_solver.oil_water_simulation = oil_water_simulation;

		// Set up the helpful field
		one_over_density.Initialize(multithreading, projection_density_field.grid, 2);

		// Face Value of levelset
		levelset_x_c.Initialize(multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_y_c.Initialize(multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_z_c.Initialize(multithreading, water_velocity_field_mac_z.grid, 2);

		// Face Value of Density
		density_half_x.Initialize(multithreading, water_velocity_field_mac_x.grid, 2);
		density_half_y.Initialize(multithreading, water_velocity_field_mac_y.grid, 2);
		density_half_z.Initialize(multithreading, water_velocity_field_mac_z.grid, 2);

		// For Cylindrical Pipe
		boundary_phi_x.Initialize(multithreading, water_velocity_field_mac_x.grid, 2);
		boundary_phi_y.Initialize(multithreading, water_velocity_field_mac_y.grid, 2);
		boundary_phi_z.Initialize(multithreading, water_velocity_field_mac_z.grid, 2);

		// Ghost levelset
		scalar_field_ghost.Initialize(multithreading, water_levelset.grid, 2);
	}
	
public: // Solver
	void Solve(const int& thread_id, const T& dt)
	{
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				water_levelset.signed_distance_field.FillGhostCellsPeriodicInYDirection(thread_id, water_levelset.arr, false);
			}
			if (is_parallel)
			{
				water_levelset.signed_distance_field.FillGhostCellsPeriodicInXDirection(thread_id, water_levelset.arr, false);
			}
		}
		
		DetermineProjectionDensity(thread_id);
		DetermineDivergence(thread_id, dt);
		DetermineJumpConditionField(thread_id, dt);
		DeterminePressure(thread_id);
		UpdateVelocity(thread_id, dt);
		
		ofstream fout;
		if (Dirichlet_Boundary_Condition)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				fout.open("Pressure_diri_x");
			
				for (int i = pressure_field.i_start; i <= pressure_field.i_end; i++)
				{
					for (int j = pressure_field.j_start; j <= pressure_field.j_end; j++)
					{
						fout << pressure_field(i, j, pressure_field.k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close(); 

				fout.open("Pressure_diri_y");
			
				for (int k = pressure_field.k_start; k <= pressure_field.k_end; k++)
				{
					for (int i = pressure_field.i_start; i <= pressure_field.i_end; i++)
					{
						fout << pressure_field(i, pressure_field.j_end/2, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
				
				fout.open("Pressure_diri_z");
			
				for (int j = pressure_field.j_start; j <= pressure_field.j_end; j++)
				{
					for (int k = pressure_field.k_start; k <= pressure_field.k_end; k++)
					{
						fout << pressure_field(pressure_field.i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
			END_HEAD_THREAD_WORK; 
		}
		if (Neumann_Boundary_Condition)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				fout.open("Pressure_neum_x");
			
				for (int i = pressure_field.i_start; i <= pressure_field.i_end; i++)
				{
					for (int j = pressure_field.j_start; j <= pressure_field.j_end; j++)
					{
						fout << pressure_field(i, j, pressure_field.k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close(); 

				fout.open("Pressure_neum_y");
			
				for (int k = pressure_field.k_start; k <= pressure_field.k_end; k++)
				{
					for (int i = pressure_field.i_start; i <= pressure_field.i_end; i++)
					{
						fout << pressure_field(i, pressure_field.j_end/2, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
				
				fout.open("Pressure_neum_z");
			
				for (int j = pressure_field.j_start; j <= pressure_field.j_end; j++)
				{
					for (int k = pressure_field.k_start; k <= pressure_field.k_end; k++)
					{
						fout << pressure_field(pressure_field.i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
			END_HEAD_THREAD_WORK;
		}
		
		//
		//// Define true solution for the test
		//FIELD_STRUCTURE_3D<T> true_solution;
		//true_solution.Initialize(multithreading, pressure_field.grid, 2);
		//BEGIN_GRID_ITERATION_3D(true_solution.partial_grids[thread_id])
		//{
		//	const T x_coor = true_solution.x_min + i*true_solution.dx, y_coor = true_solution.y_min + j*true_solution.dy, z_coor = true_solution.z_min + k*true_solution.dz;
		//	true_solution(i, j, k) = exp(x_coor)*(sin(y_coor) + sin(z_coor));
		//	if (water_levelset(i, j, k) <= 0)
		//	{
		//		true_solution(i, j, k) = exp(-POW2(x_coor) - POW2(y_coor) - POW2(z_coor));
		//	}
		//	else
		//	{
		//		true_solution(i, j, k) = 0;
		//	}
		//}
		//END_GRID_ITERATION_3D;

		//fout.open("True");
		//for (int i = pressure_field.i_start; i <= pressure_field.i_end; i++)
		//{
		//	for (int j = pressure_field.j_start; j <= pressure_field.j_end; j++)
		//	{
		//		for (int k = pressure_field.k_start; k <= pressure_field.k_end; k++)
		//		{
		//			fout << i << " " << j << " " << k << " " << true_solution(i, j, k) << endl;;
		//		}
		//	}
		//}
		//fout.close();
	}

	void DetermineProjectionDensity(const int& thread_id)
	{
		if (!use_variable_density)
		{
			return;
		}

		if (air_water_simulation)
		{
			if (air_bubble_rising)
			{
				BEGIN_GRID_ITERATION_3D(projection_density_field.partial_grids[thread_id])
				{
					assert(projection_density_field.array_for_this.i_end == water_levelset.i_end);
				
					if (water_levelset(i, j, k) <= 0)
					{
						projection_density_field(i, j, k) = air_density;
					}
					else if (water_levelset(i, j, k) > 0)
					{
						projection_density_field(i, j, k) = water_density;
					}
				}
				END_GRID_ITERATION_3D;
			}
				
			if (water_drop)
			{
				BEGIN_GRID_ITERATION_3D(projection_density_field.partial_grids[thread_id])
				{
					assert(projection_density_field.array_for_this.i_end == water_levelset.i_end);

					if (water_levelset(i, j, k) <= 0)
					{
						projection_density_field(i, j, k) = water_density;
					}
					else if (water_levelset(i, j, k) > 0)
					{
						projection_density_field(i, j, k) = air_density;
					}
				}
				END_GRID_ITERATION_3D;
			}
		}
		
		if (oil_water_simulation)
		{
			BEGIN_GRID_ITERATION_3D(projection_density_field.partial_grids[thread_id])
			{
				assert(projection_density_field.array_for_this.i_end == water_levelset.i_end);

				if (dimensionless_form)
				{
					if (water_levelset(i, j, k) <= 0)
					{
						projection_density_field(i, j, k) = 1;
					}
					else if (water_levelset(i, j, k) > 0)
					{
						projection_density_field(i, j, k) = water_density/oil_density;
					}
				}
				else
				{
					if (water_levelset(i, j, k) <= 0)
					{
						projection_density_field(i, j, k) = oil_density;
					}
					else if (water_levelset(i, j, k) > 0)
					{
						projection_density_field(i, j, k) = water_density;
					}
				}
			}
			END_GRID_ITERATION_3D;
		}

		BEGIN_GRID_ITERATION_3D(one_over_density.partial_grids[thread_id])
		{
			one_over_density(i, j, k) = 1/projection_density_field(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void DetermineDivergence(const int& thread_id, const T& dt)
	{
		// Since we use multithreading technique, non-changing variable must be used - That's why we use ghost field for calculation. 
		if (air_water_simulation)
		{
			vector_field_mac_ghost_x.FillGhostCellsZero(thread_id, water_velocity_field_mac_x.array_for_this, true);
			vector_field_mac_ghost_y.FillGhostCellsZero(thread_id, water_velocity_field_mac_y.array_for_this, true);
			vector_field_mac_ghost_z.FillGhostCellsZero(thread_id, water_velocity_field_mac_z.array_for_this, true);

			BEGIN_GRID_ITERATION_3D(divergence_field.partial_grids[thread_id])
			{
				divergence_field(i, j, k) = one_over_dx*(vector_field_mac_ghost_x(i + 1, j, k) - vector_field_mac_ghost_x(i, j, k)) + one_over_dy*(vector_field_mac_ghost_y(i, j + 1, k) - vector_field_mac_ghost_y(i, j, k)) + one_over_dz*(vector_field_mac_ghost_z(i, j, k + 1) - vector_field_mac_ghost_z(i, j, k));
			}
			END_GRID_ITERATION_3D; 
		}
		
		if (oil_water_simulation)
		{
			SetupBoundaryConditionsForVelocityDivergence(thread_id);

			BEGIN_GRID_ITERATION_3D(divergence_field.partial_grids[thread_id])
			{
				if (boundary_levelset(i, j, k) > 0)
				{
					divergence_field(i, j, k) = (T)0;
				}
				else
				{
					T val_div = one_over_dx*(vector_field_mac_ghost_x(i + 1, j, k) - vector_field_mac_ghost_x(i, j, k));
					val_div += one_over_dy*(vector_field_mac_ghost_y(i, j + 1, k) - vector_field_mac_ghost_y(i, j, k));
					val_div += one_over_dz*(vector_field_mac_ghost_z(i, j, k + 1) - vector_field_mac_ghost_z(i, j, k));
				
					divergence_field(i, j, k) = val_div;
				}
			}
			END_GRID_ITERATION_3D; 
		}
		// Scaled by time step
		const T one_over_dt = 1/dt;
		
		BEGIN_GRID_ITERATION_3D(divergence_field.partial_grids[thread_id])
		{
			divergence_field(i, j, k) *= one_over_dt;	
		}
		END_GRID_ITERATION_3D;

		//// Poisson Solver Test
		//GRID_ITERATION_3D(divergence_field.partial_grids[thread_id])
		//{
		//	const T x_coord = divergence_field.grid.x_min + i*divergence_field.grid.dx, y_coord = divergence_field.grid.y_min + j*divergence_field.dy, z_coord = divergence_field.grid.z_min + k*divergence_field.dz;
		//	if (water_levelset(i, j, k) <= 0)
		//	{
		//		divergence_field(i, j, k) = (T)8*(POW2(x_coord) + POW2(y_coord) + POW2(z_coord) - (T)3/2)*exp(-POW2(x_coord)- POW2(y_coord) - POW2(z_coord));	
		//	}
		//	else
		//	{
		//		divergence_field(i, j, k) = 0;
		//	}
		//}
		//multithreading->Sync(thread_id);
		
		// TODO : Add additional divergence control from other controllers (such as explosion particles)
		//	      void AddDivergenceFromParticles(const int& thread_id);
	}

	void DetermineJumpConditionField(const int& thread_id, const T& dt)
	{
		water_levelset.ComputeCurvatures(thread_id);

		if (use_delta_function_formulation == true)
		{
			BEGIN_GRID_ITERATION_3D(jc_on_solution.partial_grids[thread_id])
			{
				if (CSF_model)
				{
					jc_on_solution(i, j, k) = (T)0;
				}
				else
				{
					if (dimensionless_form)
					{
						jc_on_solution(i, j, k) = (T)1/We*water_levelset.curvature(i, j, k);
						//jc_on_solution(i, j, k) = (T)0;
					}
					else
					{
						jc_on_solution(i, j, k) = surface_tension*water_levelset.curvature(i, j, k);
					}
				}
				
				// For the Poisson Problem Test
				/*const T x_coord = jc_on_solution.grid.x_min + i*jc_on_solution.dx, y_coord = jc_on_solution.grid.y_min + j*jc_on_solution.dy, z_coord = jc_on_solution.grid.z_min + k*jc_on_solution.dz; 
				jc_on_solution(i, j, k) = -exp(-POW2(x_coord) - POW2(y_coord) - POW2(z_coord));
				jc_on_derivative(i, j, k) = (T)8*((T)2*POW2(x_coord) + (T)2*POW2(y_coord) + (T)2*POW2(z_coord) - x_coord - y_coord - z_coord)*exp(-POW2(x_coord) - POW2(y_coord) - POW2(z_coord));	*/
			}
			END_GRID_ITERATION_3D;

			if (oil_water_simulation)
			{
				if (is_vertical)
				{
					jc_on_solution.FillGhostCellsPeriodicInYDirection(thread_id, jc_on_solution.array_for_this, true);
					multithreading->Sync(thread_id);
				}
				if (is_parallel)
				{
					jc_on_solution.FillGhostCellsPeriodicInXDirection(thread_id, jc_on_solution.array_for_this, true);
					multithreading->Sync(thread_id);
				}
			}
		}
	}

	void DeterminePressure(const int& thread_id)
	{
		if (use_previous_pressure != true)
		{
			pressure_field.AssignAllValuesGhost(thread_id, (T)0);
		}

		SetupBoundaryConditions(thread_id, pressure_field, boundary_condition_field, boundary_levelset, water_levelset);

		water_levelset.ComputeNormals(thread_id);
		
		BEGIN_GRID_ITERATION_3D(jc_on_derivative.partial_grids[thread_id])
		{
			jc_on_derivative(i, j, k) = (T)0;
		}
		END_GRID_ITERATION_3D;
		
		if (air_water_simulation)
		{
			one_over_density.FillGhostCellsFrom(thread_id, one_over_density.array_for_this, false);
			multithreading->Sync(thread_id);
		}
		
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				one_over_density.FillGhostCellsPeriodicInYDirection(thread_id, one_over_density.array_for_this, false);
				multithreading->Sync(thread_id);
			}
			if (is_parallel)
			{
				one_over_density.FillGhostCellsPeriodicInXDirection(thread_id, one_over_density.array_for_this, false);
				multithreading->Sync(thread_id);
			} 
		}

		boundary_levelset.ComputeNormals(thread_id);

		poisson_solver.Solve(thread_id, pressure_field, projection_density_field, boundary_condition_field, divergence_field, one_over_density, water_levelset, boundary_levelset, jc_on_solution, jc_on_derivative);
	}

	void UpdateVelocity(const int& thread_id, const T& dt)
	{
		UpdateVelocityByPressureGradientVariableDensity(thread_id, water_velocity_field_mac_x, dt);
		UpdateVelocityByPressureGradientVariableDensity(thread_id, water_velocity_field_mac_y, dt);
		UpdateVelocityByPressureGradientVariableDensity(thread_id, water_velocity_field_mac_z, dt);
	}
	
public: // Member Functions
	void SetProjectionType(const PROJECTION_TYPE& projection_type_input)
	{
		projection_type = projection_type_input;
	}

	void SetupBoundaryConditions(const int& thread_id, FIELD_STRUCTURE_3D<T>& pressure_input, FIELD_STRUCTURE_3D<int>& bc_input, const LEVELSET_3D& object_levelset_input, const LEVELSET_3D& water_levelset_input)
	{
		ARRAY_3D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_3D& grid(bc_input.grid);

		if (Dirichlet_Boundary_Condition)
		{
			if (air_water_simulation)
			{
				BEGIN_GRID_ITERATION_3D(bc_input.partial_grids_ghost[thread_id])
				{
					assert(pressure_input.i_end == water_levelset.arr.i_end && pressure_input.i_end == bc_array.i_end);
				
					//const T x_coor = grid.x_min + i*grid.dx, y_coor = grid.y_min + j*grid.dy, z_coor = grid.z_min + k*grid.dz;

					if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end || k < grid.k_start || k > grid.k_end)
					{
						bc_array(i, j, k) = BC_DIR;
						//pressure_input(i, j, k) = 0;
						if (i < grid.i_start)
						{
							pressure_input(i, j, k) = pressure_input(grid.i_start, j, k);
						}
						if (i > grid.i_end)
						{
							pressure_input(i, j, k) = pressure_input(grid.i_end, j, k);
						}
						if (j < grid.j_start)
						{
							pressure_input(i, j, k) = pressure_input(i, grid.j_start, k);
						}
						if (j > grid.j_end)
						{
							pressure_input(i, j, k) = pressure_input(i, grid.j_end, k);
						}
						if (k < grid.k_start)
						{
							pressure_input(i, j, k) = pressure_input(i, j, grid.k_start);
						}
						if (k > grid.k_end)
						{
							pressure_input(i, j, k) = pressure_input(i, j, grid.k_end);
						}
					}
					else
					{
						bc_array(i, j, k) = BC_FULL;
					}
				}	
				END_GRID_ITERATION_3D;
			}

			if (oil_water_simulation)
			{
				BEGIN_GRID_ITERATION_3D(bc_input.partial_grids_ghost[thread_id])
				{
					assert(pressure_input.i_end == water_levelset.arr.i_end && pressure_input.i_end == bc_array.i_end);
				
					//const T x_coor = grid.x_min + i*grid.dx, y_coor = grid.y_min + j*grid.dy, z_coor = grid.z_min + k*grid.dz;

					if (cylindrical_boundary)
					{
						if (is_vertical)
						{
							if (boundary_levelset(i, j, k) > (T)0)
							{
								bc_array(i, j, k) = BC_DIR;
								pressure_input(i, j, k) = 0;
							}
							else if (j < grid.j_start)
							{
								bc_array(i, j, k) = BC_DIR;
								pressure_input(i, j, k) = 0;
							}
							else if (j > grid.j_end)
							{
								bc_array(i, j, k) = BC_DIR;
								pressure_input(i, j, k) = 0;
							}
							else
							{
								bc_array(i, j, k) = BC_FULL;
							} 
						}
						if (is_parallel)
						{
							if (boundary_levelset(i, j, k) > (T)0)
							{
								bc_array(i, j, k) = BC_DIR;
								pressure_input(i, j, k) = 0;
							}
							else if (i < grid.i_start)
							{
								bc_array(i, j, k) = BC_DIR;
								pressure_input(i, j, k) = 0;
							}
							else if (i > grid.i_end)
							{
								bc_array(i, j, k) = BC_DIR;
								pressure_input(i, j, k) = 0;
							}
							else
							{
								bc_array(i, j, k) = BC_FULL;
							} 
						}
					}
				
					if (regular_boundary)
					{
						if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end || k < grid.k_start || k > grid.k_end)
						{
							bc_array(i, j, k) = BC_DIR;
							//pressure_input(i, j, k) = 0;
							if (i < grid.i_start)
							{
								pressure_input(i, j, k) = pressure_input(grid.i_start, j, k);
							}
							if (i > grid.i_end)
							{
								pressure_input(i, j, k) = pressure_input(grid.i_end, j, k);
							}
							if (j < grid.j_start)
							{
								pressure_input(i, j, k) = pressure_input(i, grid.j_start, k);
							}
							if (j > grid.j_end)
							{
								pressure_input(i, j, k) = pressure_input(i, grid.j_end, k);
							}
							if (k < grid.k_start)
							{
								pressure_input(i, j, k) = pressure_input(i, j, grid.k_start);
							}
							if (k > grid.k_end)
							{
								pressure_input(i, j, k) = pressure_input(i, j, grid.k_end);
							}
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					}
				}
				END_GRID_ITERATION_3D;
			}
		}
		
		if (Neumann_Boundary_Condition)
		{
			if (air_water_simulation)
			{
				BEGIN_GRID_ITERATION_3D(bc_input.partial_grids_ghost[thread_id])
				{
					assert(pressure_input.i_end == water_levelset.arr.i_end && pressure_input.i_end == bc_array.i_end);
				
					if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end || k < grid.k_start || k > grid.k_end)
					{
						bc_array(i, j, k) = BC_NEUM;
					}
					else
					{
						bc_array(i, j, k) = BC_FULL;
					}
				}
				END_GRID_ITERATION_3D;
			}
			if (oil_water_simulation)
			{
				BEGIN_GRID_ITERATION_3D(bc_input.partial_grids_ghost[thread_id])
				{
					assert(pressure_input.i_end == water_levelset.arr.i_end && pressure_input.i_end == bc_array.i_end);
				
					//const T x_coor = grid.x_min + i*grid.dx, y_coor = grid.y_min + j*grid.dy, z_coor = grid.z_min + k*grid.dz;

					if (cylindrical_boundary)
					{
						if (is_vertical)
						{
							if (boundary_levelset(i, j, k) > 0)
							{
								bc_array(i, j, k) = BC_NEUM;
							}
							else if (j < grid.j_start)
							{
								bc_array(i, j, k) = BC_NEUM;
							}
							else if (j > grid.j_end)
							{
								bc_array(i, j, k) = BC_NEUM;
							}
							else
							{
								bc_array(i, j, k) = BC_FULL;
							} 
						}
						if (is_parallel)
						{
							if (boundary_levelset(i, j, k) > 0)
							{
								bc_array(i, j, k) = BC_NEUM;
							}
							else if (i < grid.i_start)
							{
								bc_array(i, j, k) = BC_NEUM;
							}
							else if (i > grid.i_end)
							{
								bc_array(i, j, k) = BC_NEUM;
							}
							else
							{
								bc_array(i, j, k) = BC_FULL;
							} 
						}
					}
				
					if (regular_boundary)
					{
						if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end || k < grid.k_start || k > grid.k_end)
						{
							bc_array(i, j, k) = BC_NEUM;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					}
				}
				END_GRID_ITERATION_3D;
			}
		}
	}
				
	void UpdateVelocityByPressureGradientVariableDensity(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity_field, const T dt)
	{
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInYDirection(thread_id, water_levelset.arr, true);
			}
			if (is_parallel)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInXDirection(thread_id, water_levelset.arr, true);
			}
		}
		
		if (air_water_simulation)
		{
			scalar_field_ghost.FillGhostCellsContinuousDerivativesFrom(thread_id, water_levelset.arr, true);
		}
		
		if (velocity_field.is_x_component)
		{
			BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
			{
				levelset_x_c(i, j, k) = (T)0.5*(scalar_field_ghost(i, j, k) + scalar_field_ghost(i - 1, j, k));
			}
			END_GRID_ITERATION_3D;
			
			DetermineDensityField(thread_id, levelset_x_c, (T)1.5*base_grid.dx, density_half_x); 
		}
		
		if (velocity_field.is_y_component)
		{
			BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
			{
				levelset_y_c(i, j, k) = (T)0.5*(scalar_field_ghost(i, j, k) + scalar_field_ghost(i, j - 1, k));
			}
			END_GRID_ITERATION_3D;
			
			DetermineDensityField(thread_id, levelset_y_c, (T)1.5*base_grid.dy, density_half_y); 
		}

		if (velocity_field.is_z_component)
		{
			BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
			{
				levelset_z_c(i, j, k) = (T)0.5*(scalar_field_ghost(i, j, k) + scalar_field_ghost(i, j, k - 1));
			}
			END_GRID_ITERATION_3D;
			
			DetermineDensityField(thread_id, levelset_z_c, (T)1.5*base_grid.dz, density_half_z); 	
		}
					
		if (velocity_field.is_x_component)
		{
			max_velocity_x = 0;

			BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
			{
				if (velocity_field.fixed(i, j, k) == true)
				{
					continue;
				}

				if (boundary_condition_array(i, j, k) < 0)
				{
					continue;
				}

				T& velocity_ijk = velocity_field.array_for_this(i, j, k);

				const T density_ijk = projection_density_field(i, j, k);
				const T one_over_density_ijk = 1/density_ijk;
				
				const T levelset_ijk = water_levelset(i, j, k);
				const T levelset_ijk_l = water_levelset(i - 1, j ,k);
				
				const T boundary_levelset_ijk = boundary_levelset(i, j, k);
				const T boundary_levelset_ijk_l = boundary_levelset(i - 1, j, k); 

				const T jump_condition = jc_on_solution(i, j, k);
				const T jump_condition_l = jc_on_solution(i - 1, j, k);

				const T jump_condition_gamma = (jump_condition*abs(levelset_ijk_l) + jump_condition_l*abs(levelset_ijk))/(abs(levelset_ijk_l) + abs(levelset_ijk));
				
				const T one_over_density_half_x = 1/density_half_x(i, j, k);
				
				if (is_vertical)
				{
					if (boundary_levelset_ijk*boundary_levelset_ijk_l >= 0)
					{
						if (levelset_ijk*levelset_ijk_l > 0)
						{
							velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i - 1, j , k))*one_over_density_ijk*one_over_dx;
						}
						else
						{
							if (levelset_ijk <= 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i - 1, j , k) + jump_condition_gamma)*one_over_density_half_x*one_over_dx;
							}
							else if (levelset_ijk > 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i - 1, j , k) - jump_condition_gamma)*one_over_density_half_x*one_over_dx;
							}
						} 
					} 
				}
				
				if (is_parallel)
				{
					if (boundary_levelset_ijk*boundary_levelset_ijk_l >= 0)
					{
						if (levelset_ijk*levelset_ijk_l > 0)
						{
							velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i - 1, j , k))*one_over_density_ijk*one_over_dx;
						}
						else
						{
							if (levelset_ijk <= 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i - 1, j , k) + jump_condition_gamma)*one_over_density_half_x*one_over_dx;
							}
							else if (levelset_ijk > 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i - 1, j , k) - jump_condition_gamma)*one_over_density_half_x*one_over_dx;
							}
						} 
					} 
				}

				max_velocity_x = MAX(max_velocity_x, abs(velocity_ijk));
			}
			END_GRID_ITERATION_MAX_3D(max_velocity_x);

			if (oil_water_simulation)
			{
				if (is_vertical)
				{
					BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
					{
						velocity_field(i, velocity_field.grid.j_end, k) = velocity_field(i, velocity_field.grid.j_start, k);
					}
					END_GRID_ITERATION_3D;  
				}
				if (is_parallel)
				{
					BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
					{
						velocity_field(velocity_field.grid.i_end, j, k) = velocity_field(velocity_field.grid.i_start + 1, j, k);
						velocity_field(velocity_field.grid.i_start, j, k) = velocity_field(velocity_field.grid.i_end - 1, j, k);
					}
					END_GRID_ITERATION_3D;  
				}
			}
		}
		
		if (velocity_field.is_y_component)
		{
			max_velocity_y = 0;

			BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
			{
				if (velocity_field.fixed(i, j, k) == true)
				{
					continue;
				}

				if (boundary_condition_array(i, j, k) < 0)
				{
					continue;
				}

				T& velocity_ijk = velocity_field.array_for_this(i, j, k);
				
				const T density_ijk = projection_density_field(i, j, k);
				const T one_over_density_ijk = (T)1/density_ijk;
				
				const T levelset_ijk = water_levelset(i, j, k);
				const T levelset_ijk_b = water_levelset(i, j - 1, k);
				
				const T boundary_levelset_ijk = boundary_levelset(i, j, k);
				const T boundary_levelset_ijk_b = boundary_levelset(i, j - 1, k); 

				const T jump_condition = jc_on_solution(i, j, k);
				const T jump_condition_b = jc_on_solution(i, j - 1, k);

				const T jump_condition_gamma = (jump_condition*abs(levelset_ijk_b) + jump_condition_b*abs(levelset_ijk))/(abs(levelset_ijk_b) + abs(levelset_ijk));
				
				const T one_over_density_half_y = 1/density_half_y(i, j, k);
				
				if (is_vertical)
				{
					if ((boundary_levelset_ijk*boundary_levelset_ijk_b >= 0) && (j != velocity_field.grid.j_start) && (j <= (velocity_field.grid.j_end - 1)))
					{
						if (levelset_ijk*levelset_ijk_b > 0)
						{
							velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j - 1, k))*one_over_density_ijk*one_over_dy;
						}
						else
						{
							if (levelset_ijk <= 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j - 1, k) + jump_condition_gamma)*one_over_density_half_y*one_over_dy;
							}
							else if (levelset_ijk > 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j - 1, k) - jump_condition_gamma)*one_over_density_half_y*one_over_dy;
							}
						}
					} 
				}
				
				if (is_parallel)
				{
					if (boundary_levelset_ijk*boundary_levelset_ijk_b >= 0)
					{
						if (levelset_ijk*levelset_ijk_b > 0)
						{
							velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j - 1, k))*one_over_density_ijk*one_over_dy;
						}
						else
						{
							if (levelset_ijk < 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j - 1, k) + jump_condition_gamma)*one_over_density_half_y*one_over_dy;
							}
							else if (levelset_ijk > 0)
							{
								velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j - 1, k) - jump_condition_gamma)*one_over_density_half_y*one_over_dy;
							}
						}
					} 
				}

				max_velocity_y = MAX(max_velocity_y, abs(velocity_ijk));
			}
			END_GRID_ITERATION_MAX_3D(max_velocity_y);
			
			if (oil_water_simulation)
			{
				if (is_vertical)
				{
					BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
					{
						velocity_field(i, velocity_field.grid.j_end, k) = velocity_field(i, velocity_field.grid.j_start + 1, k);
						velocity_field(i, velocity_field.grid.j_start, k) = velocity_field(i, velocity_field.grid.j_end - 1, k);
					}
					END_GRID_ITERATION_3D;  
				}
				if (is_parallel)
				{
					BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
					{
						velocity_field(velocity_field.grid.i_end, j, k) = velocity_field(velocity_field.grid.i_start, j, k);
					}
					END_GRID_ITERATION_3D;  
				}
			}
		}
		
		if (velocity_field.is_z_component)
		{
			max_velocity_z = 0;

			BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
			{
				if (velocity_field.fixed(i, j, k) == true)
				{
					continue;
				}

				if (boundary_condition_array(i, j, k) < 0)
				{
					continue;
				}

				T& velocity_ijk = velocity_field.array_for_this(i, j, k);
	
				const T density_ijk = projection_density_field(i, j, k);
				const T one_over_density_ijk = (T)1/density_ijk;
				
				const T levelset_ijk = water_levelset(i, j, k);
				const T levelset_ijk_d = water_levelset(i, j, k - 1);
			
				const T boundary_levelset_ijk = boundary_levelset(i, j, k);
				const T boundary_levelset_ijk_d = boundary_levelset(i, j, k - 1); 

				const T jump_condition = jc_on_solution(i, j, k);
				const T jump_condition_d = jc_on_solution(i, j, k - 1);

				const T jump_condition_gamma = (jump_condition*abs(levelset_ijk_d) + jump_condition_d*abs(levelset_ijk))/(abs(levelset_ijk_d) + abs(levelset_ijk));
				
				const T one_over_density_half_z = 1/density_half_z(i, j, k);
				
				if (boundary_levelset_ijk*boundary_levelset_ijk_d >= 0)
				{
					if (levelset_ijk*levelset_ijk_d > 0)
					{
						velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j, k - 1))*one_over_density_ijk*one_over_dz;
					}
					else
					{
						if (levelset_ijk <= 0)
						{
							velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j, k - 1) + jump_condition_gamma)*one_over_density_half_z*one_over_dz;
						}
						else if (levelset_ijk > 0)
						{
							velocity_ijk -= dt*(pressure_field(i, j, k) - pressure_field(i, j, k - 1) - jump_condition_gamma)*one_over_density_half_z*one_over_dz;
						}
					} 
				}
				max_velocity_z = MAX(max_velocity_z, abs(velocity_ijk));
			}
			END_GRID_ITERATION_MAX_3D(max_velocity_z);

			if (oil_water_simulation)
			{
				if (is_vertical)
				{
					BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
					{
						velocity_field(i, velocity_field.grid.j_end, k) = velocity_field(i, velocity_field.grid.j_start, k);
					}
					END_GRID_ITERATION_3D;  
				}
				if (is_parallel)
				{
					BEGIN_GRID_ITERATION_3D(velocity_field.partial_grids[thread_id])
					{
						velocity_field(velocity_field.grid.i_end, j, k) = velocity_field(velocity_field.grid.i_start, j, k);
					}
					END_GRID_ITERATION_3D;  
				}
			}
		}
	}

public: // Sub Functions
	void HeavisideFunction(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& heaviside)
	{
		T one_over_epsilon = (T)1/epsilon, one_over_pi = (T)1/PI;

		GRID_ITERATION_3D(phi.partial_grids[thread_id])
		{
			if (phi(i, j, k) < - epsilon)
			{
				heaviside(i, j, k) = 0;
			}
			else if (phi(i, j, k) > epsilon)
			{
				heaviside(i, j, k) = 1;
			}
			else
			{
				heaviside(i, j, k) = (T)0.5 + phi(i, j, k)*one_over_epsilon*(T)0.5 + (T)0.5*one_over_pi*sin(PI*phi(i, j, k)*one_over_epsilon);
			}
		}
		
		multithreading->Sync(thread_id);
	}

	void DetermineDensityField(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& density)
	{
		FIELD_STRUCTURE_3D<T> heaviside_phi;
		heaviside_phi.Initialize(multithreading, density.grid, 2);

		HeavisideFunction(thread_id, phi, epsilon, heaviside_phi);
		
		multithreading->Sync(thread_id);
		
		if (air_water_simulation)
		{
			GRID_ITERATION_3D(phi.partial_grids[thread_id])
			{
				if (air_bubble_rising)
				{
					density(i, j, k) = air_density + (water_density - air_density)*heaviside_phi(i, j, k);
				}
				if (water_drop)
				{
					density(i, j, k) = water_density + (air_density - water_density)*heaviside_phi(i, j, k);
				}
			}
		}
		
		if (oil_water_simulation)
		{
			GRID_ITERATION_3D(phi.partial_grids[thread_id])
			{
				if (dimensionless_form)
				{
					//density(i, j, k) = oil_density/water_density + ((T)1 - oil_density/water_density)*heaviside_phi(i, j, k);
					density(i, j, k) = 1 + (water_density/oil_density - (T)1)*heaviside_phi(i, j, k);
				}
				else
				{
					density(i, j, k) = oil_density + (water_density - oil_density)*heaviside_phi(i, j, k);
				}
				
			}
		}

		density.FillGhostCellsFrom(thread_id, density.array_for_this, true);

		multithreading->Sync(thread_id);
	}

	void SetupBoundaryConditionsForVelocityDivergence(const int& thread_id)
	{
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				BEGIN_GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_x.x_min + i*boundary_phi_x.dx, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
					boundary_phi_x(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_y.x_min + i*boundary_phi_y.dx, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
					boundary_phi_y(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_z.x_min + i*boundary_phi_z.dx, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
					boundary_phi_z(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
				}
				END_GRID_ITERATION_3D; 
			}
			if (is_parallel)
			{
				BEGIN_GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_x.y_min + j*boundary_phi_x.dy, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
					boundary_phi_x(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_y.y_min + j*boundary_phi_y.dy, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
					boundary_phi_y(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_z.y_min + j*boundary_phi_z.dy, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
					boundary_phi_z(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
				}
				END_GRID_ITERATION_3D; 
			}

			if (is_vertical)
			{
				vector_field_mac_ghost_x.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
				vector_field_mac_ghost_y.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
				vector_field_mac_ghost_z.FillGhostCellsPeriodicInYDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);
			
				BEGIN_GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					if (boundary_phi_x(i, j, k) > 0)
					{	
						if ((boundary_phi_x(i - 1, j, k) > 0) && (boundary_phi_x(i + 1, j, k) > 0) && (boundary_phi_x(i, j, k - 1) > 0) && (boundary_phi_x(i, j, k + 1) > 0))
						{
							vector_field_mac_ghost_x(i, j, k) = (T)0;
						}
						else if (boundary_phi_x(i - 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_x(i - 1, j, k))/(abs(boundary_phi_x(i - 1, j, k)) + abs(boundary_phi_x(i, j, k)));
							vector_field_mac_ghost_x(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_x(i - 1, j, k);
						}
						else if (boundary_phi_x(i + 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_x(i + 1, j, k))/(abs(boundary_phi_x(i + 1, j, k)) + abs(boundary_phi_x(i, j, k)));
							vector_field_mac_ghost_x(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_x(i + 1, j, k);
						}
						else
						{
							vector_field_mac_ghost_x(i, j, k) = (T)0;
						}
					}
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					if (boundary_phi_z(i, j, k) > 0)
					{
						if ((boundary_phi_z(i - 1, j, k) > 0) && (boundary_phi_z(i + 1, j, k) > 0) && (boundary_phi_z(i, j, k - 1) > 0) && (boundary_phi_z(i, j, k + 1) > 0))
						{
							vector_field_mac_ghost_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_z(i, j, k - 1) < 0)
						{
							T theta = abs(boundary_phi_z(i, j, k - 1))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
							vector_field_mac_ghost_z(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_z(i, j, k - 1);
						}
						else if (boundary_phi_z(i, j, k + 1) < 0)
						{
							T theta = abs(boundary_phi_z(i, j, k + 1))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
							vector_field_mac_ghost_z(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_z(i, j, k + 1);
						}
						else
						{
							vector_field_mac_ghost_z(i, j, k) = (T)0;
						}
					}	
				}
				END_GRID_ITERATION_3D;
			}
			if (is_parallel)
			{
				// Periodic Boundary Condition
				vector_field_mac_ghost_x.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_x.array_for_this, boundary_phi_x.array_for_this, true);
				vector_field_mac_ghost_y.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_y.array_for_this, boundary_phi_y.array_for_this, true);
				vector_field_mac_ghost_z.FillGhostCellsPeriodicInXDirection(thread_id, water_velocity_field_mac_z.array_for_this, boundary_phi_z.array_for_this, true);
			
				BEGIN_GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					if (boundary_phi_z(i, j, k) > 0)
					{
						if ((boundary_phi_z(i, j - 1, k) > 0) && (boundary_phi_z(i, j + 1, k) > 0) && (boundary_phi_z(i, j, k - 1) > 0) && (boundary_phi_z(i, j, k + 1) > 0))
						{
							vector_field_mac_ghost_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_z(i, j, k - 1) < 0)
						{
							T theta = abs(boundary_phi_z(i, j, k - 1))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
							vector_field_mac_ghost_z(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_z(i, j, k - 1);
						}
						else if (boundary_phi_z(i, j, k + 1) < 0)
						{
							T theta = abs(boundary_phi_z(i, j, k + 1))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
							vector_field_mac_ghost_z(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_z(i, j, k + 1);
						}
						else
						{
							vector_field_mac_ghost_z(i, j, k) = (T)0;
						}
					}	
				}
				END_GRID_ITERATION_3D;
			
				BEGIN_GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					if (boundary_phi_y(i, j, k) > 0)
					{
						if ((boundary_phi_y(i, j - 1, k) > 0) && (boundary_phi_y(i, j + 1, k) > 0) && (boundary_phi_y(i, j, k - 1) > 0) && (boundary_phi_y(i, j, k + 1) > 0))
						{
							vector_field_mac_ghost_y(i, j, k) = (T)0;
						}
						else if (boundary_phi_y(i, j - 1, k) < 0)
						{
							T theta = abs(boundary_phi_y(i, j - 1, k))/(abs(boundary_phi_y(i, j - 1, k)) + abs(boundary_phi_y(i, j, k)));
							vector_field_mac_ghost_y(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_y(i, j - 1, k);
						}
						else if (boundary_phi_y(i, j + 1, k) < 0)
						{
							T theta = abs(boundary_phi_y(i, j + 1, k))/(abs(boundary_phi_y(i, j + 1, k)) + abs(boundary_phi_y(i, j, k)));
							vector_field_mac_ghost_y(i, j, k) = (theta - (T)1)/theta*vector_field_mac_ghost_y(i, j + 1, k);
						}
						else
						{
							vector_field_mac_ghost_y(i, j, k) = (T)0;
						}
					}	
				}
				END_GRID_ITERATION_3D;
			}
		}
	}
};

