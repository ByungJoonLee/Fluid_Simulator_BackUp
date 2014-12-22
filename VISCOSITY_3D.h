#pragma once

#include "POISSON_SOLVER_3D.h"
#include "LEVELSET_3D.h"

class VISCOSITY_3D
{
public: 
	MULTITHREADING&				multithreading;

public:
	// Water Levelset
	LEVELSET_3D&				water_levelset;
	LEVELSET_3D&				boundary_levelset;
	FIELD_STRUCTURE_3D<T>&		water_signed_distance_field;

	GRID_STRUCTURE_3D&			base_grid;

	// For MAC Grid
	FIELD_STRUCTURE_3D<T>&		water_velocity_field_mac_x;
	FIELD_STRUCTURE_3D<T>&		water_velocity_field_mac_y;
	FIELD_STRUCTURE_3D<T>&		water_velocity_field_mac_z;
	
	FIELD_STRUCTURE_3D<T>&		scalar_field_ghost;
	FIELD_STRUCTURE_3D<T>&		velocity_field_mac_ghost_x;
	FIELD_STRUCTURE_3D<T>&		velocity_field_mac_ghost_y;
	FIELD_STRUCTURE_3D<T>&		velocity_field_mac_ghost_z;

	// For boundary condition
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_x_x;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_x_y;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_x_z;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_y_x;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_y_y;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_y_z;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_z_x;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_z_y;
	FIELD_STRUCTURE_3D<T>		velocity_field_mac_ghost_z_z;

	// Viscosity and Density
	T							water_viscosity, air_viscosity, oil_viscosity;
	T							water_density, air_density, oil_density;

	// Options for Bubble
	bool						air_bubble_rising, water_drop;
	
	// Options for Cylinder
	bool						regular_boundary, cylindrical_boundary;

	// Options for Viscosity
	bool						use_delta_function_formulation;
	
    // Options For Simulation
    bool						air_water_simulation;
    bool						oil_water_simulation;
	bool						dimensionless_form;

	// Options for Pipe
	bool						is_vertical, is_parallel;

public: // Sub Solver - For the semi-implicit method
	bool						semi_implicit_approach;

	enum POISSON_SOLVER_TYPE	poisson_solver_type;
	POISSON_SOLVER_3D			poisson_solver;

	FIELD_STRUCTURE_3D<int>		boundary_condition_field_x;
	FIELD_STRUCTURE_3D<int>		boundary_condition_field_y;
	FIELD_STRUCTURE_3D<int>		boundary_condition_field_z;

	T							tolerance;
	int							max_iteration;

	
public: // Sub Variables for Delta Function approach
	// x-component 
	// Density
	FIELD_STRUCTURE_3D<T>		density_half_x;
	
	// Levelset
	FIELD_STRUCTURE_3D<T>		levelset_x_half;
	FIELD_STRUCTURE_3D<T>		levelset_cu_x;
	FIELD_STRUCTURE_3D<T>		levelset_cd_x;
	FIELD_STRUCTURE_3D<T>		levelset_ct_x;
	FIELD_STRUCTURE_3D<T>		levelset_cb_x;
	FIELD_STRUCTURE_3D<T>		levelset_r_x;
	FIELD_STRUCTURE_3D<T>		levelset_c_x;
	FIELD_STRUCTURE_3D<T>		levelset_r_d_x;
	FIELD_STRUCTURE_3D<T>		levelset_c_d_x;
	FIELD_STRUCTURE_3D<T>		levelset_r_u_x;
	FIELD_STRUCTURE_3D<T>		levelset_c_u_x;

	// Viscosity
	FIELD_STRUCTURE_3D<T>		viscosity_cu_x;
	FIELD_STRUCTURE_3D<T>		viscosity_cd_x;
	FIELD_STRUCTURE_3D<T>		viscosity_r_x;
	FIELD_STRUCTURE_3D<T>		viscosity_c_x;
	FIELD_STRUCTURE_3D<T>		viscosity_ct_x;
	FIELD_STRUCTURE_3D<T>		viscosity_cb_x;
	FIELD_STRUCTURE_3D<T>		viscosity_r_d_x;
	FIELD_STRUCTURE_3D<T>		viscosity_c_d_x;
	FIELD_STRUCTURE_3D<T>		viscosity_r_u_x;
	FIELD_STRUCTURE_3D<T>		viscosity_c_u_x;

	// y-component
	// Density
	FIELD_STRUCTURE_3D<T>		density_half_y;
	
	// Levelset
	FIELD_STRUCTURE_3D<T>		levelset_y_half;
	FIELD_STRUCTURE_3D<T>		levelset_cu_y;
	FIELD_STRUCTURE_3D<T>		levelset_cl_y;
	FIELD_STRUCTURE_3D<T>		levelset_ct_y;
	FIELD_STRUCTURE_3D<T>		levelset_cb_y;
	FIELD_STRUCTURE_3D<T>		levelset_u_y;
	FIELD_STRUCTURE_3D<T>		levelset_c_y;

	// Viscosity
	FIELD_STRUCTURE_3D<T>		viscosity_cu_y;
	FIELD_STRUCTURE_3D<T>		viscosity_cl_y;
	FIELD_STRUCTURE_3D<T>		viscosity_ct_y;
	FIELD_STRUCTURE_3D<T>		viscosity_cb_y;
	FIELD_STRUCTURE_3D<T>		viscosity_u_y;
	FIELD_STRUCTURE_3D<T>		viscosity_c_y;

	// z-component
	// Density
	FIELD_STRUCTURE_3D<T>		density_half_z;

	// Levelset
	FIELD_STRUCTURE_3D<T>		levelset_z_half;
	FIELD_STRUCTURE_3D<T>		levelset_ct_z;
	FIELD_STRUCTURE_3D<T>		levelset_cb_z;
	FIELD_STRUCTURE_3D<T>		levelset_cu_z;
	FIELD_STRUCTURE_3D<T>		levelset_cd_z;
	FIELD_STRUCTURE_3D<T>		levelset_t_z;
	FIELD_STRUCTURE_3D<T>		levelset_c_z;

	// Viscosity
	FIELD_STRUCTURE_3D<T>		viscosity_ct_z;
	FIELD_STRUCTURE_3D<T>		viscosity_cb_z;
	FIELD_STRUCTURE_3D<T>		viscosity_cu_z;
	FIELD_STRUCTURE_3D<T>		viscosity_cd_z;
	FIELD_STRUCTURE_3D<T>		viscosity_t_z;
	FIELD_STRUCTURE_3D<T>		viscosity_c_z;

	// For Semi-Implicit Method
	FIELD_STRUCTURE_3D<T>		explicit_term_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x;
	FIELD_STRUCTURE_3D<T>		explicit_term_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y;
	FIELD_STRUCTURE_3D<T>		explicit_term_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z;

	// For Cylindrical Pipe
	FIELD_STRUCTURE_3D<T>		boundary_phi_x;
	FIELD_STRUCTURE_3D<T>		boundary_phi_y;
	FIELD_STRUCTURE_3D<T>		boundary_phi_z;

	// Speedup Variable
	T							one_over_dx, one_over_dy, one_over_dz;

public: // Dimensionless Variable
	T							Re;
	T							m, eta;

public: // Pipe Radius
	T							a;

public: // Constructor and Destructor
	VISCOSITY_3D(MULTITHREADING& multithreading_input, LEVELSET_3D& water_levelset_input, LEVELSET_3D& boundary_levelset_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_x_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_y_input, FIELD_STRUCTURE_3D<T>& water_velocity_field_mac_z_input,
				 FIELD_STRUCTURE_3D<T>& scalar_field_ghost_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_ghost_x_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_ghost_y_input, FIELD_STRUCTURE_3D<T>& velocity_field_mac_ghost_z_input)
				 : multithreading(multithreading_input), water_levelset(water_levelset_input), boundary_levelset(boundary_levelset_input), water_signed_distance_field(water_levelset_input.signed_distance_field), base_grid(water_levelset.grid), 
				 water_velocity_field_mac_x(water_velocity_field_mac_x_input), water_velocity_field_mac_y(water_velocity_field_mac_y_input), water_velocity_field_mac_z(water_velocity_field_mac_z_input), 
				 scalar_field_ghost(scalar_field_ghost_input), velocity_field_mac_ghost_x(velocity_field_mac_ghost_x_input), velocity_field_mac_ghost_y(velocity_field_mac_ghost_y_input), velocity_field_mac_ghost_z(velocity_field_mac_ghost_z_input),
                 air_density((T)0), water_density((T)0), oil_density((T)0), air_viscosity((T)0), water_viscosity((T)0), oil_viscosity((T)0), air_bubble_rising(false), water_drop(false), air_water_simulation(false), oil_water_simulation(false), semi_implicit_approach(false), dimensionless_form(false),
				 regular_boundary(false), cylindrical_boundary(false), is_vertical(false), is_parallel(false)
	{}

	~VISCOSITY_3D(void)
	{}

public: // Initialization Function
	void InitializeFromScriptBlock(const SCRIPT_BLOCK& script_block)
	{
		// Viscosity 
		if (air_water_simulation)
		{
            water_viscosity = script_block.GetFloat("water_viscosity", (T)1);
		    air_viscosity = script_block.GetFloat("air_viscosity", (T)1);
		}
		if (oil_water_simulation)
		{
            water_viscosity = script_block.GetFloat("water_viscosity", (T)1);
		    oil_viscosity = script_block.GetFloat("oil_viscosity", (T)1);
		}
		
		// Density
		if (air_water_simulation)
		{
            water_density = script_block.GetFloat("water_density", (T)1);
		    air_density = script_block.GetFloat("air_density", (T)1);
		}
		if (oil_water_simulation)
		{
            water_density = script_block.GetFloat("water_density", (T)1);
		    oil_density = script_block.GetFloat("oil_density", (T)1);
		}
		
		// Options for Bubble
		if (air_water_simulation)
		{
            air_bubble_rising = script_block.FindBlock("PROJECTION_WATER").GetBoolean("air_bubble_rising", (bool)false);
		    water_drop = script_block.FindBlock("PROJECTION_WATER").GetBoolean("water_drop", (bool)false);
		}
		
		// Options for Cylinder
		if (oil_water_simulation)
		{
			regular_boundary = script_block.FindBlock("PROJECTION_WATER").GetBoolean("regular_boundary", (bool)false);
			cylindrical_boundary = script_block.FindBlock("PROJECTION_WATER").GetBoolean("cylindrical_boundary", (bool)false);
		}

		// Options for Viscosity
		use_delta_function_formulation = script_block.GetBoolean("use_delta_function_formulation", (bool)false);
		
		// Semi-Implicit Approach
		semi_implicit_approach = script_block.GetBoolean("semi_implicit_approach", (bool)false);

		// x-component 
		// Density
		density_half_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		
		// Levelset
		levelset_x_half.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_cu_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_cd_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_ct_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_cb_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_r_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_c_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_r_d_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_c_d_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_r_u_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		levelset_c_u_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);

		// Viscosity
		viscosity_cu_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_cd_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_ct_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_cb_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_r_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_c_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_r_d_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_c_d_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_r_u_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		viscosity_c_u_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);

		// y-component
		// Density
		density_half_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);

		// Levelset
		levelset_y_half.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_cu_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_cl_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_ct_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_cb_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_u_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		levelset_c_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);

		// Viscosity
		viscosity_cu_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		viscosity_cl_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		viscosity_ct_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		viscosity_cb_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		viscosity_u_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		viscosity_c_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);

		// z-componet
		// Density
		density_half_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);

		// Levelset
		levelset_z_half.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		levelset_cu_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		levelset_cd_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		levelset_ct_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		levelset_cb_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		levelset_t_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		levelset_c_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);

		// Viscosity
		viscosity_cu_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		viscosity_cd_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		viscosity_ct_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		viscosity_cb_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		viscosity_t_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		viscosity_c_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);

		// For Cylindrical Pipe
		boundary_phi_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		boundary_phi_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		boundary_phi_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);

		// Speedup variables
		one_over_dx = base_grid.one_over_dx;
		one_over_dy = base_grid.one_over_dy;
		one_over_dz = base_grid.one_over_dz;

		// For the semi-implicit solver
		tolerance = script_block.FindBlock("PROJECTION_WATER").GetFloat("tolerance", (T)1e-4);
		max_iteration = script_block.FindBlock("PROJECTION_WATER").GetInteger("max_iteration", 30);

		const char* poisson_solver_type_input = script_block.FindBlock("PROJECTION_WATER").GetString("poisson_solver_type", "NULL");
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

		// Initialize the Poisson Solver
		poisson_solver.Initialize(&multithreading, tolerance, max_iteration);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
		
		boundary_condition_field_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 1);
		boundary_condition_field_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 1);
		boundary_condition_field_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 1);

		// Initialize the Semi-Implicit solver
		// x-component
		explicit_term_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		coef_1_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		coef_2_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		coef_3_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		coef_4_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		coef_5_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		coef_6_x.Initialize(&multithreading, water_velocity_field_mac_x.grid, 2);
		
		// y-component
		explicit_term_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		coef_1_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		coef_2_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		coef_3_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		coef_4_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		coef_5_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		coef_6_y.Initialize(&multithreading, water_velocity_field_mac_y.grid, 2);
		
		// z-component
		explicit_term_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		coef_1_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		coef_2_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		coef_3_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		coef_4_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		coef_5_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);
		coef_6_z.Initialize(&multithreading, water_velocity_field_mac_z.grid, 2);

		if (semi_implicit_approach)
		{
			cout << "--------------SEMI-IMPLICIT--------------" << endl;
			cout << "tolerance: " << tolerance << endl;
			cout << "max iteration: " << max_iteration << endl;
		}
		else
		{
			cout << "Semi-implicit for Viscosity: false" << endl;
		}
				
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

		// For boundary condition
		velocity_field_mac_ghost_x_x.Initialize(&multithreading, velocity_field_mac_ghost_x.grid, 3, false, true);
		velocity_field_mac_ghost_x_y.Initialize(&multithreading, velocity_field_mac_ghost_x.grid, 3, false, true);
		velocity_field_mac_ghost_x_z.Initialize(&multithreading, velocity_field_mac_ghost_x.grid, 3, false, true);
		velocity_field_mac_ghost_y_x.Initialize(&multithreading, velocity_field_mac_ghost_y.grid, 3, false, true);
		velocity_field_mac_ghost_y_y.Initialize(&multithreading, velocity_field_mac_ghost_y.grid, 3, false, true);
		velocity_field_mac_ghost_y_z.Initialize(&multithreading, velocity_field_mac_ghost_y.grid, 3, false, true);
		velocity_field_mac_ghost_z_x.Initialize(&multithreading, velocity_field_mac_ghost_z.grid, 3, false, true);
		velocity_field_mac_ghost_z_y.Initialize(&multithreading, velocity_field_mac_ghost_z.grid, 3, false, true);
		velocity_field_mac_ghost_z_z.Initialize(&multithreading, velocity_field_mac_ghost_z.grid, 3, false, true);
		
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

public: // Main Function
	void ApplyViscosity(const int& thread_id, const T& epsilon_for_mollification, const T& dt)
	{
		if (use_delta_function_formulation)
		{
			if (air_water_simulation)
			{
				// Speed Up Variable
				const ARRAY_3D<T>& water_array(water_levelset.arr);
			
				const ARRAY_3D<T>& velocity_x(water_velocity_field_mac_x.array_for_this);
				const ARRAY_3D<T>& velocity_y(water_velocity_field_mac_y.array_for_this);
				const ARRAY_3D<T>& velocity_z(water_velocity_field_mac_z.array_for_this);

				// Define the coefficient function -- coefficients of x-velocity
				GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
				{
					levelset_x_half(i, j, k) = (T)0.5*(water_array(i - 1, j, k) + water_array(i, j, k));
					levelset_cu_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j + 1, k) + water_array(i, j + 1, k));
					levelset_cd_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j - 1, k) + water_array(i, j - 1, k));
					levelset_ct_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k + 1) + water_array(i, j, k + 1));
					levelset_cb_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k - 1) + water_array(i, j, k - 1));
					levelset_r_x(i, j, k) = water_array(i, j, k);
					levelset_c_x(i, j, k) = water_array(i - 1, j, k);
				}
				
				multithreading.Sync(thread_id);
	
				DetermineDensityField(thread_id, levelset_x_half, epsilon_for_mollification, density_half_x);
				DetermineViscosityField(thread_id, levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
				DetermineViscosityField(thread_id, levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
				DetermineViscosityField(thread_id, levelset_r_x, epsilon_for_mollification, viscosity_r_x);
				DetermineViscosityField(thread_id, levelset_c_x, epsilon_for_mollification, viscosity_c_x);
				DetermineViscosityField(thread_id, levelset_ct_x, epsilon_for_mollification, viscosity_ct_x);
				DetermineViscosityField(thread_id, levelset_cb_x, epsilon_for_mollification, viscosity_cb_x);
				
				GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
				{
					T one_over_density_half_x = (T)1/density_half_x(i, j, k);
					T coef = dt*one_over_density_half_x;
					T first_update = (T)2*coef*one_over_dx*(viscosity_r_x(i, j, k)*one_over_dx*(velocity_x(i + 1, j, k) - velocity_x(i, j, k)) - viscosity_c_x(i, j, k)*one_over_dy*(velocity_x(i, j, k) - velocity_x(i - 1, j, k)));
					T second_update = coef*one_over_dy*(viscosity_cu_x(i, j, k)*(one_over_dy*(velocity_x(i, j + 1, k) - velocity_x(i, j, k)) + one_over_dx*(velocity_y(i, j + 1, k) - velocity_y(i - 1, j + 1, k))) - viscosity_cd_x(i, j, k)*(one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)) + one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k))));
					T third_update = coef*one_over_dz*(viscosity_ct_x(i, j, k)*(one_over_dz*(velocity_x(i, j, k + 1) - velocity_x(i, j, k)) + one_over_dx*(velocity_z(i, j, k + 1) - velocity_z(i - 1, j, k + 1))) - viscosity_cb_x(i, j, k)*(one_over_dz*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) + one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
					velocity_x(i, j, k) += first_update + second_update + third_update;
				}
				
				multithreading.Sync(thread_id);
	
				// Define for coefficient function -- coefficients of y-velocity
				GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
				{
					levelset_y_half(i, j, k) = (T)0.5*(water_array(i, j - 1, k) + water_array(i, j, k));
					levelset_cu_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k) + water_array(i + 1, j - 1, k) + water_array(i, j, k) + water_array(i + 1, j, k));
					levelset_cl_y(i, j, k) = (T)0.25*(water_array(i - 1, j - 1, k) + water_array(i, j - 1, k) + water_array(i - 1, j, k) + water_array(i, j, k));
					levelset_ct_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k + 1) + water_array(i, j - 1, k) + water_array(i, j, k + 1) + water_array(i, j, k));
					levelset_cb_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k - 1) + water_array(i, j - 1, k) + water_array(i, j, k - 1) + water_array(i, j, k));
					levelset_u_y(i, j, k) = water_array(i, j, k);
					levelset_c_y(i, j, k) = water_array(i, j - 1, k);
				}
				
				multithreading.Sync(thread_id);
	
				DetermineDensityField(thread_id, levelset_y_half, epsilon_for_mollification, density_half_y);
				DetermineViscosityField(thread_id, levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
				DetermineViscosityField(thread_id, levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
				DetermineViscosityField(thread_id, levelset_u_y, epsilon_for_mollification, viscosity_u_y);
				DetermineViscosityField(thread_id, levelset_c_y, epsilon_for_mollification, viscosity_c_y);
				DetermineViscosityField(thread_id, levelset_ct_y, epsilon_for_mollification, viscosity_ct_y);
				DetermineViscosityField(thread_id, levelset_cb_y, epsilon_for_mollification, viscosity_cb_y);
				
				GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
				{
					T one_over_density_half_y = (T)1/density_half_y(i, j, k);
					T coef = dt*one_over_density_half_y;
					T first_update = coef*one_over_dx*(viscosity_cu_y(i, j, k)*(one_over_dy*(velocity_x(i + 1, j, k) - velocity_x(i + 1, j - 1, k)) + one_over_dx*(velocity_y(i + 1, j, k) - velocity_y(i, j, k))) - viscosity_cl_y(i, j, k)*(one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)) + one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k))));
					T second_update = (T)2*coef*one_over_dy*(viscosity_u_y(i, j, k)*one_over_dy*(velocity_y(i, j + 1, k) - velocity_y(i, j, k)) - viscosity_c_y(i, j, k)*one_over_dy*(velocity_y(i, j, k) - velocity_y(i, j - 1, k)));
					T third_update = coef*one_over_dz*(viscosity_ct_y(i, j, k)*(one_over_dz*(velocity_y(i, j, k + 1) - velocity_y(i, j, k)) + one_over_dy*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))) - viscosity_cb_y(i, j, k)*(one_over_dz*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) + one_over_dy*(velocity_z(i, j, k - 1) - velocity_z(i, j - 1, k - 1))));
					velocity_y(i, j, k) += first_update + second_update + third_update;
				}
				
				multithreading.Sync(thread_id);
	
				// Define for coefficient function -- coefficients of z-velocity
				
				GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
				{
					levelset_z_half(i, j, k) = (T)0.5*(water_array(i, j, k - 1) + water_array(i, j, k));
					levelset_ct_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i + 1, j, k - 1) + water_array(i, j, k) + water_array(i + 1, j, k));
					levelset_cb_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i - 1, j, k - 1) + water_array(i, j, k) + water_array(i - 1, j, k));
					levelset_cu_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j + 1, k - 1) + water_array(i, j, k) + water_array(i, j + 1, k));
					levelset_cd_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j - 1, k - 1) + water_array(i, j, k) + water_array(i, j - 1, k));
					levelset_t_z(i, j, k) = water_array(i, j, k);
					levelset_c_z(i, j, k) = water_array(i, j, k - 1);
				}
				
				multithreading.Sync(thread_id);

				DetermineDensityField(thread_id, levelset_z_half, epsilon_for_mollification, density_half_z);
				DetermineViscosityField(thread_id, levelset_ct_z, epsilon_for_mollification, viscosity_ct_z);
				DetermineViscosityField(thread_id, levelset_cb_z, epsilon_for_mollification, viscosity_cb_z);
				DetermineViscosityField(thread_id, levelset_cu_z, epsilon_for_mollification, viscosity_cu_z);
				DetermineViscosityField(thread_id, levelset_cd_z, epsilon_for_mollification, viscosity_cd_z);
				DetermineViscosityField(thread_id, levelset_t_z, epsilon_for_mollification, viscosity_t_z);
				DetermineViscosityField(thread_id, levelset_c_z, epsilon_for_mollification, viscosity_c_z);
				
				GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
				{
					T one_over_density_half_z = (T)1/density_half_z(i, j, k);
					T coef = dt*one_over_density_half_z;
					T first_update = coef*one_over_dx*(viscosity_ct_z(i, j, k)*(one_over_dz*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) + one_over_dx*(velocity_z(i + 1, j, k) - velocity_z(i, j, k))) - viscosity_cb_z(i, j, k)*(one_over_dz*(velocity_x(i - 1, j, k) - velocity_x(i - 1, j, k - 1)) + one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
					T second_update = coef*one_over_dy*(viscosity_cu_z(i, j, k)*(one_over_dz*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) + one_over_dy*(velocity_z(i, j + 1, k) - velocity_z(i, j, k))) - viscosity_cd_z(i, j, k)*(one_over_dz*(velocity_y(i, j - 1, k) - velocity_y(i, j - 1, k - 1)) + one_over_dy*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))));
					T third_update = (T)2*coef*one_over_dz*(viscosity_t_z(i, j, k)*one_over_dz*(velocity_z(i, j, k + 1) - velocity_z(i, j, k)) - viscosity_c_z(i, j, k)*one_over_dz*(velocity_z(i, j, k) - velocity_z(i, j, k - 1)));
					velocity_z(i, j, k) += first_update + second_update + third_update;
				}
			
				multithreading.Sync(thread_id);
			}

			if (oil_water_simulation)
			{
				if (semi_implicit_approach)
				{
					// Speed Up Variable
					//const ARRAY_3D<T>& water_array(water_levelset.arr);
					const ARRAY_3D<T>& water_array(scalar_field_ghost.array_for_this);
			
					SetupBoundaryConditionsForVelocity(thread_id);

					const ARRAY_3D<T>& velocity_x(velocity_field_mac_ghost_x.array_for_this);
					const ARRAY_3D<T>& velocity_y(velocity_field_mac_ghost_y.array_for_this);
					const ARRAY_3D<T>& velocity_z(velocity_field_mac_ghost_z.array_for_this);
					
					if (dimensionless_form)
					{
						// Define the coefficient function -- coefficients of x-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
						{
							levelset_x_half(i, j, k) = (T)0.5*(water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_cu_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j + 1, k) + water_array(i, j + 1, k));
							levelset_cd_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j - 1, k) + water_array(i, j - 1, k));
							levelset_ct_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k + 1) + water_array(i, j, k + 1));
							levelset_cb_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k - 1) + water_array(i, j, k - 1));
							levelset_r_x(i, j, k) = water_array(i, j, k);
							levelset_c_x(i, j, k) = water_array(i - 1, j, k);	
						}
				
						multithreading.Sync(thread_id);
						
						DetermineDensityField(thread_id, levelset_x_half, epsilon_for_mollification, density_half_x);
						DetermineViscosityField(thread_id, levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
						DetermineViscosityField(thread_id, levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
						DetermineViscosityField(thread_id, levelset_r_x, epsilon_for_mollification, viscosity_r_x);
						DetermineViscosityField(thread_id, levelset_c_x, epsilon_for_mollification, viscosity_c_x);
						DetermineViscosityField(thread_id, levelset_ct_x, epsilon_for_mollification, viscosity_ct_x);
						DetermineViscosityField(thread_id, levelset_cb_x, epsilon_for_mollification, viscosity_cb_x);

						GRID_ITERATION_3D(coef_1_x.partial_grids[thread_id])
						{
							coef_1_x(i, j, k) = (T)2*dt/(Re*density_half_x(i, j, k))*viscosity_r_x(i, j, k);	
							coef_2_x(i, j, k) = (T)2*dt/(Re*density_half_x(i, j, k))*viscosity_c_x(i, j, k);
							coef_3_x(i, j, k) = dt/(Re*density_half_x(i, j, k))*viscosity_cu_x(i, j, k);
							coef_4_x(i, j, k) = dt/(Re*density_half_x(i, j, k))*viscosity_cd_x(i, j, k);
							coef_5_x(i, j, k) = dt/(Re*density_half_x(i, j, k))*viscosity_ct_x(i, j, k);
							coef_6_x(i, j, k) = dt/(Re*density_half_x(i, j, k))*viscosity_cb_x(i, j, k);
						}

						multithreading.Sync(thread_id);

						if (is_vertical)
						{
							GRID_ITERATION_3D(explicit_term_x.partial_grids[thread_id])
							{
								explicit_term_x(i, j, k) = velocity_x(i, j, k) + dt/(Re*density_half_x(i, j, k))*(one_over_dy*(viscosity_cu_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j + 1, k) - velocity_field_mac_ghost_y_x(i - 1, j + 1, k)) - viscosity_cd_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j, k) - velocity_field_mac_ghost_y_x(i - 1, j, k)))) + dt/(Re*density_half_x(i, j, k))*(one_over_dz*(viscosity_ct_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_z_x(i, j, k + 1) - velocity_field_mac_ghost_z_x(i - 1, j, k + 1)) - viscosity_cb_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_z_x(i, j, k) - velocity_field_mac_ghost_z_x(i - 1, j, k))));
							}

							multithreading.Sync(thread_id); 
						}
						
						if (is_parallel)
						{
							GRID_ITERATION_3D(explicit_term_x.partial_grids[thread_id])
							{
								explicit_term_x(i, j, k) = velocity_x(i, j, k) + dt/(Re*density_half_x(i, j, k))*(one_over_dy*(viscosity_cu_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j + 1, k) - velocity_field_mac_ghost_y_x(i - 1, j + 1, k)) - viscosity_cd_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j, k) - velocity_field_mac_ghost_y_x(i - 1, j, k)))) + dt/(Re*density_half_x(i, j, k))*(one_over_dz*(viscosity_ct_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_z_x(i, j, k + 1) - velocity_field_mac_ghost_z_x(i - 1, j, k + 1)) - viscosity_cb_x(i, j, k)*one_over_dx*(velocity_field_mac_ghost_z_x(i, j, k) - velocity_field_mac_ghost_z_x(i - 1, j, k))));
							}

							multithreading.Sync(thread_id); 
						}

						//ofstream fout;
						//fout.open("explicit_term_x");
						//for (int i = explicit_term_x.i_start; i <= explicit_term_x.i_end; i++)
						//{
						//	for (int j = explicit_term_x.j_start; j <= explicit_term_x.j_end; j++)
						//	{
						//		fout << explicit_term_x.array_for_this(i, j, explicit_term_x.k_end/2) << " ";
						//	}
						//	fout << "\n";
						//}
						//fout.close();

						//fout.open("density_half_x");
						//for (int i = density_half_x.i_start; i <= density_half_x.i_end; i++)
						//{
						//	for (int j = density_half_x.j_start; j <= density_half_x.j_end; j++)
						//	{
						//		fout << density_half_x.array_for_this(i, j, density_half_x.k_end/2) << " ";
						//	}
						//	fout << "\n";
						//}
						//fout.close();

						//fout.open("viscosity_cu_x");
						//for (int i = viscosity_cu_x.i_start; i <= viscosity_cu_x.i_end; i++)
						//{
						//	for (int j = viscosity_cu_x.j_start; j <= viscosity_cu_x.j_end; j++)
						//	{
						//		fout << viscosity_cu_x.array_for_this(i, j, viscosity_cu_x.k_end/2) << " ";
						//	}
						//	fout << "\n";
						//}
						//fout.close();

						//fout.open("first_term_x_1");
						//for (int i = explicit_term_x.i_start; i <= explicit_term_x.i_end; i++)
						//{
						//	for (int j = explicit_term_x.j_start; j <= explicit_term_x.j_end; j++)
						//	{
						//		//fout << dt/(Re*density_half_x(i, j, explicit_term_x.k_end/2))*(one_over_dy*(viscosity_cu_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j + 1, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j + 1, explicit_term_x.k_end/2)) - viscosity_cd_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j, explicit_term_x.k_end/2)))) << " ";
						//		fout << dt/(Re*density_half_x(i, j, explicit_term_x.k_end/2))*(one_over_dy*(viscosity_cu_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j + 1, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j + 1, explicit_term_x.k_end/2)))) << " ";
						//		//fout << velocity_field_mac_ghost_y_x(i, j, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j, explicit_term_x.k_end/2) << " ";
						//	}
						//	fout << "\n";
						//}
						//fout.close();

						//fout.open("first_term_x_2");
						//for (int i = explicit_term_x.i_start; i <= explicit_term_x.i_end; i++)
						//{
						//	for (int j = explicit_term_x.j_start; j <= explicit_term_x.j_end; j++)
						//	{
						//		//fout << dt/(Re*density_half_x(i, j, explicit_term_x.k_end/2))*(one_over_dy*(viscosity_cu_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j + 1, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j + 1, explicit_term_x.k_end/2)) - viscosity_cd_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j, explicit_term_x.k_end/2)))) << " ";
						//		//fout << one_over_dy*(viscosity_cu_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j + 1, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j + 1, explicit_term_x.k_end/2))) << " ";
						//		fout << dt/(Re*density_half_x(i, j, explicit_term_x.k_end/2))*(one_over_dy*(viscosity_cd_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_y_x(i, j, explicit_term_x.k_end/2) - velocity_field_mac_ghost_y_x(i - 1, j, explicit_term_x.k_end/2)))) << " ";
						//	}
						//	fout << "\n";
						//}
						//fout.close();

						//fout.open("second_term_x");
						//for (int i = explicit_term_x.i_start; i <= explicit_term_x.i_end; i++)
						//{
						//	for (int j = explicit_term_x.j_start; j <= explicit_term_x.j_end; j++)
						//	{
						//		fout << dt/(Re*density_half_x(i, j, explicit_term_x.k_end/2))*(one_over_dz*(viscosity_ct_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_z_x(i, j, explicit_term_x.k_end/2 + 1) - velocity_field_mac_ghost_z_x(i - 1, j, explicit_term_x.k_end/2 + 1)) - viscosity_cb_x(i, j, explicit_term_x.k_end/2)*one_over_dx*(velocity_field_mac_ghost_z_x(i, j, explicit_term_x.k_end/2) - velocity_field_mac_ghost_z_x(i - 1, j, explicit_term_x.k_end/2)))) << " ";
						//	}
						//	fout << "\n";
						//}
						//fout.close();

						SetupBoundaryConditionsForVelocity(thread_id, water_velocity_field_mac_x, boundary_condition_field_x);

						if (is_vertical)
						{
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x, boundary_condition_field_x, explicit_term_x);
							poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_z, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x, boundary_condition_field_x, explicit_term_x);
							//poisson_solver.SolveForViscosityYPeriodic(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x, boundary_condition_field_x, explicit_term_x);
						}
						if (is_parallel)
						{
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_z, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x, boundary_condition_field_x, explicit_term_x, false);
							poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_y, velocity_field_mac_ghost_x_z, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x, boundary_condition_field_x, explicit_term_x, false);
						}
												
						// Define for coefficient function -- coefficients of y-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
						{
							levelset_y_half(i, j, k) = (T)0.5*(water_array(i, j - 1, k) + water_array(i, j, k));
							levelset_cu_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k) + water_array(i + 1, j - 1, k) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cl_y(i, j, k) = (T)0.25*(water_array(i - 1, j - 1, k) + water_array(i, j - 1, k) + water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_ct_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k + 1) + water_array(i, j - 1, k) + water_array(i, j, k + 1) + water_array(i, j, k));
							levelset_cb_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k - 1) + water_array(i, j - 1, k) + water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_u_y(i, j, k) = water_array(i, j, k);
							levelset_c_y(i, j, k) = water_array(i, j - 1, k);
						}
				
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_y_half, epsilon_for_mollification, density_half_y);
						DetermineViscosityField(thread_id, levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
						DetermineViscosityField(thread_id, levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
						DetermineViscosityField(thread_id, levelset_u_y, epsilon_for_mollification, viscosity_u_y);
						DetermineViscosityField(thread_id, levelset_c_y, epsilon_for_mollification, viscosity_c_y);
						DetermineViscosityField(thread_id, levelset_ct_y, epsilon_for_mollification, viscosity_ct_y);
						DetermineViscosityField(thread_id, levelset_cb_y, epsilon_for_mollification, viscosity_cb_y);
	
						GRID_ITERATION_3D(coef_1_y.partial_grids[thread_id])
						{
							coef_1_y(i, j, k) = dt/(Re*density_half_y(i, j, k))*viscosity_cu_y(i, j, k);
							coef_2_y(i, j, k) = dt/(Re*density_half_y(i, j, k))*viscosity_cl_y(i, j, k);
							coef_3_y(i, j, k) = (T)2*dt/(Re*density_half_y(i, j, k))*viscosity_u_y(i, j, k);
							coef_4_y(i, j, k) = (T)2*dt/(Re*density_half_y(i, j, k))*viscosity_c_y(i, j, k);
							coef_5_y(i, j, k) = dt/(Re*density_half_y(i, j, k))*viscosity_ct_y(i, j, k);
							coef_6_y(i, j, k) = dt/(Re*density_half_y(i, j, k))*viscosity_cb_y(i, j, k);
						}
						
						multithreading.Sync(thread_id);

						if (is_vertical)
						{
							GRID_ITERATION_3D(explicit_term_y.partial_grids[thread_id])
							{
								//explicit_term_y(i, j, k) = velocity_y(i, j, k) + dt/(Re*density_half_y(i, j, k))*(one_over_dx*(viscosity_cu_y(i, j, k)*one_over_dy*(velocity_x(i + 1, j, k) - velocity_x(i + 1, j - 1, k)) - viscosity_cl_y(i, j, k)*one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)))) + dt/(Re*density_half_y(i, j ,k))*(one_over_dz*(viscosity_ct_y(i, j, k)*one_over_dy*(velocity_z(i, j, k + 1) - velocity_z(i, j - 1, k + 1)) - viscosity_cb_y(i, j, k)*one_over_dx*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))));
								explicit_term_y(i, j, k) = velocity_y(i, j, k) + dt/(Re*density_half_y(i, j, k))*(one_over_dx*(viscosity_cu_y(i, j, k)*one_over_dy*(velocity_field_mac_ghost_x_x(i + 1, j, k) - velocity_field_mac_ghost_x_x(i + 1, j - 1, k)) - viscosity_cl_y(i, j, k)*one_over_dy*(velocity_field_mac_ghost_x_x(i, j, k) - velocity_field_mac_ghost_x_x(i, j - 1, k)))) + dt/(Re*density_half_y(i, j ,k))*(one_over_dz*(viscosity_ct_y(i, j, k)*one_over_dy*(velocity_field_mac_ghost_z_z(i, j, k + 1) - velocity_field_mac_ghost_z_z(i, j - 1, k + 1)) - viscosity_cb_y(i, j, k)*one_over_dx*(velocity_field_mac_ghost_z_z(i, j, k) - velocity_field_mac_ghost_z_z(i, j - 1, k))));
							}
							multithreading.Sync(thread_id); 
						}
						if (is_parallel)
						{
							GRID_ITERATION_3D(explicit_term_y.partial_grids[thread_id])
							{
								explicit_term_y(i, j, k) = velocity_y(i, j, k) + dt/(Re*density_half_y(i, j, k))*(one_over_dx*(viscosity_cu_y(i, j, k)*one_over_dy*(velocity_field_mac_ghost_x_y(i + 1, j, k) - velocity_field_mac_ghost_x_y(i + 1, j - 1, k)) - viscosity_cl_y(i, j, k)*one_over_dy*(velocity_field_mac_ghost_x_y(i, j, k) - velocity_field_mac_ghost_x_y(i, j - 1, k)))) + dt/(Re*density_half_y(i, j ,k))*(one_over_dz*(viscosity_ct_y(i, j, k)*one_over_dy*(velocity_field_mac_ghost_z_y(i, j, k + 1) - velocity_field_mac_ghost_z_y(i, j - 1, k + 1)) - viscosity_cb_y(i, j, k)*one_over_dx*(velocity_field_mac_ghost_z_y(i, j, k) - velocity_field_mac_ghost_z_y(i, j - 1, k))));
							}
							multithreading.Sync(thread_id); 
						}

						SetupBoundaryConditionsForVelocity(thread_id, water_velocity_field_mac_y, boundary_condition_field_y);
						
						if (is_vertical)
						{
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y); 
							poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_z, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y);
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y, boundary_phi_y);
							//poisson_solver.SolveForViscosityYPeriodic(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y);
						}
						if (is_parallel)
						{
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_x, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_z, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y, false);
							poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_y_y, velocity_field_mac_ghost_y_z, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y, false); 
						}
						
						// Define for coefficient function -- coefficients of z-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
						{
							levelset_z_half(i, j, k) = (T)0.5*(water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_ct_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i + 1, j, k - 1) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cb_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i - 1, j, k - 1) + water_array(i, j, k) + water_array(i - 1, j, k));
							levelset_cu_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j + 1, k - 1) + water_array(i, j, k) + water_array(i, j + 1, k));
							levelset_cd_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j - 1, k - 1) + water_array(i, j, k) + water_array(i, j - 1, k));
							levelset_t_z(i, j, k) = water_array(i, j, k);
							levelset_c_z(i, j, k) = water_array(i, j, k - 1);
						}
						
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_z_half, epsilon_for_mollification, density_half_z);
						DetermineViscosityField(thread_id, levelset_ct_z, epsilon_for_mollification, viscosity_ct_z);
						DetermineViscosityField(thread_id, levelset_cb_z, epsilon_for_mollification, viscosity_cb_z);
						DetermineViscosityField(thread_id, levelset_cu_z, epsilon_for_mollification, viscosity_cu_z);
						DetermineViscosityField(thread_id, levelset_cd_z, epsilon_for_mollification, viscosity_cd_z);
						DetermineViscosityField(thread_id, levelset_t_z, epsilon_for_mollification, viscosity_t_z);
						DetermineViscosityField(thread_id, levelset_c_z, epsilon_for_mollification, viscosity_c_z);
					
						GRID_ITERATION_3D(coef_1_z.partial_grids[thread_id])
						{
							coef_1_z(i, j, k) = dt/(Re*density_half_z(i, j, k))*viscosity_ct_z(i, j, k);
							coef_2_z(i, j, k) = dt/(Re*density_half_z(i, j, k))*viscosity_cb_z(i, j, k);
							coef_3_z(i, j, k) = dt/(Re*density_half_z(i, j, k))*viscosity_cu_z(i, j, k);
							coef_4_z(i, j, k) = dt/(Re*density_half_z(i, j, k))*viscosity_cd_z(i, j, k);
							coef_5_z(i, j, k) = (T)2*dt/(Re*density_half_z(i, j, k))*viscosity_t_z(i, j, k);
							coef_6_z(i, j, k) = (T)2*dt/(Re*density_half_z(i, j, k))*viscosity_c_z(i, j, k);
						}
	
						multithreading.Sync(thread_id);
	
						if (is_vertical)
						{
							GRID_ITERATION_3D(explicit_term_z.partial_grids[thread_id])
							{
								//explicit_term_z(i, j, k) = velocity_z(i, j, k) + dt/(Re*density_half_z(i, j, k))*(one_over_dx*(one_over_dz*(viscosity_ct_z(i, j, k)*(velocity_x(i, j, k) - velocity_x(i, j, k - 1))) - one_over_dz*(viscosity_cb_z(i, j, k)*(velocity_x(i - 1, j, k) - velocity_x(i - 1, j, k - 1))))) + dt/(Re*density_half_z(i, j, k))*(one_over_dy*(one_over_dz*(viscosity_cu_z(i, j, k)*(velocity_y(i, j, k) - velocity_y(i, j, k - 1))) - one_over_dx*(viscosity_cd_z(i, j, k)*(velocity_y(i, j - 1, k) - velocity_y(i, j - 1, k - 1)))));
								explicit_term_z(i, j, k) = velocity_z(i, j, k) + dt/(Re*density_half_z(i, j, k))*(one_over_dx*(one_over_dz*(viscosity_ct_z(i, j, k)*(velocity_field_mac_ghost_x_z(i + 1, j, k) - velocity_field_mac_ghost_x_z(i + 1, j, k - 1))) - one_over_dz*(viscosity_cb_z(i, j, k)*(velocity_field_mac_ghost_x_z(i, j, k) - velocity_field_mac_ghost_x_z(i, j, k - 1))))) + dt/(Re*density_half_z(i, j, k))*(one_over_dy*(one_over_dz*(viscosity_cu_z(i, j, k)*(velocity_field_mac_ghost_y_z(i, j + 1, k) - velocity_field_mac_ghost_y_z(i, j + 1, k - 1))) - one_over_dz*(viscosity_cd_z(i, j, k)*(velocity_field_mac_ghost_y_z(i, j, k) - velocity_field_mac_ghost_y_z(i, j, k - 1)))));
							}
							multithreading.Sync(thread_id); 
						}
						if (is_parallel)
						{
							GRID_ITERATION_3D(explicit_term_z.partial_grids[thread_id])
							{
								explicit_term_z(i, j, k) = velocity_z(i, j, k) + dt/(Re*density_half_z(i, j, k))*(one_over_dx*(one_over_dz*(viscosity_ct_z(i, j, k)*(velocity_field_mac_ghost_x_z(i + 1, j, k) - velocity_field_mac_ghost_x_z(i + 1, j, k - 1))) - one_over_dz*(viscosity_cb_z(i, j, k)*(velocity_field_mac_ghost_x_z(i, j, k) - velocity_field_mac_ghost_x_z(i, j, k - 1))))) + dt/(Re*density_half_z(i, j, k))*(one_over_dy*(one_over_dz*(viscosity_cu_z(i, j, k)*(velocity_field_mac_ghost_y_z(i, j + 1, k) - velocity_field_mac_ghost_y_z(i, j + 1, k - 1))) - one_over_dx*(viscosity_cd_z(i, j, k)*(velocity_field_mac_ghost_y_z(i, j, k) - velocity_field_mac_ghost_y_z(i, j, k - 1)))));
							}
							multithreading.Sync(thread_id); 
						}

						SetupBoundaryConditionsForVelocity(thread_id, water_velocity_field_mac_z, boundary_condition_field_z);
	
						if (is_vertical)
						{
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z, boundary_condition_field_z, explicit_term_z); 
							poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_x, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z, boundary_condition_field_z, explicit_term_z);
							//poisson_solver.SolveForViscosityYPeriodic(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z, boundary_condition_field_z, explicit_term_z); 
						}
						if (is_parallel)
						{
							//poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_x, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z, boundary_condition_field_z, explicit_term_z, false);
							poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z, velocity_field_mac_ghost_z_y, velocity_field_mac_ghost_z_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z, boundary_condition_field_z, explicit_term_z, false); 
						}
					}
					else
					{
						// Define the coefficient function -- coefficients of x-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
						{
							levelset_x_half(i, j, k) = (T)0.5*(water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_cu_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j + 1, k) + water_array(i, j + 1, k));
							levelset_cd_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j - 1, k) + water_array(i, j - 1, k));
							levelset_ct_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k + 1) + water_array(i, j, k + 1));
							levelset_cb_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k - 1) + water_array(i, j, k - 1));
							levelset_r_x(i, j, k) = water_array(i, j, k);
							levelset_c_x(i, j, k) = water_array(i - 1, j, k);
						}
				
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_x_half, epsilon_for_mollification, density_half_x);
						DetermineViscosityField(thread_id, levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
						DetermineViscosityField(thread_id, levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
						DetermineViscosityField(thread_id, levelset_r_x, epsilon_for_mollification, viscosity_r_x);
						DetermineViscosityField(thread_id, levelset_c_x, epsilon_for_mollification, viscosity_c_x);
						DetermineViscosityField(thread_id, levelset_ct_x, epsilon_for_mollification, viscosity_ct_x);
						DetermineViscosityField(thread_id, levelset_cb_x, epsilon_for_mollification, viscosity_cb_x);
	
						GRID_ITERATION_3D(coef_1_x.partial_grids[thread_id])
						{
							coef_1_x(i, j, k) = (T)2*dt/density_half_x(i, j, k)*viscosity_r_x(i, j, k);	
							coef_2_x(i, j, k) = (T)2*dt/density_half_x(i, j, k)*viscosity_c_x(i, j, k);
							coef_3_x(i, j, k) = dt/density_half_x(i, j, k)*viscosity_cu_x(i, j, k);
							coef_4_x(i, j, k) = dt/density_half_x(i, j, k)*viscosity_cd_x(i, j, k);
							coef_5_x(i, j, k) = dt/density_half_x(i, j, k)*viscosity_ct_x(i, j, k);
							coef_6_x(i, j, k) = dt/density_half_x(i, j, k)*viscosity_cb_x(i, j, k);
						}
	
						multithreading.Sync(thread_id);
	
						GRID_ITERATION_3D(explicit_term_x.partial_grids[thread_id])
						{
							explicit_term_x(i, j, k) = velocity_x(i, j, k) + dt/density_half_x(i, j, k)*(one_over_dy*(viscosity_cu_x(i, j, k)*one_over_dx*(velocity_y(i, j + 1, k) - velocity_y(i - 1, j + 1, k)) - viscosity_cd_x(i, j, k)*one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k)))) + dt/density_half_x(i, j, k)*(one_over_dz*(viscosity_ct_x(i, j, k)*one_over_dx*(velocity_z(i, j, k + 1) - velocity_z(i - 1, j, k + 1)) - viscosity_cb_x(i, j, k)*one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
						}
	
						multithreading.Sync(thread_id);
	
						SetupBoundaryConditionsForVelocity(thread_id, water_velocity_field_mac_x, boundary_condition_field_x);
	
						poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, coef_6_x, boundary_condition_field_x, explicit_term_x);
	
						// Define for coefficient function -- coefficients of y-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
						{
							levelset_y_half(i, j, k) = (T)0.5*(water_array(i, j - 1, k) + water_array(i, j, k));
							levelset_cu_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k) + water_array(i + 1, j - 1, k) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cl_y(i, j, k) = (T)0.25*(water_array(i - 1, j - 1, k) + water_array(i, j - 1, k) + water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_ct_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k + 1) + water_array(i, j - 1, k) + water_array(i, j, k + 1) + water_array(i, j, k));
							levelset_cb_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k - 1) + water_array(i, j - 1, k) + water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_u_y(i, j, k) = water_array(i, j, k);
							levelset_c_y(i, j, k) = water_array(i, j - 1, k);
						}
					
						multithreading.Sync(thread_id);
		
						DetermineDensityField(thread_id, levelset_y_half, epsilon_for_mollification, density_half_y);
						DetermineViscosityField(thread_id, levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
						DetermineViscosityField(thread_id, levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
						DetermineViscosityField(thread_id, levelset_u_y, epsilon_for_mollification, viscosity_u_y);
						DetermineViscosityField(thread_id, levelset_c_y, epsilon_for_mollification, viscosity_c_y);
						DetermineViscosityField(thread_id, levelset_ct_y, epsilon_for_mollification, viscosity_ct_y);
						DetermineViscosityField(thread_id, levelset_cb_y, epsilon_for_mollification, viscosity_cb_y);
						
						GRID_ITERATION_3D(coef_1_y.partial_grids[thread_id])
						{
							coef_1_y(i, j, k) = dt/density_half_y(i, j, k)*viscosity_cu_y(i, j, k);
							coef_2_y(i, j, k) = dt/density_half_y(i, j, k)*viscosity_cl_y(i, j, k);
							coef_3_y(i, j, k) = (T)2*dt/density_half_y(i, j, k)*viscosity_u_y(i, j, k);
							coef_4_y(i, j, k) = (T)2*dt/density_half_y(i, j, k)*viscosity_c_y(i, j, k);
							coef_5_y(i, j, k) = dt/density_half_y(i, j, k)*viscosity_ct_y(i, j, k);
							coef_6_y(i, j, k) = dt/density_half_y(i, j, k)*viscosity_cb_y(i, j, k);
						}
						
						multithreading.Sync(thread_id);
	
						GRID_ITERATION_3D(explicit_term_y.partial_grids[thread_id])
						{
							explicit_term_y(i, j, k) = velocity_y(i, j, k) + dt/density_half_y(i, j, k)*(one_over_dx*(viscosity_cu_y(i, j, k)*one_over_dy*(velocity_x(i + 1, j, k) - velocity_x(i + 1, j - 1, k)) - viscosity_cl_y(i, j, k)*one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)))) + dt/density_half_y(i, j ,k)*(one_over_dz*(viscosity_ct_y(i, j, k)*one_over_dx*(velocity_z(i, j, k + 1) - velocity_z(i - 1, j, k + 1)) - viscosity_cb_y(i, j, k)*one_over_dx*(velocity_z(i, j, k - 1) - velocity_z(i, j - 1, k - 1))));
						}
	
						multithreading.Sync(thread_id);
	
						SetupBoundaryConditionsForVelocity(thread_id, water_velocity_field_mac_y, boundary_condition_field_y);
						
						poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, coef_6_y, boundary_condition_field_y, explicit_term_y);

						// Define for coefficient function -- coefficients of z-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
						{
							levelset_z_half(i, j, k) = (T)0.5*(water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_ct_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i + 1, j, k - 1) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cb_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i - 1, j, k - 1) + water_array(i, j, k) + water_array(i - 1, j, k));
							levelset_cu_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j + 1, k - 1) + water_array(i, j, k) + water_array(i, j + 1, k));
							levelset_cd_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j - 1, k - 1) + water_array(i, j, k) + water_array(i, j - 1, k));
							levelset_t_z(i, j, k) = water_array(i, j, k);
							levelset_c_z(i, j, k) = water_array(i, j, k - 1);
						}
						
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_z_half, epsilon_for_mollification, density_half_z);
						DetermineViscosityField(thread_id, levelset_ct_z, epsilon_for_mollification, viscosity_ct_z);
						DetermineViscosityField(thread_id, levelset_cb_z, epsilon_for_mollification, viscosity_cb_z);
						DetermineViscosityField(thread_id, levelset_cu_z, epsilon_for_mollification, viscosity_cu_z);
						DetermineViscosityField(thread_id, levelset_cd_z, epsilon_for_mollification, viscosity_cd_z);
						DetermineViscosityField(thread_id, levelset_t_z, epsilon_for_mollification, viscosity_t_z);
						DetermineViscosityField(thread_id, levelset_c_z, epsilon_for_mollification, viscosity_c_z);
					
						GRID_ITERATION_3D(coef_1_z.partial_grids[thread_id])
						{
							coef_1_z(i, j, k) = dt/density_half_z(i, j, k)*viscosity_ct_z(i, j, k);
							coef_2_z(i, j, k) = dt/density_half_z(i, j, k)*viscosity_cb_z(i, j, k);
							coef_3_z(i, j, k) = dt/density_half_z(i, j, k)*viscosity_cu_z(i, j, k);
							coef_4_z(i, j, k) = dt/density_half_z(i, j, k)*viscosity_cd_z(i, j, k);
							coef_5_z(i, j, k) = (T)2*dt/density_half_z(i, j, k)*viscosity_t_z(i, j, k);
							coef_6_z(i, j, k) = (T)2*dt/density_half_z(i, j, k)*viscosity_c_z(i, j, k);
						}
	
						multithreading.Sync(thread_id);

						GRID_ITERATION_3D(explicit_term_z.partial_grids[thread_id])
						{
							explicit_term_z(i, j, k) = water_velocity_field_mac_z(i, j, k) + dt/density_half_z(i, j, k)*(one_over_dx*one_over_dz*(viscosity_ct_z(i, j, k)*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) - viscosity_cb_z(i, j, k)*(velocity_x(i - 1, j, k) - velocity_x(i - 1, j, k - 1)))) + dt/density_half_z(i, j, k)*(one_over_dy*one_over_dz*(viscosity_cu_z(i, j, k)*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) - viscosity_cd_z(i, j, k)*(velocity_y(i, j - 1, k) - velocity_y(i, j - 1, k - 1))));
						}
						
						multithreading.Sync(thread_id);
					
						SetupBoundaryConditionsForVelocity(thread_id, water_velocity_field_mac_z, boundary_condition_field_z);
	
						poisson_solver.SolveForViscosity(thread_id, water_velocity_field_mac_z, coef_1_z, coef_2_z, coef_3_z, coef_4_z, coef_5_z, coef_6_z, boundary_condition_field_z, explicit_term_z);
					}
				}
				else
				{
					// Speed Up Variable
					const ARRAY_3D<T>& water_array(water_levelset.arr);
			
					const ARRAY_3D<T>& velocity_x(water_velocity_field_mac_x.array_for_this);
					const ARRAY_3D<T>& velocity_y(water_velocity_field_mac_y.array_for_this);
					const ARRAY_3D<T>& velocity_z(water_velocity_field_mac_z.array_for_this);
					
					SetupBoundaryConditionsForVelocity(thread_id);

					if (dimensionless_form)
					{
						// Define the coefficient function -- coefficients of x-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
						{
							levelset_x_half(i, j, k) = (T)0.5*(water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_cu_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j + 1, k) + water_array(i, j + 1, k));
							levelset_cd_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j - 1, k) + water_array(i, j - 1, k));
							levelset_ct_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k + 1) + water_array(i, j, k + 1));
							levelset_cb_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k - 1) + water_array(i, j, k - 1));
							levelset_r_x(i, j, k) = water_array(i, j, k);
							levelset_c_x(i, j, k) = water_array(i - 1, j, k);
						}
				
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_x_half, epsilon_for_mollification, density_half_x);
						DetermineViscosityField(thread_id, levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
						DetermineViscosityField(thread_id, levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
						DetermineViscosityField(thread_id, levelset_r_x, epsilon_for_mollification, viscosity_r_x);
						DetermineViscosityField(thread_id, levelset_c_x, epsilon_for_mollification, viscosity_c_x);
						DetermineViscosityField(thread_id, levelset_ct_x, epsilon_for_mollification, viscosity_ct_x);
						DetermineViscosityField(thread_id, levelset_cb_x, epsilon_for_mollification, viscosity_cb_x);
						
						GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
						{
							if (water_velocity_field_mac_x.fixed(i, j, k) == true)
							{
								continue;
							}

							T one_over_density_half_x = (T)1/density_half_x(i, j, k);
							T coef = dt/Re*one_over_density_half_x;
							T first_update = (T)2*coef*one_over_dx*(viscosity_r_x(i, j, k)*one_over_dx*(velocity_x(i + 1, j, k) - velocity_x(i, j, k)) - viscosity_c_x(i, j, k)*one_over_dy*(velocity_x(i, j, k) - velocity_x(i - 1, j, k)));
							T second_update = coef*one_over_dy*(viscosity_cu_x(i, j, k)*(one_over_dy*(velocity_x(i, j + 1, k) - velocity_x(i, j, k)) + one_over_dx*(velocity_y(i, j + 1, k) - velocity_y(i - 1, j + 1, k))) - viscosity_cd_x(i, j, k)*(one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)) + one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k))));
							T third_update = coef*one_over_dz*(viscosity_ct_x(i, j, k)*(one_over_dz*(velocity_x(i, j, k + 1) - velocity_x(i, j, k)) + one_over_dx*(velocity_z(i, j, k + 1) - velocity_z(i - 1, j, k + 1))) - viscosity_cb_x(i, j, k)*(one_over_dz*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) + one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
							velocity_x(i, j, k) += first_update + second_update + third_update;
						}
						
						multithreading.Sync(thread_id);
			
						// Define for coefficient function -- coefficients of y-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
						{
							levelset_y_half(i, j, k) = (T)0.5*(water_array(i, j - 1, k) + water_array(i, j, k));
							levelset_cu_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k) + water_array(i + 1, j - 1, k) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cl_y(i, j, k) = (T)0.25*(water_array(i - 1, j - 1, k) + water_array(i, j - 1, k) + water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_ct_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k + 1) + water_array(i, j - 1, k) + water_array(i, j, k + 1) + water_array(i, j, k));
							levelset_cb_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k - 1) + water_array(i, j - 1, k) + water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_u_y(i, j, k) = water_array(i, j, k);
							levelset_c_y(i, j, k) = water_array(i, j - 1, k);
						}
						
						multithreading.Sync(thread_id);
			
						DetermineDensityField(thread_id, levelset_y_half, epsilon_for_mollification, density_half_y);
						DetermineViscosityField(thread_id, levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
						DetermineViscosityField(thread_id, levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
						DetermineViscosityField(thread_id, levelset_u_y, epsilon_for_mollification, viscosity_u_y);
						DetermineViscosityField(thread_id, levelset_c_y, epsilon_for_mollification, viscosity_c_y);
						DetermineViscosityField(thread_id, levelset_ct_y, epsilon_for_mollification, viscosity_ct_y);
						DetermineViscosityField(thread_id, levelset_cb_y, epsilon_for_mollification, viscosity_cb_y);
						
						GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
						{
							if (water_velocity_field_mac_y.fixed(i, j, k) == true)
							{
								continue;
							}

							T one_over_density_half_y = (T)1/density_half_y(i, j, k);
							T coef = dt/Re*one_over_density_half_y;
							T first_update = coef*one_over_dx*(viscosity_cu_y(i, j, k)*(one_over_dy*(velocity_x(i + 1, j, k) - velocity_x(i + 1, j - 1, k)) + one_over_dx*(velocity_y(i + 1, j, k) - velocity_y(i, j, k))) - viscosity_cl_y(i, j, k)*(one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)) + one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k))));
							T second_update = (T)2*coef*one_over_dy*(viscosity_u_y(i, j, k)*one_over_dy*(velocity_y(i, j + 1, k) - velocity_y(i, j, k)) - viscosity_c_y(i, j, k)*one_over_dy*(velocity_y(i, j, k) - velocity_y(i, j - 1, k)));
							T third_update = coef*one_over_dz*(viscosity_ct_y(i, j, k)*(one_over_dz*(velocity_y(i, j, k + 1) - velocity_y(i, j, k)) + one_over_dy*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))) - viscosity_cb_y(i, j, k)*(one_over_dz*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) + one_over_dy*(velocity_z(i, j, k - 1) - velocity_z(i, j - 1, k - 1))));
							velocity_y(i, j, k) += first_update + second_update + third_update;
						}
						
						multithreading.Sync(thread_id);
			
						// Define for coefficient function -- coefficients of z-velocity
						
						GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
						{
							levelset_z_half(i, j, k) = (T)0.5*(water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_ct_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i + 1, j, k - 1) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cb_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i - 1, j, k - 1) + water_array(i, j, k) + water_array(i - 1, j, k));
							levelset_cu_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j + 1, k - 1) + water_array(i, j, k) + water_array(i, j + 1, k));
							levelset_cd_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j - 1, k - 1) + water_array(i, j, k) + water_array(i, j - 1, k));
							levelset_t_z(i, j, k) = water_array(i, j, k);
							levelset_c_z(i, j, k) = water_array(i, j, k - 1);
						}
						
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_z_half, epsilon_for_mollification, density_half_z);
						DetermineViscosityField(thread_id, levelset_ct_z, epsilon_for_mollification, viscosity_ct_z);
						DetermineViscosityField(thread_id, levelset_cb_z, epsilon_for_mollification, viscosity_cb_z);
						DetermineViscosityField(thread_id, levelset_cu_z, epsilon_for_mollification, viscosity_cu_z);
						DetermineViscosityField(thread_id, levelset_cd_z, epsilon_for_mollification, viscosity_cd_z);
						DetermineViscosityField(thread_id, levelset_t_z, epsilon_for_mollification, viscosity_t_z);
						DetermineViscosityField(thread_id, levelset_c_z, epsilon_for_mollification, viscosity_c_z);
						
						GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
						{
							if (water_velocity_field_mac_z.fixed(i, j, k) == true)
							{
								continue;
							}

							T one_over_density_half_z = (T)1/density_half_z(i, j, k);
							T coef = dt/Re*one_over_density_half_z;
							T first_update = coef*one_over_dx*(viscosity_ct_z(i, j, k)*(one_over_dz*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) + one_over_dx*(velocity_z(i + 1, j, k) - velocity_z(i, j, k))) - viscosity_cb_z(i, j, k)*(one_over_dz*(velocity_x(i - 1, j, k) - velocity_x(i - 1, j, k - 1)) + one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
							T second_update = coef*one_over_dy*(viscosity_cu_z(i, j, k)*(one_over_dz*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) + one_over_dy*(velocity_z(i, j + 1, k) - velocity_z(i, j, k))) - viscosity_cd_z(i, j, k)*(one_over_dz*(velocity_y(i, j - 1, k) - velocity_y(i, j - 1, k - 1)) + one_over_dy*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))));
							T third_update = (T)2*coef*one_over_dz*(viscosity_t_z(i, j, k)*one_over_dz*(velocity_z(i, j, k + 1) - velocity_z(i, j, k)) - viscosity_c_z(i, j, k)*one_over_dz*(velocity_z(i, j, k) - velocity_z(i, j, k - 1)));
							velocity_z(i, j, k) += first_update + second_update + third_update;
						}
				
						multithreading.Sync(thread_id);
					}
					else
					{
						// Define the coefficient function -- coefficients of x-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
						{
							levelset_x_half(i, j, k) = (T)0.5*(water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_cu_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j + 1, k) + water_array(i, j + 1, k));
							levelset_cd_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j - 1, k) + water_array(i, j - 1, k));
							levelset_ct_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k + 1) + water_array(i, j, k + 1));
							levelset_cb_x(i, j, k) = (T)0.25*(water_array(i - 1, j, k) + water_array(i, j, k) + water_array(i - 1, j, k - 1) + water_array(i, j, k - 1));
							levelset_r_x(i, j, k) = water_array(i, j, k);
							levelset_c_x(i, j, k) = water_array(i - 1, j, k);
						}
					
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_x_half, epsilon_for_mollification, density_half_x);
						DetermineViscosityField(thread_id, levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
						DetermineViscosityField(thread_id, levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
						DetermineViscosityField(thread_id, levelset_r_x, epsilon_for_mollification, viscosity_r_x);
						DetermineViscosityField(thread_id, levelset_c_x, epsilon_for_mollification, viscosity_c_x);
						DetermineViscosityField(thread_id, levelset_ct_x, epsilon_for_mollification, viscosity_ct_x);
						DetermineViscosityField(thread_id, levelset_cb_x, epsilon_for_mollification, viscosity_cb_x);
					
						GRID_ITERATION_3D(water_velocity_field_mac_x.partial_grids[thread_id])
						{
							T one_over_density_half_x = (T)1/density_half_x(i, j, k);
							T coef = dt*one_over_density_half_x;
							T first_update = (T)2*coef*one_over_dx*(viscosity_r_x(i, j, k)*one_over_dx*(velocity_x(i + 1, j, k) - velocity_x(i, j, k)) - viscosity_c_x(i, j, k)*one_over_dy*(velocity_x(i, j, k) - velocity_x(i - 1, j, k)));
							T second_update = coef*one_over_dy*(viscosity_cu_x(i, j, k)*(one_over_dy*(velocity_x(i, j + 1, k) - velocity_x(i, j, k)) + one_over_dx*(velocity_y(i, j + 1, k) - velocity_y(i - 1, j + 1, k))) - viscosity_cd_x(i, j, k)*(one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)) + one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k))));
							T third_update = coef*one_over_dz*(viscosity_ct_x(i, j, k)*(one_over_dz*(velocity_x(i, j, k + 1) - velocity_x(i, j, k)) + one_over_dx*(velocity_z(i, j, k + 1) - velocity_z(i - 1, j, k + 1))) - viscosity_cb_x(i, j, k)*(one_over_dz*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) + one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
							velocity_x(i, j, k) += first_update + second_update + third_update;
						}
					
						multithreading.Sync(thread_id);
		
						// Define for coefficient function -- coefficients of y-velocity
						GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
						{
							levelset_y_half(i, j, k) = (T)0.5*(water_array(i, j - 1, k) + water_array(i, j, k));
							levelset_cu_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k) + water_array(i + 1, j - 1, k) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cl_y(i, j, k) = (T)0.25*(water_array(i - 1, j - 1, k) + water_array(i, j - 1, k) + water_array(i - 1, j, k) + water_array(i, j, k));
							levelset_ct_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k + 1) + water_array(i, j - 1, k) + water_array(i, j, k + 1) + water_array(i, j, k));
							levelset_cb_y(i, j, k) = (T)0.25*(water_array(i, j - 1, k - 1) + water_array(i, j - 1, k) + water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_u_y(i, j, k) = water_array(i, j, k);
							levelset_c_y(i, j, k) = water_array(i, j - 1, k);
						}
					
						multithreading.Sync(thread_id);
		
						DetermineDensityField(thread_id, levelset_y_half, epsilon_for_mollification, density_half_y);
						DetermineViscosityField(thread_id, levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
						DetermineViscosityField(thread_id, levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
						DetermineViscosityField(thread_id, levelset_u_y, epsilon_for_mollification, viscosity_u_y);
						DetermineViscosityField(thread_id, levelset_c_y, epsilon_for_mollification, viscosity_c_y);
						DetermineViscosityField(thread_id, levelset_ct_y, epsilon_for_mollification, viscosity_ct_y);
						DetermineViscosityField(thread_id, levelset_cb_y, epsilon_for_mollification, viscosity_cb_y);	
					
						GRID_ITERATION_3D(water_velocity_field_mac_y.partial_grids[thread_id])
						{
							T one_over_density_half_y = (T)1/density_half_y(i, j, k);
							T coef = dt*one_over_density_half_y;
							T first_update = coef*one_over_dx*(viscosity_cu_y(i, j, k)*(one_over_dy*(velocity_x(i + 1, j, k) - velocity_x(i + 1, j - 1, k)) + one_over_dx*(velocity_y(i + 1, j, k) - velocity_y(i, j, k))) - viscosity_cl_y(i, j, k)*(one_over_dy*(velocity_x(i, j, k) - velocity_x(i, j - 1, k)) + one_over_dx*(velocity_y(i, j, k) - velocity_y(i - 1, j, k))));
							T second_update = (T)2*coef*one_over_dy*(viscosity_u_y(i, j, k)*one_over_dy*(velocity_y(i, j + 1, k) - velocity_y(i, j, k)) - viscosity_c_y(i, j, k)*one_over_dy*(velocity_y(i, j, k) - velocity_y(i, j - 1, k)));
							T third_update = coef*one_over_dz*(viscosity_ct_y(i, j, k)*(one_over_dz*(velocity_y(i, j, k + 1) - velocity_y(i, j, k)) + one_over_dy*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))) - viscosity_cb_y(i, j, k)*(one_over_dz*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) + one_over_dy*(velocity_z(i, j, k - 1) - velocity_z(i, j - 1, k - 1))));
							velocity_y(i, j, k) += first_update + second_update + third_update;
						}
					
						multithreading.Sync(thread_id);
		
						// Define for coefficient function -- coefficients of z-velocity
					
						GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
						{
							levelset_z_half(i, j, k) = (T)0.5*(water_array(i, j, k - 1) + water_array(i, j, k));
							levelset_ct_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i + 1, j, k - 1) + water_array(i, j, k) + water_array(i + 1, j, k));
							levelset_cb_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i - 1, j, k - 1) + water_array(i, j, k) + water_array(i - 1, j, k));
							levelset_cu_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j + 1, k - 1) + water_array(i, j, k) + water_array(i, j + 1, k));
							levelset_cd_z(i, j, k) = (T)0.25*(water_array(i, j, k - 1) + water_array(i, j - 1, k - 1) + water_array(i, j, k) + water_array(i, j - 1, k));
							levelset_t_z(i, j, k) = water_array(i, j, k);
							levelset_c_z(i, j, k) = water_array(i, j, k - 1);
						}
						
						multithreading.Sync(thread_id);
	
						DetermineDensityField(thread_id, levelset_z_half, epsilon_for_mollification, density_half_z);
						DetermineViscosityField(thread_id, levelset_ct_z, epsilon_for_mollification, viscosity_ct_z);
						DetermineViscosityField(thread_id, levelset_cb_z, epsilon_for_mollification, viscosity_cb_z);
						DetermineViscosityField(thread_id, levelset_cu_z, epsilon_for_mollification, viscosity_cu_z);
						DetermineViscosityField(thread_id, levelset_cd_z, epsilon_for_mollification, viscosity_cd_z);
						DetermineViscosityField(thread_id, levelset_t_z, epsilon_for_mollification, viscosity_t_z);
						DetermineViscosityField(thread_id, levelset_c_z, epsilon_for_mollification, viscosity_c_z);
						
						GRID_ITERATION_3D(water_velocity_field_mac_z.partial_grids[thread_id])
						{
							T one_over_density_half_z = (T)1/density_half_z(i, j, k);
							T coef = dt*one_over_density_half_z;
							T first_update = coef*one_over_dx*(viscosity_ct_z(i, j, k)*(one_over_dz*(velocity_x(i, j, k) - velocity_x(i, j, k - 1)) + one_over_dx*(velocity_z(i + 1, j, k) - velocity_z(i, j, k))) - viscosity_cb_z(i, j, k)*(one_over_dz*(velocity_x(i - 1, j, k) - velocity_x(i - 1, j, k - 1)) + one_over_dx*(velocity_z(i, j, k) - velocity_z(i - 1, j, k))));
							T second_update = coef*one_over_dy*(viscosity_cu_z(i, j, k)*(one_over_dz*(velocity_y(i, j, k) - velocity_y(i, j, k - 1)) + one_over_dy*(velocity_z(i, j + 1, k) - velocity_z(i, j, k))) - viscosity_cd_z(i, j, k)*(one_over_dz*(velocity_y(i, j - 1, k) - velocity_y(i, j - 1, k - 1)) + one_over_dy*(velocity_z(i, j, k) - velocity_z(i, j - 1, k))));
							T third_update = (T)2*coef*one_over_dz*(viscosity_t_z(i, j, k)*one_over_dz*(velocity_z(i, j, k + 1) - velocity_z(i, j, k)) - viscosity_c_z(i, j, k)*one_over_dz*(velocity_z(i, j, k) - velocity_z(i, j, k - 1)));
							velocity_z(i, j, k) += first_update + second_update + third_update;
						}
					
						multithreading.Sync(thread_id);
					}
						
				}
			}
		}
	}

public: // Member Function
	void HeavisideFunction(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& heaviside)
	{
		T one_over_epsilon = (T)1/epsilon, one_over_pi = (T)1/PI;

		GRID_ITERATION_3D(phi.partial_grids[thread_id])
		{
			if (phi(i, j, k) < - epsilon)
			{
				heaviside(i, j, k) = (T)0;
			}
			else if (phi(i, j, k) > epsilon)
			{
				heaviside(i, j, k) = (T)1;
			}
			else
			{
				//heaviside(i, j, k) = (T)0.5 + phi(i, j, k)*one_over_epsilon*(T)0.5 + (T)0.5*one_over_pi*sin(PI*phi(i, j, k)*one_over_epsilon);
				heaviside(i, j, k) = (T)0.5*((T)1 + phi(i, j, k)*one_over_epsilon + one_over_pi*sin(PI*phi(i, j, k)*one_over_epsilon));
			}
		}
		
		multithreading.Sync(thread_id);
	}
	
	void DetermineViscosityField(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& viscosity)
	{
		FIELD_STRUCTURE_3D<T> heaviside_phi;
		heaviside_phi.Initialize(&multithreading, viscosity.grid, 2);

		HeavisideFunction(thread_id, phi, epsilon, heaviside_phi);

		multithreading.Sync(thread_id);

		if (air_water_simulation)
		{
            GRID_ITERATION_3D(phi.partial_grids[thread_id])
		    {
		    	if (air_bubble_rising)
		    	{
		    		viscosity(i, j, k) = air_viscosity + (water_viscosity - air_viscosity)*heaviside_phi(i, j, k);
		    	}
		    	if (water_drop)
		    	{
		    		viscosity(i, j, k) = water_viscosity + (air_viscosity - water_viscosity)*heaviside_phi(i, j, k);
		    	}
		    }
		}
		
		if (oil_water_simulation)
		{
            GRID_ITERATION_3D(phi.partial_grids[thread_id])
		    {
				if (dimensionless_form)
				{
					//viscosity(i, j, k) = oil_viscosity/water_viscosity + ((T)1 - oil_viscosity/water_viscosity)*heaviside_phi(i, j, k);
					viscosity(i, j, k) = (T)1 + (m - (T)1)*heaviside_phi(i, j, k);
				}
				else
				{
					viscosity(i, j, k) = oil_viscosity + (water_viscosity - oil_viscosity)*heaviside_phi(i, j, k);
				}
			}
		}

		viscosity.FillGhostCellsFrom(thread_id, viscosity.array_for_this, true);

		multithreading.Sync(thread_id);
	}

	void DetermineDensityField(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& density)
	{
		FIELD_STRUCTURE_3D<T> heaviside_phi;
		heaviside_phi.Initialize(&multithreading, density.grid, 2);

		HeavisideFunction(thread_id, phi, epsilon, heaviside_phi);
		
		multithreading.Sync(thread_id);
		
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
					density(i, j, k) = (T)1 + (eta - (T)1)*heaviside_phi(i, j, k);
				}
				else
				{
					density(i, j, k) = oil_density + (water_density - oil_density)*heaviside_phi(i, j, k);
				}
				
			}
		}

		density.FillGhostCellsFrom(thread_id, density.array_for_this, true);

		multithreading.Sync(thread_id);
	}

public: // Member Functions
	void SetupBoundaryConditionsForVelocity(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity_input, FIELD_STRUCTURE_3D<int>& bc_input)
	{
		ARRAY_3D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_3D& grid(bc_input.grid);

		if (is_vertical)
		{
			GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
			{
				T x_coor = boundary_phi_x.x_min + i*boundary_phi_x.dx, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
				boundary_phi_x(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
			{
				T x_coor = boundary_phi_y.x_min + i*boundary_phi_y.dx, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
				boundary_phi_y(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
			{
				T x_coor = boundary_phi_z.x_min + i*boundary_phi_z.dx, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
				boundary_phi_z(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
			}
			multithreading.Sync(thread_id); 
		}
		if (is_parallel)
		{
			GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
			{
				T y_coor = boundary_phi_x.y_min + j*boundary_phi_x.dy, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
				boundary_phi_x(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
			{
				T y_coor = boundary_phi_y.y_min + j*boundary_phi_y.dy, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
				boundary_phi_y(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
			}
			multithreading.Sync(thread_id);

			GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
			{
				T y_coor = boundary_phi_z.y_min + j*boundary_phi_z.dy, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
				boundary_phi_z(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
			}
			multithreading.Sync(thread_id); 
		}

		GRID_ITERATION_3D(bc_input.partial_grids_ghost[thread_id])
		{
			assert(velocity_input.i_end == water_levelset.arr.i_end && velocity_input.i_end == bc_array.i_end);
				
			//const T x_coor = grid.x_min + i*grid.dx, y_coor = grid.y_min + j*grid.dy, z_coor = grid.z_min + k*grid.dz;
			if (regular_boundary)
			{
				if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end || k < grid.k_start || k > grid.k_end)
				{
					bc_array(i, j, k) = BC_DIR;
					velocity_input(i, j, k) = 0;
				}
				else
				{
					bc_array(i, j, k) = BC_FULL;
				}
			}
			
			if (cylindrical_boundary)
			{
				if (is_vertical)
				{
					if (velocity_input.is_x_component)
					{
						if (boundary_phi_x(i, j, k) > 0)
						{
							bc_array(i, j, k) = BC_DIR;
						}
						else if (j < boundary_phi_x.grid.j_start)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else if (j > boundary_phi_x.grid.j_end)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					}
					if (velocity_input.is_y_component)
					{
						if (boundary_phi_y(i, j, k) > 0)
						{
							bc_array(i, j, k) = BC_DIR;
						}
						else if (j < boundary_phi_y.grid.j_start)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else if (j >= boundary_phi_y.grid.j_end)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					}
					if (velocity_input.is_z_component)
					{
						if (boundary_phi_z(i, j, k) > 0)
						{
							bc_array(i, j, k) = BC_DIR;
						}
						else if (j < boundary_phi_z.grid.j_start)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else if (j > boundary_phi_z.grid.j_end)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					} 
				}
				if (is_parallel)
				{
					if (velocity_input.is_x_component)
					{
						if (boundary_phi_x(i, j, k) > 0)
						{
							bc_array(i, j, k) = BC_DIR;
						}
						else if (i < boundary_phi_x.grid.i_start)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else if (i >= boundary_phi_x.grid.i_end)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					}
					if (velocity_input.is_y_component)
					{
						if (boundary_phi_y(i, j, k) > 0)
						{
							bc_array(i, j, k) = BC_DIR;
						}
						else if (i < boundary_phi_y.grid.i_start)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else if (i > boundary_phi_y.grid.i_end)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					}
					if (velocity_input.is_z_component)
					{
						if (boundary_phi_z(i, j, k) > 0)
						{
							bc_array(i, j, k) = BC_DIR;
						}
						else if (i < boundary_phi_z.grid.i_start)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else if (i > boundary_phi_z.grid.i_end)
						{
							bc_array(i, j, k) = BC_PER;
						}
						else
						{
							bc_array(i, j, k) = BC_FULL;
						}
					} 
				}
			}
		}
		
		multithreading.Sync(thread_id);
	}

	void SetupBoundaryConditionsForVelocity(const int& thread_id)
	{
		if (oil_water_simulation)
		{
			if (is_vertical)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInYDirection(thread_id, water_levelset.arr, true);

				GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_x.x_min + i*boundary_phi_x.dx, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
					boundary_phi_x(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_y.x_min + i*boundary_phi_y.dx, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
					boundary_phi_y(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					T x_coor = boundary_phi_z.x_min + i*boundary_phi_z.dx, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
					boundary_phi_z(i, j, k) = sqrt(POW2(x_coor) + POW2(z_coor)) - a;
				}
				multithreading.Sync(thread_id); 
			}

			if (is_parallel)
			{
				scalar_field_ghost.FillGhostCellsPeriodicInXDirection(thread_id, water_levelset.arr, true);

				GRID_ITERATION_3D(boundary_phi_x.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_x.y_min + j*boundary_phi_x.dy, z_coor = boundary_phi_x.z_min + k*boundary_phi_x.dz;
					boundary_phi_x(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_y.y_min + j*boundary_phi_y.dy, z_coor = boundary_phi_y.z_min + k*boundary_phi_y.dz;
					boundary_phi_y(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					T y_coor = boundary_phi_z.y_min + j*boundary_phi_z.dy, z_coor = boundary_phi_z.z_min + k*boundary_phi_z.dz;
					boundary_phi_z(i, j, k) = sqrt(POW2(y_coor) + POW2(z_coor)) - a;
				}
				multithreading.Sync(thread_id); 
			}

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
					if (boundary_phi_x(i, j, k) > 0)
					{	
						if ((boundary_phi_x(i - 1, j, k) > 0) && (boundary_phi_x(i + 1, j, k) > 0) && (boundary_phi_x(i, j, k - 1) > 0) && (boundary_phi_x(i, j, k + 1) > 0))
						{
							
							velocity_field_mac_ghost_x_x(i, j, k) = (T)0;
							velocity_field_mac_ghost_x_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_x(i - 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_x(i - 1, j, k))/(abs(boundary_phi_x(i - 1, j, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_x(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_x(i - 1, j, k);
						}
						else if (boundary_phi_x(i + 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_x(i + 1, j, k))/(abs(boundary_phi_x(i + 1, j, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_x(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_x(i + 1, j, k);								
						}
						else if (boundary_phi_x(i, j, k - 1) < 0)
						{
							T theta = abs(boundary_phi_x(i, j, k - 1))/(abs(boundary_phi_x(i, j, k - 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_x(i, j, k - 1);
						}
						else if (boundary_phi_x(i, j, k + 1) < 0)
						{
							T theta = abs(boundary_phi_x(i, j, k + 1))/(abs(boundary_phi_x(i, j, k + 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_x(i, j, k + 1);
						}
						else
						{
							velocity_field_mac_ghost_x_x(i, j, k) = (T)0;
							velocity_field_mac_ghost_x_z(i, j, k) = (T)0;
						}
					}
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					if (boundary_phi_z(i, j, k) > 0)
					{
						if ((boundary_phi_z(i - 1, j, k) > 0) && (boundary_phi_z(i + 1, j, k) > 0) && (boundary_phi_z(i, j, k - 1) > 0) && (boundary_phi_z(i, j, k + 1) > 0))
						{
							velocity_field_mac_ghost_z_x(i, j, k) = (T)0;
							velocity_field_mac_ghost_z_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_z(i - 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_z(i - 1, j, k))/(abs(boundary_phi_z(i - 1, j, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_x(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_z(i - 1, j, k);
						}
						else if (boundary_phi_z(i + 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_z(i + 1, j, k))/(abs(boundary_phi_z(i + 1, j, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_x(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_z(i + 1, j, k);
						}
						else if (boundary_phi_z(i, j, k - 1) < 0)
						{
							T theta = abs(boundary_phi_z(i, j, k - 1))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_z(i, j, k - 1);
						}
						else if (boundary_phi_z(i, j, k + 1) < 0)
						{
							T theta = abs(boundary_phi_z(i, j, k + 1))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k) = (theta - (T)1)/theta*velocity_field_mac_ghost_z(i, j, k + 1);
						}
						else
						{
							velocity_field_mac_ghost_z_x(i, j, k) = (T)0;
							velocity_field_mac_ghost_z_z(i, j, k) = (T)0;
						}
					}	
				}
				multithreading.Sync(thread_id);
			
				GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					if (boundary_phi_y(i, j, k) > 0)
					{
						if ((boundary_phi_y(i - 1, j, k) > 0) && (boundary_phi_y(i + 1, j, k) > 0) && (boundary_phi_y(i, j, k - 1) > 0) && (boundary_phi_y(i, j, k + 1) > 0))
						{
							velocity_field_mac_ghost_y_x(i, j, k) = (T)0;
							velocity_field_mac_ghost_y_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_y(i - 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_y(i - 1, j, k))/(abs(boundary_phi_y(i - 1, j, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_x(i, j, k) = (theta - (T)1)/theta*water_velocity_field_mac_y(i - 1, j, k);
						}
						else if (boundary_phi_y(i + 1, j, k) < 0)
						{
							T theta = abs(boundary_phi_y(i + 1, j, k))/(abs(boundary_phi_y(i + 1, j, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_x(i, j, k) = (theta - (T)1)/theta*water_velocity_field_mac_y(i + 1, j, k);
						}
						else if (boundary_phi_y(i, j, k - 1) < 0)
						{
							T theta = abs(boundary_phi_y(i, j, k - 1))/(abs(boundary_phi_y(i, j, k - 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k) = (theta - (T)1)/theta*water_velocity_field_mac_y(i, j, k - 1);
						}
						else if (boundary_phi_y(i, j, k + 1) < 0)
						{
							T theta = abs(boundary_phi_y(i, j, k + 1))/(abs(boundary_phi_y(i, j, k + 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k) = (theta - (T)1)/theta*water_velocity_field_mac_y(i, j, k + 1);
						}
						else
						{
							velocity_field_mac_ghost_y_x(i, j, k) = (T)0;
							velocity_field_mac_ghost_y_z(i, j, k) = (T)0;
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
					if (boundary_phi_x(i, j, k) > 0)
					{	
						if ((boundary_phi_x(i, j - 1, k) > 0) && (boundary_phi_x(i, j + 1, k) > 0) && (boundary_phi_x(i, j, k - 1) > 0) && (boundary_phi_x(i, j, k + 1) > 0))
						{
							velocity_field_mac_ghost_x_y(i, j, k) = (T)0;
							velocity_field_mac_ghost_x_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_x(i, j - 1, k) < -0)
						{
							T theta = abs(boundary_phi_x(i, j - 1, k))/(abs(boundary_phi_x(i, j - 1, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j - 1, k);
						}
						else if (boundary_phi_x(i, j + 1, k) < -0)
						{
							T theta = abs(boundary_phi_x(i, j + 1, k))/(abs(boundary_phi_x(i, j + 1, k)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j + 1, k);								
						}
						else if (boundary_phi_x(i, j, k - 1) < -0)
						{
							T theta = abs(boundary_phi_x(i, j, k - 1))/(abs(boundary_phi_x(i, j, k - 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k - 1);
						}
						else if (boundary_phi_x(i, j, k + 1) < -0)
						{
							T theta = abs(boundary_phi_x(i, j, k + 1))/(abs(boundary_phi_x(i, j, k + 1)) + abs(boundary_phi_x(i, j, k)));
							velocity_field_mac_ghost_x_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_x(i, j, k + 1);
						}
						else
						{
							velocity_field_mac_ghost_x_y(i, j, k) = (T)0;
							velocity_field_mac_ghost_x_z(i, j, k) = (T)0;
						}
					}
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_3D(boundary_phi_z.partial_grids[thread_id])
				{
					if (boundary_phi_z(i, j, k) > 0)
					{
						if ((boundary_phi_z(i, j - 1, k) > 0) && (boundary_phi_z(i, j + 1, k) > 0) && (boundary_phi_z(i, j, k - 1) > 0) && (boundary_phi_z(i, j, k + 1) > 0))
						{
							velocity_field_mac_ghost_z_y(i, j, k) = (T)0;
							velocity_field_mac_ghost_z_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_z(i, j - 1, k) < -0)
						{
							T theta = abs(boundary_phi_z(i, j - 1, k))/(abs(boundary_phi_z(i, j - 1, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j - 1, k);
						}
						else if (boundary_phi_z(i, j + 1, k) < -0)
						{
							T theta = abs(boundary_phi_z(i, j + 1, k))/(abs(boundary_phi_z(i, j + 1, k)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j + 1, k);
						}
						else if (boundary_phi_z(i, j, k - 1) < -0)
						{
							T theta = abs(boundary_phi_z(i, j, k - 1))/(abs(boundary_phi_z(i, j, k - 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k - 1);
						}
						else if (boundary_phi_z(i, j, k + 1) < -0)
						{
							T theta = abs(boundary_phi_z(i, j, k + 1))/(abs(boundary_phi_z(i, j, k + 1)) + abs(boundary_phi_z(i, j, k)));
							velocity_field_mac_ghost_z_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_z(i, j, k + 1);
						}
						else
						{
							velocity_field_mac_ghost_z_y(i, j, k) = (T)0;
							velocity_field_mac_ghost_z_z(i, j, k) = (T)0;
						}
					}	
				}
				multithreading.Sync(thread_id);
			
				GRID_ITERATION_3D(boundary_phi_y.partial_grids[thread_id])
				{
					if (boundary_phi_y(i, j, k) > 0)
					{
						if ((boundary_phi_y(i, j - 1, k) > 0) && (boundary_phi_y(i, j + 1, k) > 0) && (boundary_phi_y(i, j, k - 1) > 0) && (boundary_phi_y(i, j, k + 1) > 0))
						{
							velocity_field_mac_ghost_y_y(i, j, k) = (T)0;
							velocity_field_mac_ghost_y_z(i, j, k) = (T)0;
						}
						else if (boundary_phi_y(i, j - 1, k) < 0)
						{
							T theta = abs(boundary_phi_y(i, j - 1, k))/(abs(boundary_phi_y(i, j - 1, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j - 1, k);
						}
						else if (boundary_phi_y(i, j + 1, k) < 0)
						{
							T theta = abs(boundary_phi_y(i, j + 1, k))/(abs(boundary_phi_y(i, j + 1, k)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_y(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j + 1, k);
						}
						else if (boundary_phi_y(i, j, k - 1) < 0)
						{
							T theta = abs(boundary_phi_y(i, j, k - 1))/(abs(boundary_phi_y(i, j, k - 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k - 1);
						}
						else if (boundary_phi_y(i, j, k + 1) < 0)
						{
							T theta = abs(boundary_phi_y(i, j, k + 1))/(abs(boundary_phi_y(i, j, k + 1)) + abs(boundary_phi_y(i, j, k)));
							velocity_field_mac_ghost_y_z(i, j, k) = (theta - 1)/theta*velocity_field_mac_ghost_y(i, j, k + 1);
						}
						else
						{
							velocity_field_mac_ghost_y_y(i, j, k) = (T)0;
							velocity_field_mac_ghost_y_z(i, j, k) = (T)0;
						}
					}	
				}
				multithreading.Sync(thread_id);
			}
		}
	}
};