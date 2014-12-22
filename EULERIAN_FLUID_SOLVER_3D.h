#pragma once

#include "COMMON_DEFINITIONS.h"
#include "LEVELSET_3D.h"
#include "WORLD_DISCRETIZATION_3D.h"

// Subsolvers
#include "ADVECTION_3D.h"
#include "PROJECTION_3D.h"
#include "VISCOSITY_3D.h"

class EULERIAN_FLUID_SOLVER_3D
{
public: // Grids
	GRID_STRUCTURE_3D			base_grid;
	GRID_STRUCTURE_3D			ghost_grid;

	ARRAY<GRID_STRUCTURE_3D>	partial_base_grids;
	ARRAY<GRID_STRUCTURE_3D>	partial_ghost_grids;

public: // Dynamic data
	// Velocity Field for MAC grid
	FIELD_STRUCTURE_3D<T>*		water_velocity_field_mac_x;
	FIELD_STRUCTURE_3D<T>*		water_velocity_field_mac_y;
	FIELD_STRUCTURE_3D<T>*		water_velocity_field_mac_z;
	
	// Velocity Field for RK3rd
	FIELD_STRUCTURE_3D<T>		water_velocity_field_x_rk;
	FIELD_STRUCTURE_3D<T>		water_velocity_field_y_rk;
	FIELD_STRUCTURE_3D<T>		water_velocity_field_z_rk;

public: // Face value for levelset
	FIELD_STRUCTURE_3D<T>		levelset_x_c;
	FIELD_STRUCTURE_3D<T>		levelset_y_c;
	FIELD_STRUCTURE_3D<T>		levelset_z_c;

public: // Face value for density
	FIELD_STRUCTURE_3D<T>		density_half_x;
	FIELD_STRUCTURE_3D<T>		density_half_y;
	FIELD_STRUCTURE_3D<T>		density_half_z;

	// Density and Viscosity Field
	FIELD_STRUCTURE_3D<T>		density_field, viscosity_field;

	// Levelset for RK3rd
	FIELD_STRUCTURE_3D<T>		temp_for_levelset;

	// Water levelset
	LEVELSET_3D*				water_levelset;
	LEVELSET_3D*				boundary_levelset;

	LEVELSET_3D					object_levelset;

	WORLD_DISCRETIZATION_3D*	world_discretization;

public: // Ghost fields
	FIELD_STRUCTURE_3D<T>		scalar_field_ghost;
	FIELD_STRUCTURE_3D<T>		vector_field_mac_ghost_x;
	FIELD_STRUCTURE_3D<T>		vector_field_mac_ghost_y;
	FIELD_STRUCTURE_3D<T>		vector_field_mac_ghost_z;

public: // Simulation Options
	bool						air_water_simulation;
	bool						oil_water_simulation;
	bool						dimensionless_form;
	bool						is_vertical, is_parallel;
	bool						CSF_model;
	bool						regular_boundary, cylindrical_boundary;

public: // Simulation properties
	bool						use_water_solver;
	
	bool						fastsweeping_reinitialization;
	bool                        sussmanstyle_reinitialization;

	T							water_density, air_density, oil_density;
	T							water_viscosity, air_viscosity, oil_viscosity;
	VT							gravity;
	T							surface_tension;
	T							epsilon_for_mollification;

	bool						use_density_buoyancy;
	T	 						density_buoyancy_coefficient;

	int&						frame;

	T							water_height;

public: // Subsolvers
	ADVECTION_3D*				advecting_field_variables;
	PROJECTION_3D*				water_projection;
	VISCOSITY_3D*				viscosity_solver;

	// VERTEX_PARTICLE_METHOD		vertex_particle_method;

public: // Multithreading
	MULTITHREADING*				multithreading;						// multithreading should be initialized outside of this class for more efficient use of threads

public: // Option For Time Marching
	int							order_for_time_advancing;

public: // For debugging
	bool						is_Levelset_advection_active;
	bool						is_Velocity_advection_active;
	bool						is_external_force_active;
	bool						is_sourcing_active;
	bool						is_projection_active;

public: // Option for Viscosity Term
	bool						use_delta_function_formulation;

public: // Dimensionless Variables
	T							R1, R2, m, eta, K, g, f, a, A, V0, J, Re, We;
	T							one_over_We;

public: // Variables for Surface Tension Term
	FIELD_STRUCTURE_3D<T>		curv_x;
	FIELD_STRUCTURE_3D<T>		curv_y;
	FIELD_STRUCTURE_3D<T>		curv_z;
	FIELD_STRUCTURE_3D<T>		phi_c;
	FIELD_STRUCTURE_3D<T>		phi_l;
	FIELD_STRUCTURE_3D<T>		phi_d;
	FIELD_STRUCTURE_3D<T>		phi_b;
	FIELD_STRUCTURE_3D<T>		density_c;
	FIELD_STRUCTURE_3D<T>		density_l;
	FIELD_STRUCTURE_3D<T>		density_d;
	FIELD_STRUCTURE_3D<T>		density_b;
	FIELD_STRUCTURE_3D<T>		delta_c;
	FIELD_STRUCTURE_3D<T>		delta_l;
	FIELD_STRUCTURE_3D<T>		delta_d;
	FIELD_STRUCTURE_3D<T>		delta_b;
	
public: // For Boundary Condition of Cylindrical Pipe
	FIELD_STRUCTURE_3D<T>		boundary_phi_x;
	FIELD_STRUCTURE_3D<T>		boundary_phi_y;
	FIELD_STRUCTURE_3D<T>		boundary_phi_z;

public: // For CFL TimeStep Debugging
	T							c_f, g_f, s_f, v_f;

public: // Speedup Variable
	int							i_start_x, i_end_x, j_start_x, j_end_x, k_start_x, k_end_x;
	int							i_start_y, i_end_y, j_start_y, j_end_y, k_start_y, k_end_y;
	int							i_start_z, i_end_z, j_start_z, j_end_z, k_start_z, k_end_z;

public: // Scaling number for Reinitialization 
	int							scaling_number_for_reinitialization;
	int							iteration_number_for_reinitialization;

public: // Sign Function for Reinitialization
	FIELD_STRUCTURE_3D<T>		sign_function;

public: // Consructors and Destructor
	EULERIAN_FLUID_SOLVER_3D(int& frame_input)
		: frame(frame_input), multithreading(0), water_height(-numeric_limits<T>::max()), advecting_field_variables(0), water_projection(0), viscosity_solver(0),
		water_levelset(0), boundary_levelset(0), water_velocity_field_mac_x(0), water_velocity_field_mac_y(0), water_velocity_field_mac_z(0),
		use_water_solver(true), fastsweeping_reinitialization(false), sussmanstyle_reinitialization(false),
		is_Levelset_advection_active(false), is_Velocity_advection_active(false), is_external_force_active(false), is_sourcing_active(false), is_projection_active(false), order_for_time_advancing((int)1),
		use_delta_function_formulation(false), air_water_simulation(false), oil_water_simulation(false), dimensionless_form(false),
		is_vertical(false), is_parallel(false), CSF_model(false), c_f((T)0), g_f((T)0), s_f((T)0), v_f((T)0),
		regular_boundary(false), cylindrical_boundary(false)
	{}

	~EULERIAN_FLUID_SOLVER_3D(void)
	{
		// Note: Don't delete object_list here because is it initialized in SIMULATION_WORLD
		DELETE_POINTER(advecting_field_variables);

		if(use_water_solver)
		{
			DELETE_POINTER(water_projection);
			DELETE_POINTER(water_velocity_field_mac_x);
			DELETE_POINTER(water_velocity_field_mac_y);
			DELETE_POINTER(water_velocity_field_mac_z);
			DELETE_POINTER(water_levelset);
			DELETE_POINTER(boundary_levelset);
			DELETE_POINTER(viscosity_solver);
		}
	}

public: // Initialization Functions
	void InitializeFromSciptBlock(MULTITHREADING* multithreading_input, const SCRIPT_BLOCK& script_block)
	{
		// Initialize grids from outside
		SCRIPT_BLOCK grid_sb = script_block.FindBlock("GRID_STRUCTURE_3D");

		base_grid.InitializeFromBlock(grid_sb);

		InitializeFromScriptBlock(multithreading_input, base_grid, script_block);
	}

	void InitializeFromScriptBlock(MULTITHREADING* multithreading_input, const GRID_STRUCTURE_3D& world_grid, const SCRIPT_BLOCK& fluid_solver_block)
	{
		// Simulation properties from script
		use_water_solver = fluid_solver_block.GetBoolean("use_water_solver", (bool)true);
		
		int ghost_width = fluid_solver_block.GetInteger("ghost_width", 3);
		surface_tension = fluid_solver_block.GetFloat("surface_tension", (T)0);
		
		// Reinitialization
		fastsweeping_reinitialization = fluid_solver_block.GetBoolean("fastsweeping_reinitialization", (bool)false);
		sussmanstyle_reinitialization = fluid_solver_block.GetBoolean("sussmanstyle_reinitialization", (bool)false);
		scaling_number_for_reinitialization = fluid_solver_block.GetInteger("scaling_number_for_reinitialization", (int)5);
		iteration_number_for_reinitialization = fluid_solver_block.GetInteger("iteration_number_for_reinitialization", (int)10);
		
		// Option for Viscosity Term
		use_delta_function_formulation = fluid_solver_block.GetBoolean("use_delta_function_formulation", (bool)false);

		// Option for Surface Tension Term
		CSF_model = fluid_solver_block.GetBoolean("CSF_model", (bool)false);

		// Viscosity
		if (air_water_simulation)
		{
			water_viscosity = fluid_solver_block.GetFloat("water_viscosity", (T)1);
			air_viscosity = fluid_solver_block.GetFloat("air_viscosity", (T)0);
		}

		if (oil_water_simulation)
		{
			water_viscosity = fluid_solver_block.GetFloat("water_viscosity", (T)1);
			oil_viscosity = fluid_solver_block.GetFloat("oil_viscosity", (T)0);
		}
		
		gravity = fluid_solver_block.GetVector3("gravity", VT((T)0, (T)-9.8, (T)0));
		
		// Density
		if (air_water_simulation)
		{
			water_density = fluid_solver_block.GetFloat("water_density", 1);
			air_density = fluid_solver_block.GetFloat("air_density", 1);
		}
		if (oil_water_simulation)
		{
			water_density = fluid_solver_block.GetFloat("water_density", 1);
			oil_density = fluid_solver_block.GetFloat("oil_density", 1);
		}
		
		// Dimensionless Form
		if (oil_water_simulation)
		{
			dimensionless_form = fluid_solver_block.GetBoolean("dimensionless_form", (bool)false);
		}

		// For the pipe location
		if (oil_water_simulation)
		{
			is_vertical = world_discretization->is_vertical;
			is_parallel = world_discretization->is_parallel;
		}
		
		// Time marching order
		order_for_time_advancing = fluid_solver_block.GetInteger("order_for_time_advancing", (int)1);

		// Debugging 
		is_Levelset_advection_active = fluid_solver_block.GetBoolean("is_Levelset_advection_active", (bool)false);
		is_Velocity_advection_active = fluid_solver_block.GetBoolean("is_Velocity_advection_active", (bool)false);
		is_external_force_active = fluid_solver_block.GetBoolean("is_external_force_active", (bool)false);
		is_sourcing_active = fluid_solver_block.GetBoolean("is_sourcing_active", (bool)false);
		is_projection_active = fluid_solver_block.GetBoolean("is_projection_active", (bool)false);

		cout << "---------------DEBUGGING FOR SOLVER---------------" << endl;
		if (is_Levelset_advection_active)
		{
			cout << "Leveset Advection : true" << endl;
		}
		else
		{
			cout << "Leveset Advection : false" << endl;
		}

		if (is_Velocity_advection_active)
		{
			cout << "Velocity Advection : true" << endl;
		}
		else
		{
			cout << "Velocity Advection : false" << endl;
		}

		if (is_sourcing_active)
		{
			cout << "Sourcing : true" << endl;
		}
		else
		{
			cout << "Sourcing : false" << endl;
		}

		if (CSF_model)
		{
			cout << "---- CSF model: true" << endl;
		}
		else
		{
			cout << "---- CSF model: false" << endl;
		}

		if (is_projection_active)
		{
			cout << "Projection : true" << endl;
		}
		else
		{
			cout << "Projection : false" << endl;
		}

		cout << "--------------FLUID SOLVER VARIABLES--------------" << endl;
		if (use_water_solver)
		{
			cout << "use water solver : true" << endl;
		}
		else
		{
			cout << "use water solver : false" << endl;
		}

		cout << "ghost width : " << ghost_width << endl;
		cout << "gravity: (" << gravity.i << " ," << gravity.j << " ," << gravity.k << ") " << endl;
		
		if (air_water_simulation)
		{
			cout << "water density: " << water_density << endl;
			cout << "air density: " << air_density << endl;
			cout << "water viscosity: " << water_viscosity << endl;
			cout << "air viscosity: " << air_viscosity << endl;
			cout << "surface tension: " << surface_tension << endl;
		}
		if (oil_water_simulation)
		{
			cout << "water density: " << water_density << endl;
			cout << "oil density: " << oil_density << endl;
			cout << "water viscosity: " << water_viscosity << endl;
			cout << "oil viscosity: " << oil_viscosity << endl;
			cout << "surface tension: " << surface_tension << endl;
			
			if (dimensionless_form)
			{
				cout << "Dimensionless Form: true" << endl;
			}
			else
			{
				cout << "Dimensionless Form: false" << endl;
			}
		}

		switch (order_for_time_advancing)
		{
			case 1:
				cout << "Time advancing: Forward Euler" << endl;

			case 2:
				cout << "Time advancing: Runge-Kutta 2nd" << endl;
		
			case 3:
				cout << "Time advancing: Runge-Kutta 3rd" << endl;
		
			default:
				cout << "Time advancing: Forward Euler" << endl;
		}

		if (fastsweeping_reinitialization)
		{
			cout << "Reinitialized by Fast Sweeping Method: true" << endl;
		}
		else
		{
			cout << "Reinitialized by Fast Sweeping Method: false" << endl;
		}

		if (sussmanstyle_reinitialization)
		{
			cout << "Reinitialized by SussmanStyle: true" << endl;
		}
		else
		{
			cout << "Reinitialized by SussmanStyle: false" << endl;
		}

		if (air_water_simulation)
		{
			if (world_discretization->large_bubble)
			{
				cout << "Large Bubble Simulation is activated!" << endl;
			}
			if (world_discretization->small_bubble)
			{
				cout << "Small Bubble Simulation is activated!" << endl;
			}
		}
		
		base_grid.Initialize(world_grid, fluid_solver_block.GetFloat("resolution_scale", (T)1));
		ghost_grid.Initialize(base_grid.Enlarged(ghost_width));

		// Epsilon For Mollification -- Assume that the given grid is uniform
		epsilon_for_mollification = 1.5*base_grid.dx;

		// Multithreading 
		multithreading = multithreading_input;
		base_grid.SplitInZDirection(multithreading->num_threads, partial_base_grids);
		ghost_grid.SplitInZDirection(multithreading->num_threads, partial_ghost_grids);

		// Initialize data fields
		if (use_water_solver)
		{
			// Water Velocity Field for MAC grid
			DELETE_POINTER(water_velocity_field_mac_x);
			water_velocity_field_mac_x = new FIELD_STRUCTURE_3D<T>();
			water_velocity_field_mac_x->Initialize(multithreading, base_grid.i_res + 1, base_grid.j_res, base_grid.k_res, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min - (T)0.5*base_grid.dx, base_grid.y_min, base_grid.z_min, base_grid.x_max + (T)0.5*base_grid.dx, base_grid.y_max, base_grid.z_max, 2, false, true);
			water_velocity_field_mac_x->is_x_component = true;
			water_velocity_field_mac_x->is_y_component = false;
			water_velocity_field_mac_x->is_z_component = false;
			i_start_x = water_velocity_field_mac_x->i_start, i_end_x = water_velocity_field_mac_x->i_end, j_start_x = water_velocity_field_mac_x->j_start, j_end_x = water_velocity_field_mac_x->j_end, k_start_x = water_velocity_field_mac_x->k_start, k_end_x = water_velocity_field_mac_x->k_end;
			
			DELETE_POINTER(water_velocity_field_mac_y);
			water_velocity_field_mac_y = new FIELD_STRUCTURE_3D<T>();
			water_velocity_field_mac_y->Initialize(multithreading, base_grid.i_res, base_grid.j_res + 1, base_grid.k_res, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min, base_grid.y_min - (T)0.5*base_grid.dy, base_grid.z_min, base_grid.x_max, base_grid.y_max + (T)0.5*base_grid.dy, base_grid.z_max, 2, false, true);
			water_velocity_field_mac_y->is_x_component = false;
			water_velocity_field_mac_y->is_y_component = true;
			water_velocity_field_mac_y->is_z_component = false;
			i_start_y = water_velocity_field_mac_y->i_start, i_end_y = water_velocity_field_mac_y->i_end, j_start_y = water_velocity_field_mac_y->j_start, j_end_y = water_velocity_field_mac_y->j_end, k_start_y = water_velocity_field_mac_y->k_start, k_end_y = water_velocity_field_mac_y->k_end;

			DELETE_POINTER(water_velocity_field_mac_z);
			water_velocity_field_mac_z = new FIELD_STRUCTURE_3D<T>();
			water_velocity_field_mac_z->Initialize(multithreading, base_grid.i_res, base_grid.j_res, base_grid.k_res + 1, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min, base_grid.y_min, base_grid.z_min - (T)0.5*base_grid.dz, base_grid.x_max, base_grid.y_max, base_grid.z_max + (T)0.5*base_grid.dz, 2, false, true);
			water_velocity_field_mac_z->is_x_component = false;
			water_velocity_field_mac_z->is_y_component = false;
			water_velocity_field_mac_z->is_z_component = true;
			i_start_z = water_velocity_field_mac_z->i_start, i_end_z = water_velocity_field_mac_z->i_end, j_start_z = water_velocity_field_mac_z->j_start, j_end_z = water_velocity_field_mac_z->j_end, k_start_z = water_velocity_field_mac_z->k_start, k_end_z = water_velocity_field_mac_z->k_end;

			water_velocity_field_x_rk.Initialize(multithreading, base_grid.i_res + 1, base_grid.j_res, base_grid.k_res, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min - (T)0.5*base_grid.dx, base_grid.y_min, base_grid.z_min, base_grid.x_max + (T)0.5*base_grid.dx, base_grid.y_max, base_grid.z_max, 2, false, true);
			water_velocity_field_y_rk.Initialize(multithreading, base_grid.i_res, base_grid.j_res + 1, base_grid.k_res, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min, base_grid.y_min - (T)0.5*base_grid.dy, base_grid.z_min, base_grid.x_max, base_grid.y_max + (T)0.5*base_grid.dy, base_grid.z_max, 2, false, true);
			water_velocity_field_z_rk.Initialize(multithreading, base_grid.i_res, base_grid.j_res, base_grid.k_res + 1, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min, base_grid.y_min, base_grid.z_min - (T)0.5*base_grid.dz, base_grid.x_max, base_grid.y_max, base_grid.z_max + (T)0.5*base_grid.dz, 2, false, true);
		}
		
		object_levelset.Initialize(multithreading, base_grid, 2);				// Ghost cell is not required if we place wall conditions inside of simulation domain
		scalar_field_ghost.Initialize(multithreading, base_grid, 3, true, false);
		vector_field_mac_ghost_x.Initialize(multithreading, water_velocity_field_mac_x->grid, 3, false, true);
		vector_field_mac_ghost_y.Initialize(multithreading, water_velocity_field_mac_y->grid, 3, false, true);
		vector_field_mac_ghost_z.Initialize(multithreading, water_velocity_field_mac_z->grid, 3, false, true);
		
		scalar_field_ghost.is_scalar = true;
		
		vector_field_mac_ghost_x.is_x_component = true;
		vector_field_mac_ghost_x.is_y_component = false;
		vector_field_mac_ghost_x.is_z_component = false;
		
		vector_field_mac_ghost_y.is_x_component = false;
		vector_field_mac_ghost_y.is_y_component = true;
		vector_field_mac_ghost_y.is_z_component = false;

		vector_field_mac_ghost_z.is_x_component = false;
		vector_field_mac_ghost_z.is_y_component = false;
		vector_field_mac_ghost_z.is_z_component = true;

		// Water levelset
		if (use_water_solver)
		{
			DELETE_POINTER(water_levelset);
			water_levelset = new LEVELSET_3D();
			water_levelset->Initialize(multithreading, base_grid, 2);
			water_levelset->is_periodic = fluid_solver_block.GetBoolean("is_periodic", false);
			
			// For the boundary control
			DELETE_POINTER(boundary_levelset);
			boundary_levelset = new LEVELSET_3D();			
			boundary_levelset->Initialize(multithreading, base_grid, 2);
			boundary_levelset->is_periodic = fluid_solver_block.GetBoolean("is_periodic", false);

			if (oil_water_simulation)
			{
				if (is_vertical)
				{
					water_levelset->periodic_num_x = 1;
					water_levelset->periodic_num_y = fluid_solver_block.GetInteger("periodic_num", (int)1);
					boundary_levelset->periodic_num_x = 1;
					boundary_levelset->periodic_num_y = fluid_solver_block.GetInteger("periodic_num", (int)1);
				}
				if (is_parallel)
				{
					water_levelset->periodic_num_x = fluid_solver_block.GetInteger("periodic_num", (int)1);
					water_levelset->periodic_num_y = 1;
					boundary_levelset->periodic_num_x = fluid_solver_block.GetInteger("periodic_num", (int)1);
					boundary_levelset->periodic_num_y = 1;
				}
			}
		}

		// Initial values
		if (use_water_solver)
		{
			if (air_water_simulation)
			{
				water_velocity_field_mac_x->array_for_this.AssignAllValues(T());
				water_velocity_field_mac_y->array_for_this.AssignAllValues(T());
				water_velocity_field_mac_z->array_for_this.AssignAllValues(T());
			}
			if (oil_water_simulation)
			{
				water_velocity_field_mac_x->array_for_this.AssignAllValues(T());
				water_velocity_field_mac_y->array_for_this.AssignAllValues(T());
				water_velocity_field_mac_z->array_for_this.AssignAllValues(T());
				
				R1 = R2/a;
				m = fluid_solver_block.GetFloat("m", (T)1);
				eta = fluid_solver_block.GetFloat("eta", (T)1);
				K = fluid_solver_block.GetFloat("K", (T)1);
				g = fluid_solver_block.GetFloat("g", (T)1);
				f = (T)1/(K-(T)1)*g*(oil_density - K*water_density);
				A = (POW2(a) - (T)1) + (m*K + (T)2*(K - 1)*log(a));
				V0 = (f + water_density*g)*POW2(R1)/((T)4*water_viscosity)*A;
				J = surface_tension*R1*oil_density/POW2(oil_viscosity);
				//J = fluid_solver_block.GetFloat("J", (T)1);
				Re = (oil_density*R1*abs(V0))/oil_viscosity;
				//Re = fluid_solver_block.GetFloat("Re", (T)1.0);
				//V0 = Re*oil_viscosity/(oil_density*R1);
				We = ((T)1/J)*POW2(Re);
				//We = 11.3221;
				//We = oil_density*R1*POW2(V0)/surface_tension;
				one_over_We = (T)1/We;

				cout << "--------------PROPETIES--------------" << endl;
				cout << "a (R2/R1) = " << a << endl;
				cout << "R2 = " << R2 << endl;
				cout << "R1 = " << R1 << endl;
				cout << "m (mu2/mu1) = " << m << endl;
				cout << "K (f + rho_1*g)/(f + rho_2*g) = " << K << endl;
				cout << "g = " << g << endl;
				cout << "V0 = " << V0 << endl;
				cout << "f = " << f << endl;
				cout << "J = " << J << endl;
				cout << "We = " << We << endl;
				cout << "Re = " << Re << endl;
					
				int ghost_width_x;
				ghost_width_x = water_velocity_field_mac_x->ghost_width;
				
				int ghost_width_y;
				ghost_width_y = water_velocity_field_mac_y->ghost_width;

				if (dimensionless_form)
				{
					if (is_parallel)
					{
						for (int k = k_start_x - ghost_width_x; k <= k_end_x + ghost_width_x; k++)
						{
							for (int j = j_start_x - ghost_width_x; j <= j_end_x + ghost_width_x; j++)
							{
								for (int i = i_start_x - ghost_width_x; i <= i_end_x + ghost_width_x; i++)
								{
									T x_coor = water_velocity_field_mac_x->grid.x_min + i*water_velocity_field_mac_x->grid.dx, y_coor = water_velocity_field_mac_x->grid.y_min + j*water_velocity_field_mac_x->grid.dy, z_coor = water_velocity_field_mac_x->grid.z_min + k*water_velocity_field_mac_x->grid.dz;
									T magnitude = sqrt(POW2(y_coor) + POW2(z_coor));
									if ((magnitude >= 1) && (magnitude <= a))
									{
										water_velocity_field_mac_x->array_for_this(i, j, k) = (T)1/A*(POW2(a) - (POW2(y_coor) + POW2(z_coor)) - (T)2*(K - 1)*log(sqrt(POW2(y_coor) + POW2(z_coor))/a));
									}								
									else if (magnitude < 1)
									{
										water_velocity_field_mac_x->array_for_this(i, j, k) = (T)1 - (T)1/A*(m*(POW2(y_coor) + POW2(z_coor))*K);
									}
									else
									{
										water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0;
									}
								}
							}
						}
					}
					
					if (is_vertical)
					{
						for (int k = k_start_y; k <= k_end_y; k++)
						{
							for (int j = j_start_y; j <= j_end_y; j++)
							{
								for (int i = i_start_y; i <= i_end_y; i++)
								{
									T x_coor = water_velocity_field_mac_y->grid.x_min + i*water_velocity_field_mac_y->grid.dx, z_coor = water_velocity_field_mac_y->grid.z_min + k*water_velocity_field_mac_y->grid.dz;
									T magnitude = sqrt(POW2(x_coor) + POW2(z_coor));

									if ((magnitude >= 1) && (magnitude <= a))
									{
										if (V0 > 0)
										{
											water_velocity_field_mac_y->array_for_this(i, j, k) = (T)1/A*(POW2(a) - (POW2(x_coor) + POW2(z_coor)) - (T)2*(K - 1)*log(sqrt(POW2(x_coor) + POW2(z_coor))/a));
										}
										else
										{
											water_velocity_field_mac_y->array_for_this(i, j, k) = -(T)1/A*(POW2(a) - (POW2(x_coor) + POW2(z_coor)) - (T)2*(K - 1)*log(sqrt(POW2(x_coor) + POW2(z_coor))/a));
										}
									}								
									else if (magnitude < 1)
									{
										if (V0 > 0)
										{
											water_velocity_field_mac_y->array_for_this(i, j, k) = (T)1 - (T)1/A*(m*(POW2(x_coor) + POW2(z_coor))*K);
										}
										else
										{
											water_velocity_field_mac_y->array_for_this(i, j, k) = -((T)1 - (T)1/A*(m*(POW2(x_coor) + POW2(z_coor))*K));
										}
									}
									/*if (abs(magnitude) < 1e-8)
									{
										water_velocity_field_mac_y->array_for_this(i, j, k) = (T)1;
									}
									else
									{
										water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0;
									}*/
								}
							}
						}
												
						water_velocity_field_mac_x->array_for_this.AssignAllValues(T());
						water_velocity_field_mac_z->array_for_this.AssignAllValues(T());
					}
					
				}
				else
				{
					for (int k = k_start_x; k <= k_end_x; k++)
					{
						for (int j = j_start_x; j <= j_end_x; j++)
						{
							for (int i = i_start_x; i <= i_end_x; i++)
							{
								T x_coor = water_velocity_field_mac_x->grid.x_min + i*water_velocity_field_mac_x->grid.dx, y_coor = water_velocity_field_mac_x->grid.y_min + j*water_velocity_field_mac_x->grid.dy, z_coor = water_velocity_field_mac_x->grid.z_min + k*water_velocity_field_mac_x->grid.dz;
								if ((sqrt(POW2(y_coor) + POW2(z_coor)) >= R1) && (sqrt(POW2(y_coor) + POW2(z_coor)) <= R2))
								{
									water_velocity_field_mac_x->array_for_this(i, j, k) = (f + water_density*g)*(T)1/((T)4*water_viscosity)*(POW2(R2) - (POW2(y_coor) + POW2(z_coor)) - (T)2*POW2(R1)*(K - (T)1)*log(sqrt(POW2(y_coor) + POW2(z_coor))/R2));
								}								
								else if (sqrt(POW2(y_coor) + POW2(z_coor)) < R1)
								{
									water_velocity_field_mac_x->array_for_this(i, j, k) = (f + water_density*g)*(T)1/((T)4*water_viscosity)*(POW2(R1)*A - m*(POW2(y_coor) + POW2(z_coor))*K);
								}
							}
						}
					}
					water_velocity_field_mac_y->array_for_this.AssignAllValues(T());
					water_velocity_field_mac_z->array_for_this.AssignAllValues(T());
				}
			}
		}
		
		// Density and Viscosity Field
		viscosity_field.Initialize(multithreading, base_grid, 2, true, false);
		density_field.Initialize(multithreading, base_grid, 2, true, false);

		// Initialize advection solver
		DELETE_POINTER(advecting_field_variables);
		advecting_field_variables = new ADVECTION_3D(*multithreading, *water_levelset, scalar_field_ghost, *boundary_levelset, *water_velocity_field_mac_x, *water_velocity_field_mac_y, *water_velocity_field_mac_z, vector_field_mac_ghost_x, vector_field_mac_ghost_y, vector_field_mac_ghost_z);

		// Note : Three lines below should come before advecting_field_variables->InitializeFromBlock()
		advecting_field_variables->use_water_solver	= use_water_solver;
		advecting_field_variables->InitializeFromBlock(fluid_solver_block);
		advecting_field_variables->air_water_simulation = air_water_simulation;
		advecting_field_variables->oil_water_simulation = oil_water_simulation;
		
		if (oil_water_simulation)
		{
			advecting_field_variables->is_vertical = is_vertical;
			advecting_field_variables->is_parallel = is_parallel;
			advecting_field_variables->a = a;
		}
				
		// Initialize incompressibility solver
		if (use_water_solver)
		{
			DELETE_POINTER(water_projection);
			water_projection = new PROJECTION_3D(multithreading, *water_levelset, *water_velocity_field_mac_x, *water_velocity_field_mac_y, *water_velocity_field_mac_z, vector_field_mac_ghost_x, vector_field_mac_ghost_y, vector_field_mac_ghost_z, *boundary_levelset, frame, use_water_solver);
			
			if (air_water_simulation)
			{
				water_projection->air_water_simulation = true;
				water_projection->oil_water_simulation = false;
				water_projection->InitializeFromBlock(fluid_solver_block.FindBlock("PROJECTION_WATER"));
				water_projection->water_density = water_density;
				water_projection->air_density = air_density;
				water_projection->surface_tension = fluid_solver_block.GetFloat("surface_tension", (T)1);
				water_projection->use_delta_function_formulation = use_delta_function_formulation;
			}
			if (oil_water_simulation)
			{
				water_projection->air_water_simulation = false;
				water_projection->oil_water_simulation = true;
				water_projection->InitializeFromBlock(fluid_solver_block.FindBlock("PROJECTION_WATER"));
				water_projection->a = a;

				if (dimensionless_form)
				{
					water_projection->dimensionless_form = true;
					water_projection->We = We;
					water_projection->R1 = R1;
					water_projection->m = m;
					water_projection->eta = eta;
				}
				else
				{
					water_projection->dimensionless_form = false;
					water_projection->We = (T)0;
					water_projection->m = (T)1;
					water_projection->eta = (T)1;
				}

				water_projection->water_density = water_density;
				water_projection->oil_density = oil_density;
				water_projection->surface_tension = fluid_solver_block.GetFloat("surface_tension", (T)1);
				water_projection->use_delta_function_formulation = use_delta_function_formulation;
				regular_boundary = water_projection->regular_boundary;
				cylindrical_boundary = water_projection->cylindrical_boundary;
				water_projection->is_vertical = is_vertical;
				water_projection->is_parallel = is_parallel;
			}
		}

		// Initialize viscosity solver
		if (use_water_solver)
		{
			DELETE_POINTER(viscosity_solver);
			viscosity_solver = new VISCOSITY_3D(*multithreading, *water_levelset, *boundary_levelset, *water_velocity_field_mac_x, *water_velocity_field_mac_y, 
												*water_velocity_field_mac_z, scalar_field_ghost, vector_field_mac_ghost_x, vector_field_mac_ghost_y, vector_field_mac_ghost_z); 
			if (air_water_simulation)
			{
				viscosity_solver->air_water_simulation = true;
				viscosity_solver->oil_water_simulation = false;
			}
			if (oil_water_simulation)
			{
				viscosity_solver->air_water_simulation = false;
				viscosity_solver->oil_water_simulation = true;
				
				viscosity_solver->is_vertical = is_vertical;
				viscosity_solver->is_parallel = is_parallel;

				viscosity_solver->Re = Re;
				viscosity_solver->a = a;

				if (dimensionless_form)
				{
					viscosity_solver->dimensionless_form = true;
					viscosity_solver->m = m;
					viscosity_solver->eta = eta;
				}
				else
				{
					viscosity_solver->dimensionless_form = false;
					viscosity_solver->m = (T)1;
					viscosity_solver->eta = (T)1;
				}
			}

			viscosity_solver->InitializeFromScriptBlock(fluid_solver_block);
		}

		if (use_water_solver)
		{
			water_levelset->AssignAllValuesLevelsetThreaded(base_grid.dx*(T)3);
			water_levelset->FillGhostCellsFromThreaded(&(water_levelset->phi), false);
			water_levelset->curvature_by_normal_vector = fluid_solver_block.GetBoolean("curvature_by_normal_vector", false);
			water_levelset->curvature_by_levelset = fluid_solver_block.GetBoolean("curvature_by_levelset", false);
			boundary_levelset->AssignAllValuesLevelsetThreaded(base_grid.dx*(T)3);
			boundary_levelset->FillGhostCellsFromThreaded(&(boundary_levelset->phi), false);
		}

		// For CSF Model
		curv_x.Initialize(multithreading, water_velocity_field_mac_x->grid, 2);
		curv_y.Initialize(multithreading, water_velocity_field_mac_y->grid, 2);
		curv_z.Initialize(multithreading, water_velocity_field_mac_z->grid, 2);

		phi_c.Initialize(multithreading, water_levelset->grid, 2);
		phi_l.Initialize(multithreading, water_levelset->grid, 2);
		phi_d.Initialize(multithreading, water_levelset->grid, 2);
		phi_b.Initialize(multithreading, water_levelset->grid, 2);

		density_c.Initialize(multithreading, water_levelset->grid, 2);
		density_l.Initialize(multithreading, water_levelset->grid, 2);
		density_d.Initialize(multithreading, water_levelset->grid, 2);
		density_b.Initialize(multithreading, water_levelset->grid, 2);
		
		delta_c.Initialize(multithreading, water_levelset->grid, 2);
		delta_l.Initialize(multithreading, water_levelset->grid, 2);
		delta_d.Initialize(multithreading, water_levelset->grid, 2);
		delta_b.Initialize(multithreading, water_levelset->grid, 2);

		// For boundary condition of cylindrical pipe
		boundary_phi_x.Initialize(multithreading, water_velocity_field_mac_x->grid, 2);
		boundary_phi_y.Initialize(multithreading, water_velocity_field_mac_y->grid, 2);
		boundary_phi_z.Initialize(multithreading, water_velocity_field_mac_z->grid, 2);
		
		sign_function.Initialize(multithreading, water_levelset->signed_distance_field.grid, 3);

		// Face Value of levelset
		levelset_x_c.Initialize(multithreading, water_velocity_field_mac_x->grid, 2);
		levelset_y_c.Initialize(multithreading, water_velocity_field_mac_y->grid, 2);
		levelset_z_c.Initialize(multithreading, water_velocity_field_mac_z->grid, 2);

		// Face Value of Density
		density_half_x.Initialize(multithreading, water_velocity_field_mac_x->grid, 2);
		density_half_y.Initialize(multithreading, water_velocity_field_mac_y->grid, 2);
		density_half_z.Initialize(multithreading, water_velocity_field_mac_z->grid, 2);
	}

public: // Advancing
	void AdvanceOneTimeStepThread(const int& thread_id, const T& dt)
	{
		ofstream fout;
		if (order_for_time_advancing == 1)
		{
			Solve(thread_id, dt);
		}
		if (order_for_time_advancing == 2)
		{
			if (is_Levelset_advection_active)
			{
				BEGIN_HEAD_THREAD_WORK
				{
					temp_for_levelset.Initialize(multithreading, base_grid, 3);
				}
				END_HEAD_THREAD_WORK;

				BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
				{
					T temp = water_levelset->arr(i, j, k);
					temp_for_levelset(i, j, k) = temp;
				}
				END_GRID_ITERATION_3D;

				Levelset_Advection(thread_id, dt);
				Levelset_Advection(thread_id, dt);

				BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
				{
					water_levelset->arr(i, j, k) = (T)0.5*(temp_for_levelset(i, j, k) + water_levelset->arr(i, j, k));
				}
				END_GRID_ITERATION_3D;
			}

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				T temp_x = water_velocity_field_mac_x->array_for_this(i, j, k);
				water_velocity_field_x_rk(i, j, k) = temp_x;
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				T temp_y = water_velocity_field_mac_y->array_for_this(i, j, k);
				water_velocity_field_y_rk(i, j, k) = temp_y;
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				T temp_z = water_velocity_field_mac_z->array_for_this(i, j, k);
				water_velocity_field_z_rk(i, j, k) = temp_z;
			}
			END_GRID_ITERATION_3D;

			SolveForRK(thread_id, dt);
			SolveForRK(thread_id, dt);

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0.5*water_velocity_field_x_rk(i, j, k) + (T)0.5*water_velocity_field_mac_x->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0.5*water_velocity_field_y_rk(i, j, k) + (T)0.5*water_velocity_field_mac_y->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0.5*water_velocity_field_z_rk(i, j, k) + (T)0.5*water_velocity_field_mac_z->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;
			
			// Reinitialization by Forward Euler
			//Reinitialization(thread_id, dt);

			// Reinitialization using RK 2nd
			BEGIN_HEAD_THREAD_WORK
			{
				temp_for_levelset.Initialize(multithreading, base_grid, 3);
			}
			END_HEAD_THREAD_WORK;

			BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
			{
				T temp = water_levelset->arr(i, j, k);
				temp_for_levelset(i, j, k) = temp;
			}
			END_GRID_ITERATION_3D;

			Reinitialization(thread_id, dt);
			Reinitialization(thread_id, dt);

			BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
			{
				water_levelset->arr(i, j, k) = (T)0.5*(temp_for_levelset(i, j, k) + water_levelset->arr(i, j, k));
			}
			END_GRID_ITERATION_3D;

			if (is_vertical)
			{
				fout.open("water_levelset");
				for (int k = water_levelset->grid.k_start; k <= water_levelset->grid.k_end; k++)
				{
					for (int j = water_levelset->grid.j_start; j <= water_levelset->grid.j_end; j++)
					{
						fout << water_levelset->arr(water_levelset->grid.i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
			if (is_parallel)
			{
				fout.open("water_levelset");
				for (int i = water_levelset->grid.i_start; i <= water_levelset->grid.i_end; i++)
				{
					for (int j = water_levelset->grid.j_start; j <= water_levelset->grid.j_end; j++)
					{
						fout << water_levelset->arr(i, j, water_levelset->grid.k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
		}
		if (order_for_time_advancing == 3)
		{
			if (is_Levelset_advection_active)
			{
				Levelset_Advection(thread_id, dt);
			}

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				T temp_x = water_velocity_field_mac_x->array_for_this(i, j, k);
				water_velocity_field_x_rk(i, j, k) = temp_x;
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				T temp_y = water_velocity_field_mac_y->array_for_this(i, j, k);
				water_velocity_field_y_rk(i, j, k) = temp_y;
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				T temp_z = water_velocity_field_mac_z->array_for_this(i, j, k);
				water_velocity_field_z_rk(i, j, k) = temp_z;
			}
			END_GRID_ITERATION_3D;

			SolveForRK(thread_id, dt);
			SolveForRK(thread_id, dt);

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0.75*water_velocity_field_x_rk(i, j, k) + (T)0.25*water_velocity_field_mac_x->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0.75*water_velocity_field_y_rk(i, j, k) + (T)0.25*water_velocity_field_mac_y->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0.75*water_velocity_field_z_rk(i, j, k) + (T)0.25*water_velocity_field_mac_z->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;
			
			SolveForRK(thread_id, dt);

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_field_mac_x->array_for_this(i, j, k) = (T)1/(T)3*water_velocity_field_x_rk(i, j, k) + (T)2/(T)3*water_velocity_field_mac_x->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				water_velocity_field_mac_y->array_for_this(i, j, k) = (T)1/(T)3*water_velocity_field_y_rk(i, j, k) + (T)2/(T)3*water_velocity_field_mac_y->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				water_velocity_field_mac_z->array_for_this(i, j, k) = 1/(T)3*water_velocity_field_z_rk(i, j, k) + 2/(T)3*water_velocity_field_mac_z->array_for_this(i, j, k);
			}
			END_GRID_ITERATION_3D;
			
			// Reinitialization by Forward Euler
			//Reinitialization(thread_id, dt);

			// Reinitialization using RK 3rd
			BEGIN_HEAD_THREAD_WORK
			{
				temp_for_levelset.Initialize(multithreading, base_grid, 3);
			}
			END_HEAD_THREAD_WORK;

			BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
			{
				T temp = water_levelset->arr(i, j, k);
				temp_for_levelset(i, j, k) = temp;
			}
			END_GRID_ITERATION_3D;

			Reinitialization(thread_id, dt);
			Reinitialization(thread_id, dt);

			BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
			{
				water_levelset->arr(i, j, k) = (T)0.75*temp_for_levelset(i, j, k) + (T)0.25*water_levelset->arr(i, j, k);
			}
			END_GRID_ITERATION_3D;
			
			Reinitialization(thread_id, dt);

			BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
			{
				water_levelset->arr(i, j, k) = (T)1/(T)3*temp_for_levelset(i, j, k) + (T)2/(T)3*water_levelset->arr(i, j, k);
			}
			END_GRID_ITERATION_3D;
		}
	}

	void Solve(const int& thread_id, const T& dt)
	{
		ofstream fout;
		if (air_water_simulation)
		{
			if (is_Levelset_advection_active)
			{
				Levelset_Advection(thread_id, dt);
			}
		
			if (is_Velocity_advection_active)
			{
				Velocity_Advection(thread_id, dt);
			}

			if (is_external_force_active)
			{
				ApplyGravity(thread_id, dt);
			}
			
			if (is_sourcing_active)
			{
				Sourcing(thread_id, dt);
			}
		
			SetupBoundaryConditionForVelocity(thread_id);

			if (is_projection_active)
			{
				Projection(thread_id, dt);
			}
		}
		
		if (oil_water_simulation)
		{
			if (is_Levelset_advection_active)
			{
				Levelset_Advection(thread_id, dt);
			}

			if (is_Velocity_advection_active)
			{
				Velocity_Advection(thread_id, dt);
			}
						
			if (is_vertical)
			{
				fout.open("velocity_x_after_advection");
				for (int i = water_velocity_field_mac_x->i_start; i <= water_velocity_field_mac_x->i_end; i++)
				{
					for (int j = water_velocity_field_mac_x->j_start; j <= water_velocity_field_mac_x->j_end; j++)
					{
						fout << water_velocity_field_mac_x->array_for_this(i, j, water_velocity_field_mac_x->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_y_after_advection");
				for (int i = water_velocity_field_mac_y->i_start; i <= water_velocity_field_mac_y->i_end; i++)
				{
					for (int j = water_velocity_field_mac_y->j_start; j <= water_velocity_field_mac_y->j_end; j++)
					{
						fout << water_velocity_field_mac_y->array_for_this(i, j, water_velocity_field_mac_y->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_z_after_advection");
				for (int j = water_velocity_field_mac_z->j_start; j <= water_velocity_field_mac_z->j_end; j++)
				{
					for (int k = water_velocity_field_mac_z->k_start; k <= water_velocity_field_mac_z->k_end; k++)
					{
						fout << water_velocity_field_mac_z->array_for_this(water_velocity_field_mac_z->i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
			if (is_parallel)
			{
				fout.open("velocity_x_after_advection");
				for (int i = water_velocity_field_mac_x->i_start; i <= water_velocity_field_mac_x->i_end; i++)
				{
					for (int j = water_velocity_field_mac_x->j_start; j <= water_velocity_field_mac_x->j_end; j++)
					{
						fout << water_velocity_field_mac_x->array_for_this(i, j, water_velocity_field_mac_x->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_y_after_advection");
				for (int i = water_velocity_field_mac_y->i_start; i <= water_velocity_field_mac_y->i_end; i++)
				{
					for (int j = water_velocity_field_mac_y->j_start; j <= water_velocity_field_mac_y->j_end; j++)
					{
						fout << water_velocity_field_mac_y->array_for_this(i, j, water_velocity_field_mac_y->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_z_after_advection");
				for (int i = water_velocity_field_mac_z->i_start; i <= water_velocity_field_mac_z->i_end; i++)
				{
					for (int j = water_velocity_field_mac_z->j_start; j <= water_velocity_field_mac_z->j_end; j++)
					{
						fout << water_velocity_field_mac_z->array_for_this(i, j, water_velocity_field_mac_z->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}

			if (is_external_force_active)
			{
				ApplyGravity(thread_id, dt); 
			}

			if (is_sourcing_active)
			{
				Sourcing(thread_id, dt);
			}
			
			if (is_vertical)
			{
				fout.open("velocity_x_after_sourcing");
				for (int i = water_velocity_field_mac_x->i_start; i <= water_velocity_field_mac_x->i_end; i++)
				{
					for (int j = water_velocity_field_mac_x->j_start; j <= water_velocity_field_mac_x->j_end; j++)
					{
						fout << water_velocity_field_mac_x->array_for_this(i, j, water_velocity_field_mac_x->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_y_after_sourcing");
				for (int i = water_velocity_field_mac_y->i_start; i <= water_velocity_field_mac_y->i_end; i++)
				{
					for (int j = water_velocity_field_mac_y->j_start; j <= water_velocity_field_mac_y->j_end; j++)
					{
						fout << water_velocity_field_mac_y->array_for_this(i, j, water_velocity_field_mac_y->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_z_after_sourcing");
				for (int j = water_velocity_field_mac_z->j_start; j <= water_velocity_field_mac_z->j_end; j++)
				{
					for (int k = water_velocity_field_mac_z->k_start; k <= water_velocity_field_mac_z->k_end; k++)
					{
						fout << water_velocity_field_mac_z->array_for_this(water_velocity_field_mac_z->i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
			if (is_parallel)
			{
				fout.open("velocity_x_after_sourcing");
				for (int i = water_velocity_field_mac_x->i_start; i <= water_velocity_field_mac_x->i_end; i++)
				{
					for (int j = water_velocity_field_mac_x->j_start; j <= water_velocity_field_mac_x->j_end; j++)
					{
						fout << water_velocity_field_mac_x->array_for_this(i, j, water_velocity_field_mac_x->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_y_after_sourcing");
				for (int i = water_velocity_field_mac_y->i_start; i <= water_velocity_field_mac_y->i_end; i++)
				{
					for (int j = water_velocity_field_mac_y->j_start; j <= water_velocity_field_mac_y->j_end; j++)
					{
						fout << water_velocity_field_mac_y->array_for_this(i, j, water_velocity_field_mac_y->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_z_after_sourcing");
				for (int i = water_velocity_field_mac_z->i_start; i <= water_velocity_field_mac_z->i_end; i++)
				{
					for (int j = water_velocity_field_mac_z->j_start; j <= water_velocity_field_mac_z->j_end; j++)
					{
						fout << water_velocity_field_mac_z->array_for_this(i, j, water_velocity_field_mac_z->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}

			//SetupBoundaryConditionForVelocity(thread_id);

			if (is_projection_active)
			{
				Projection(thread_id, dt);
			}

			if (is_vertical)
			{
				fout.open("velocity_x_after_projection");
				for (int i = water_velocity_field_mac_x->i_start; i <= water_velocity_field_mac_x->i_end; i++)
				{
					for (int j = water_velocity_field_mac_x->j_start; j <= water_velocity_field_mac_x->j_end; j++)
					{
						fout << water_velocity_field_mac_x->array_for_this(i, j, water_velocity_field_mac_x->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_y_after_projection");
				for (int j = water_velocity_field_mac_y->j_start; j <= water_velocity_field_mac_y->j_end; j++)
				{
					for (int k = water_velocity_field_mac_y->k_start; k <= water_velocity_field_mac_y->k_end; k++)
					{
						fout << water_velocity_field_mac_y->array_for_this(water_velocity_field_mac_y->i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_z_after_projection");
				for (int j = water_velocity_field_mac_z->j_start; j <= water_velocity_field_mac_z->j_end; j++)
				{
					for (int k = water_velocity_field_mac_z->k_start; k <= water_velocity_field_mac_z->k_end; k++)
					{
						fout << water_velocity_field_mac_z->array_for_this(water_velocity_field_mac_z->i_end/2, j, k) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}

			if (is_parallel)
			{
				fout.open("velocity_x_after_projection");
				for (int i = water_velocity_field_mac_x->i_start; i <= water_velocity_field_mac_x->i_end; i++)
				{
					for (int j = water_velocity_field_mac_x->j_start; j <= water_velocity_field_mac_x->j_end; j++)
					{
						fout << water_velocity_field_mac_x->array_for_this(i, j, water_velocity_field_mac_x->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_y_after_projection");
				for (int i = water_velocity_field_mac_y->i_start; i <= water_velocity_field_mac_y->i_end; i++)
				{
					for (int j = water_velocity_field_mac_y->j_start; j <= water_velocity_field_mac_y->j_end; j++)
					{
						fout << water_velocity_field_mac_y->array_for_this(i, j, water_velocity_field_mac_y->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("velocity_z_after_projection");
				for (int i = water_velocity_field_mac_z->i_start; i <= water_velocity_field_mac_z->i_end; i++)
				{
					for (int j = water_velocity_field_mac_z->j_start; j <= water_velocity_field_mac_z->j_end; j++)
					{
						fout << water_velocity_field_mac_z->array_for_this(i, j, water_velocity_field_mac_z->k_end/2) << " ";
					}
					fout << "\n";
				}
				fout.close(); 
			}
		}
		// Reinitialization by Forward Euler
		//Reinitialization(thread_id, dt);

		// Reinitialization using RK 3rd
		BEGIN_HEAD_THREAD_WORK
		{
			temp_for_levelset.Initialize(multithreading, base_grid, 3);
		}
		END_HEAD_THREAD_WORK;

		BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
		{
			T temp = water_levelset->arr(i, j, k);
			temp_for_levelset(i, j, k) = temp;
		}
		END_GRID_ITERATION_3D;

		Reinitialization(thread_id, dt);
		Reinitialization(thread_id, dt);

		BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
		{
			water_levelset->arr(i, j, k) = (T)0.75*temp_for_levelset(i, j, k) + (T)0.25*water_levelset->arr(i, j, k);
		}
		END_GRID_ITERATION_3D;
		//water_levelset->FillGhostCellsContinuousDerivatesFromPointer(thread_id, &(water_levelset->phi), false);
		
		Reinitialization(thread_id, dt);

		BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
		{
			water_levelset->arr(i, j, k) = (T)1/(T)3*temp_for_levelset(i, j, k) + (T)2/(T)3*water_levelset->arr(i, j, k);
		}
		END_GRID_ITERATION_3D;
		//water_levelset->FillGhostCellsContinuousDerivatesFromPointer(thread_id, &(water_levelset->phi), false);
		
		if (is_vertical)
		{
			fout.open("water_levelset");
			for (int k = water_levelset->grid.k_start; k <= water_levelset->grid.k_end; k++)
			{
				for (int j = water_levelset->grid.j_start; j <= water_levelset->grid.j_end; j++)
				{
					fout << water_levelset->arr(water_levelset->grid.i_end/2, j, k) << " ";
				}
				fout << "\n";
			}
			fout.close(); 
		}
		if (is_parallel)
		{
			fout.open("water_levelset");
			for (int i = water_levelset->grid.i_start; i <= water_levelset->grid.i_end; i++)
			{
				for (int j = water_levelset->grid.j_start; j <= water_levelset->grid.j_end; j++)
				{
					fout << water_levelset->arr(i, j, water_levelset->grid.k_end/2) << " ";
				}
				fout << "\n";
			}
			fout.close(); 
		}
	}

	void SolveForRK(const int& thread_id, const T& dt)
	{
		if (is_Velocity_advection_active)
		{
			Velocity_Advection(thread_id, dt);
		}
		
		if (is_external_force_active)
		{
			ApplyGravity(thread_id, dt);		 
		}

		if (is_sourcing_active)
		{
			Sourcing(thread_id, dt);
		}
				
		//SetupBoundaryConditionForVelocity(thread_id);

		if (is_projection_active)
		{
			Projection(thread_id, dt);
		}
	}
	
public: // Simulation steps
	void Levelset_Advection(const int& thread_id, const T& dt)
	{
		LOG::Begin(thread_id, "Levelset Advection");
		
		advecting_field_variables->Solve_Levelset(thread_id, dt);

		LOG::End(thread_id);
	}
	
	void Velocity_Advection(const int& thread_id, const T& dt)
	{
		LOG::Begin(thread_id, "Velocity Advection");

		advecting_field_variables->Solve_Velocity(thread_id, dt);

		LOG::End(thread_id);
	}

	void DetermineDensity(const int& thread_id)
	{
		if (oil_water_simulation)
		{
			BEGIN_GRID_ITERATION_3D(density_field.partial_grids[thread_id])
			{
				if (dimensionless_form)
				{
					if (water_levelset->arr(i, j, k) <= 0)
					{
						density_field(i, j, k) = 1;
					}
					else if (water_levelset->arr(i, j, k) > 0)
					{
						density_field(i, j, k) = eta;
					}
				}
				else
				{
					if (water_levelset->arr(i, j, k) <= 0)
					{
						density_field(i, j, k) = oil_density;
					}
					else if (water_levelset->arr(i, j, k) > 0)
					{
						density_field(i, j, k) = water_density;
					}
				}
			}

			density_field.FillGhostCellsFrom(thread_id, density_field.array_for_this, true);
			
			END_GRID_ITERATION_3D;
		}
	}
	
	void SetupBoundaryConditionForVelocity(void)
	{
		if (air_water_simulation)
		{
			if (water_projection->air_bubble_rising)
			{
				//// Front and Rear Face
				//BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
				//{
				//	water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)0;
				//	water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)0;
				//}
				//END_GRID_ITERATION_3D;

				//BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
				//{
				//	water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
				//	water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
				//}
				//END_GRID_ITERATION_3D;

				//BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
				//{
				//	water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
				//	water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
				//}
				//END_GRID_ITERATION_3D;

				for (int j = j_start_x; j <= j_end_x; j++)
				{
					for (int k = k_start_x; k <= k_end_x; k++)
					{
						water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)0;
						water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)0;
					}
				}

				// Right and Left Face
				for (int i = i_start_y; i <= i_end_y; i++)
				{
					for (int k = k_start_y; k <= k_end_y; k++)
					{
						water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
						water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
					}
				}

				// Upper and Bottom Face
				for (int i = i_start_z; i <= i_end_z; i++)
				{
					for (int j = j_start_z; j <= j_end_z; j++)
					{
						water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
						water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
					}
				}
			}
			if (oil_water_simulation)
			{
				if (is_parallel)
				{
					if (regular_boundary)
					{
						// Front and Rear Face
						for (int j = j_start_x; j <= j_end_x; j++)
						{
							for (int k = k_start_x; k <= k_end_x; k++)
							{
								water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_start_x + 1, j, k) - water_velocity_field_mac_x->array_for_this(i_start_x + 2, j, k);
								water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_end_x - 1, j, k) - water_velocity_field_mac_x->array_for_this(i_end_x - 2, j, k);
							}
						}

						// Right and Left Face
						for (int i = i_start_y; i <= i_end_y; i++)
						{
							for (int k = k_start_y; k <= k_end_y; k++)
							{
								water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
								water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
							}
						}

						// Upper and Bottom Face
						for (int i = i_start_z; i <= i_end_z; i++)
						{
							for (int j = j_start_z; j <= j_end_z; j++)
							{
								water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
								water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
							}
						}
					}
					if (cylindrical_boundary)
					{
						GRID_ITERATION_3D(boundary_phi_x.grid)
						{
							if (boundary_phi_x(i, j, k) > 0)
							{
								water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0;
							}
							else
							{
								water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_start_x + 1, j, k) - water_velocity_field_mac_x->array_for_this(i_start_x + 2, j, k);
								water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_end_x - 1, j, k) - water_velocity_field_mac_x->array_for_this(i_end_x - 2, j, k);
							}
						}

						GRID_ITERATION_3D(boundary_phi_y.grid)
						{
							if (boundary_phi_y(i, j, k) > 0)
							{
								water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0;
							}
							else
							{
								water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
								water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
							}
						}
						
						GRID_ITERATION_3D(boundary_phi_z.grid)
						{
							if (boundary_phi_z(i, j, k) > 0)
							{
								water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0;
							}
							else
							{
								water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
								water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
							}
						}
					}
				}

				if (is_vertical)
				{
					if (regular_boundary)
					{
						// Front and Rear Face
						for (int j = j_start_x; j <= j_end_x; j++)
						{
							for (int k = k_start_x; k <= k_end_x; k++)
							{
								water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)0;
								water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)0;
							}
						}
	
						// Right and Left Face
						for (int i = i_start_y; i <= i_end_y; i++)
						{
							for (int k = k_start_y; k <= k_end_y; k++)
							{
								water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_start_y + 1, k) - water_velocity_field_mac_y->array_for_this(i, j_start_y + 2, k);
								water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_end_y - 1, k) - water_velocity_field_mac_y->array_for_this(i, j_end_y - 2, k);
							}
						}
	
						// Upper and Bottom Face
						for (int i = i_start_z; i <= i_end_z; i++)
						{
							for (int j = j_start_z; j <= j_end_z; j++)
							{
								water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
								water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
							}
						}
					}
					if (cylindrical_boundary)
					{
						GRID_ITERATION_3D(boundary_phi_x.grid)
						{
							if (boundary_phi_x(i, j, k) > 0)
							{
								water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0;
							}
							else
							{
								water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)0;
								water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)0;
							}
						}

						GRID_ITERATION_3D(boundary_phi_y.grid)
						{
							if (boundary_phi_y(i, j, k) > 0)
							{
								water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0;
							}
							else
							{
								T x_coor = water_velocity_field_mac_y->grid.x_min + i*water_velocity_field_mac_y->grid.dx, y_coor = water_velocity_field_mac_y->grid.y_min + j*water_velocity_field_mac_y->grid.dy, z_coor = water_velocity_field_mac_y->grid.z_min + k*water_velocity_field_mac_y->grid.dz;
								if ((sqrt(POW2(x_coor) + POW2(z_coor)) >= 1) && (sqrt(POW2(x_coor) + POW2(z_coor)) <= a))
								{
									water_velocity_field_mac_y->array_for_this(i, j_start_y, k) =(T)1/A*(POW2(a) - (POW2(x_coor) + POW2(z_coor)) - (T)2*(K - 1)*log(sqrt(POW2(x_coor) + POW2(z_coor))/a));
									water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_end_y - 1, k) - water_velocity_field_mac_y->array_for_this(i, j_end_y - 2, k);
								}								
								else if (sqrt(POW2(x_coor) + POW2(z_coor)) < 1)
								{
									water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)1 - (T)1/A*(m*(POW2(x_coor) + POW2(z_coor))*K);
									water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_end_y - 1, k) - water_velocity_field_mac_y->array_for_this(i, j_end_y - 2, k);
								}
								else
								{
									water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
									water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_end_y - 1, k) - water_velocity_field_mac_y->array_for_this(i, j_end_y - 2, k);
								}
							}
						}
						
						GRID_ITERATION_3D(boundary_phi_z.grid)
						{
							if (boundary_phi_z(i, j, k) > 0)
							{
								water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0;
							}
							else
							{
								water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
								water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
							}
						}
					}
				}
			}
		}
	}

	void SetupBoundaryConditionForVelocity(const int& thread_id)
	{
		if (air_water_simulation)
		{
			if (water_projection->air_bubble_rising)
			{
				// Front and Rear Face
				BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
				{
					water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)0;
					water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)0;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
				{
					water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
					water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
				{
					water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
					water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
				}
				END_GRID_ITERATION_3D;
			}
		}
		
		if (oil_water_simulation)
		{
			if (is_parallel)
			{
				if (regular_boundary)
				{
					// Front and Rear Face
					for (int j = j_start_x; j <= j_end_x; j++)
					{
						for (int k = k_start_x; k <= k_end_x; k++)
						{
							water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_start_x + 1, j, k) - water_velocity_field_mac_x->array_for_this(i_start_x + 2, j, k);
							water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_end_x - 1, j, k) - water_velocity_field_mac_x->array_for_this(i_end_x - 2, j, k);
						}
					}
					
					// Right and Left Face
					for (int i = i_start_y; i <= i_end_y; i++)
					{
						for (int k = k_start_y; k <= k_end_y; k++)
						{
							water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
							water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
						}
					}

					// Upper and Bottom Face
					for (int i = i_start_z; i <= i_end_z; i++)
					{
						for (int j = j_start_z; j <= j_end_z; j++)
						{
							water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
							water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
						}
					}
				}
				if (cylindrical_boundary)
				{
					GRID_ITERATION_3D(boundary_phi_x.grid)
					{
						if (boundary_phi_x(i, j, k) > 0)
						{
							water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0;
						}
						else
						{
							water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_start_x + 1, j, k) - water_velocity_field_mac_x->array_for_this(i_start_x + 2, j, k);
							water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)2*water_velocity_field_mac_x->array_for_this(i_end_x - 1, j, k) - water_velocity_field_mac_x->array_for_this(i_end_x - 2, j, k);
						}
					}

					GRID_ITERATION_3D(boundary_phi_y.grid)
					{
						if (boundary_phi_y(i, j, k) > 0)
						{
							water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0;
						}
						else
						{
							water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
							water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
						}
					}
						
					GRID_ITERATION_3D(boundary_phi_z.grid)
					{
						if (boundary_phi_z(i, j, k) > 0)
						{
							water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0;
						}
						else
						{
							water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)0;
							water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)0;
						}
					}
				}
			}

			if (is_vertical)
			{
				if (regular_boundary)
				{
					// Front and Rear Face
					for (int j = j_start_x; j <= j_end_x; j++)
					{
						for (int k = k_start_x; k <= k_end_x; k++)
						{
							water_velocity_field_mac_x->array_for_this(i_start_x, j, k) = (T)0;
							water_velocity_field_mac_x->array_for_this(i_end_x, j, k) = (T)0;
						}
					}

					// Right and Left Face
					for (int i = i_start_y; i <= i_end_y; i++)
					{
						for (int k = k_start_y; k <= k_end_y; k++)
						{
							water_velocity_field_mac_y->array_for_this(i, j_start_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_start_y + 1, k) - water_velocity_field_mac_y->array_for_this(i, j_start_y + 2, k);
							water_velocity_field_mac_y->array_for_this(i, j_end_y, k) = (T)2*water_velocity_field_mac_y->array_for_this(i, j_end_y - 1, k) - water_velocity_field_mac_y->array_for_this(i, j_end_y - 2, k);
						}
					}
	
					// Upper and Bottom Face
					for (int i = i_start_z; i <= i_end_z; i++)
					{
						for (int j = j_start_z; j <= j_end_z; j++)
						{
							water_velocity_field_mac_z->array_for_this(i, j, k_start_z) = (T)0;
							water_velocity_field_mac_z->array_for_this(i, j, k_end_z) = (T)0;
						}
					}
				}
				if (cylindrical_boundary)
				{
					if (is_vertical)
					{
						BEGIN_GRID_ITERATION_3D(boundary_levelset->partial_grids[thread_id])
						{
							if (boundary_levelset->arr(i, j, k) > 0)
							{
								if (boundary_levelset->arr(i - 1, j, k) < 0)
								{
									water_velocity_field_mac_x->array_for_this(i, j, k) = (T)0;
								}
								if (boundary_levelset->arr(i + 1, j, k) < 0)
								{
									water_velocity_field_mac_x->array_for_this(i + 1, j, k) = (T)0;
								}
								if (boundary_levelset->arr(i, j, k - 1) < 0)
								{
									water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0;
								}
								if (boundary_levelset->arr(i, j, k + 1) < 0)
								{
									water_velocity_field_mac_z->array_for_this(i, j, k + 1) = (T)0;
								}
							}
						}
						END_GRID_ITERATION_3D; 
					}
					if (is_parallel)
					{
						BEGIN_GRID_ITERATION_3D(boundary_levelset->partial_grids[thread_id])
						{
							if (boundary_levelset->arr(i, j, k) > 0)
							{
								if (boundary_levelset->arr(i, j - 1, k) < 0)
								{
									water_velocity_field_mac_y->array_for_this(i, j, k) = (T)0;
								}
								if (boundary_levelset->arr(i, j + 1, k) < 0)
								{
									water_velocity_field_mac_z->array_for_this(i + 1, j, k) = (T)0;
								}
								if (boundary_levelset->arr(i, j, k - 1) < 0)
								{
									water_velocity_field_mac_z->array_for_this(i, j, k) = (T)0;
								}
								if (boundary_levelset->arr(i, j, k + 1) < 0)
								{
									water_velocity_field_mac_z->array_for_this(i, j, k + 1) = (T)0;
								}
							}
						}
						END_GRID_ITERATION_3D; 
					}
				}
			}
		}
	}
	
	void Sourcing(const int& thread_id, const T& dt)
	{
		LOG::Begin(thread_id, "Sourcing");

		DetermineDensity(thread_id);

		viscosity_solver->ApplyViscosity(thread_id, epsilon_for_mollification, dt);
		
		if (world_discretization->object_list.num_of_elements > 0)
		{
			SourcingFromObjectList(thread_id);
		}

		if (CSF_model)
		{
			ApplySurfaceTension(thread_id, dt);
			water_projection->CSF_model = true;
		}

		//if (sussmanstyle_reinitialization)
		//{
		//	water_levelset->ReinitializeBySussmanStyle(thread_id);
		//}
		LOG::End(thread_id);
	}
	
	void ApplySurfaceTension(const int& thread_id, const T& dt)
	{
		water_levelset->ComputeCurvatures(thread_id);
		
		if (oil_water_simulation)
		{
			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				curv_x(i, j, k) = (T)0.5*(water_levelset->curvature(i, j, k) + water_levelset->curvature(i - 1, j, k));
			}
			
			curv_x.FillGhostCellsContinuousDerivativesFrom(thread_id, curv_x.array_for_this, true);
			
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				curv_y(i, j, k) = (T)0.5*(water_levelset->curvature(i, j, k) + water_levelset->curvature(i, j - 1, k));
			}
			
			curv_y.FillGhostCellsContinuousDerivativesFrom(thread_id, curv_y.array_for_this, true);

			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				curv_z(i, j, k) = (T)0.5*(water_levelset->curvature(i, j, k) + water_levelset->curvature(i, j, k - 1));
			}
			
			curv_z.FillGhostCellsContinuousDerivativesFrom(thread_id, curv_z.array_for_this, true);

			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_levelset->partial_grids[thread_id])
			{
				phi_c(i, j, k) = water_levelset->arr(i, j, k);
				phi_l(i, j, k) = water_levelset->arr(i - 1, j, k);
				phi_d(i, j, k) = water_levelset->arr(i, j - 1, k);
				phi_b(i, j, k) = water_levelset->arr(i, j, k - 1);
				density_c(i, j, k) = density_field(i, j, k);
				density_l(i, j, k) = density_field(i - 1, j, k);
				density_d(i, j, k) = density_field(i, j - 1, k);
				density_b(i, j, k) = density_field(i, j, k - 1);
			}
			
			phi_c.FillGhostCellsContinuousDerivativesFrom(thread_id, phi_c.array_for_this, true);
			phi_l.FillGhostCellsContinuousDerivativesFrom(thread_id, phi_l.array_for_this, true);
			phi_d.FillGhostCellsContinuousDerivativesFrom(thread_id, phi_d.array_for_this, true);
			phi_b.FillGhostCellsContinuousDerivativesFrom(thread_id, phi_b.array_for_this, true);
			density_c.FillGhostCellsFrom(thread_id, density_c.array_for_this, true);
			density_l.FillGhostCellsFrom(thread_id, density_l.array_for_this, true);
			density_d.FillGhostCellsFrom(thread_id, density_d.array_for_this, true);
			density_b.FillGhostCellsFrom(thread_id, density_b.array_for_this, true);

			END_GRID_ITERATION_3D;

			DeltaFunction(thread_id, phi_c, epsilon_for_mollification, delta_c);
			DeltaFunction(thread_id, phi_l, epsilon_for_mollification, delta_l);
			DeltaFunction(thread_id, phi_d, epsilon_for_mollification, delta_d);
			DeltaFunction(thread_id, phi_b, epsilon_for_mollification, delta_b);

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				T delta_x = (T)0.5*(delta_c(i, j, k) + delta_l(i, j, k));
				T density_x = (T)0.5*(density_c(i, j, k) + density_l(i, j, k));

				T surf_x = one_over_We*curv_x(i, j, k)*delta_x*water_velocity_field_mac_x->one_over_dx*(phi_c(i, j, k) - phi_l(i, j, k))/density_x;
				
				water_velocity_field_mac_x->array_for_this(i, j, k) += dt*surf_x;
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				T delta_y = (T)0.5*(delta_c(i, j, k) + delta_d(i, j, k));
				T density_y = (T)0.5*(density_c(i, j, k) + density_d(i, j, k));

				T surf_y = one_over_We*curv_y(i, j, k)*delta_y*water_velocity_field_mac_y->one_over_dy*(phi_c(i, j, k) - phi_d(i, j, k))/density_y;

				water_velocity_field_mac_y->array_for_this(i, j, k) += dt*surf_y;
			}
			END_GRID_ITERATION_3D;

			BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
			{
				T delta_z = (T)0.5*(delta_c(i, j, k) + delta_b(i, j, k));
				T density_z = (T)0.5*(density_c(i, j, k) + density_b(i, j, k));

				T surf_z = one_over_We*curv_z(i, j, k)*delta_z*water_velocity_field_mac_z->one_over_dz*(phi_c(i, j, k) - phi_b(i, j, k))/density_z;

				water_velocity_field_mac_z->array_for_this(i, j, k) += dt*surf_z;
			}
			END_GRID_ITERATION_3D;
		}
	}

	void Reinitialization(const int& thread_id, const T& dt)
	{
		LOG::Begin(thread_id, "Reinitialization");

		if (fastsweeping_reinitialization)
		{
			water_levelset->ReinitializeByFastSweeping(thread_id);
		}
		if (sussmanstyle_reinitialization)
		{
			T fictious_dt = base_grid.dx/scaling_number_for_reinitialization;
			int iter(0);

			BEGIN_GRID_ITERATION_3D(sign_function.partial_grids[thread_id])
			{
				sign_function(i, j, k) = water_levelset->signed_distance_field(i, j, k)/sqrt(POW2(water_levelset->signed_distance_field(i, j, k)) + POW2(water_levelset->grid.dx));
			}
			END_GRID_ITERATION_3D;

			while (iter < iteration_number_for_reinitialization)
			{
				advecting_field_variables->ReinitializationBySussman(thread_id, fictious_dt, sign_function);
				iter = iter + 1;
			}
		}

		LOG::End(thread_id);
	}
	void Projection(const int& thread_id, const T& dt)
	{
		if (use_water_solver)
		{
			water_projection->surface_tension = surface_tension;
			water_projection->Solve(thread_id, dt);
			
			// water_projection->FastSweepingExtrapolation<VT>(thread_id, *water_velocity_field, vector_field_ghost, object_levelset);
		}
	}

	void ApplyGravity(const int& thread_id, const T& dt)
	{
		const VT gravity_dt(dt*gravity);
		
		/*DetermineDensityField(thread_id, water_levelset->signed_distance_field, 1.5*base_grid.dx, density_field); 
		multithreading->Sync(thread_id);*/

		// Gravity
		if (use_water_solver)
		{
			if (air_water_simulation)
			{
				ARRAY_3D<T>& water_velocity_x(water_velocity_field_mac_x->array_for_this);
				ARRAY_3D<T>& water_velocity_y(water_velocity_field_mac_y->array_for_this);
				ARRAY_3D<T>& water_velocity_z(water_velocity_field_mac_z->array_for_this);

				BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
				{
					water_velocity_x(i, j, k) += gravity_dt.x;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
				{
					water_velocity_y(i, j, k) += gravity_dt.y;
				}
				END_GRID_ITERATION_3D;

				BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
				{
					water_velocity_z(i, j, k) += gravity_dt.z;
				}
				END_GRID_ITERATION_3D;
			}
			if (oil_water_simulation)
			{
				ARRAY_3D<T>& water_velocity_x(water_velocity_field_mac_x->array_for_this);
				ARRAY_3D<T>& water_velocity_y(water_velocity_field_mac_y->array_for_this);
				ARRAY_3D<T>& water_velocity_z(water_velocity_field_mac_z->array_for_this);
								
				if (dimensionless_form)
				{
					// Set up the cell face value of boundary levelset
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
						const VT driving_pressure(dt*VT(0, f, 0));

						scalar_field_ghost.FillGhostCellsPeriodicInYDirection(thread_id, water_levelset->arr, true);

						// Determine the cell face levelset and density
						BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
						{
							/*levelset_y_c(i, j, k) = (T)0.5*(water_levelset->arr(i, j, k) + water_levelset->arr(i, j - 1, k));*/
							levelset_y_c(i, j, k) = (T)0.5*(scalar_field_ghost(i, j, k) + scalar_field_ghost(i, j - 1, k));
						}
						END_GRID_ITERATION_3D;

						DetermineDensityField(thread_id, levelset_y_c, (T)1.5*base_grid.dy, density_half_y); 

						// Apply gravity only inside the cylinder domain
						BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
						{
							if (water_velocity_field_mac_y->fixed(i, j, k) == true)
							{
								continue;
							}
							if (boundary_phi_y(i, j, k) < 0)
							{
								// External force including driving pressure
								water_velocity_y(i, j, k) += R1/POW2(V0)*(gravity_dt.y + (T)1/(oil_density*density_half_y(i, j, k))*driving_pressure.y);
							}
							else
							{
								continue;
							}
						}
						END_GRID_ITERATION_3D;
					}
					if (is_parallel)
					{
						const VT driving_pressure(dt*VT(f, 0, 0));

						scalar_field_ghost.FillGhostCellsPeriodicInXDirection(thread_id, water_levelset->arr, true);

						// Determine the cell face levelset and density
						BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
						{
							levelset_x_c(i, j, k) = (T)0.5*(scalar_field_ghost(i, j, k) + scalar_field_ghost(i - 1, j, k));
						}
						END_GRID_ITERATION_3D;

						DetermineDensityField(thread_id, levelset_x_c, (T)1.5*base_grid.dx, density_half_x); 

						BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
						{
							if (water_velocity_field_mac_x->fixed(i, j, k) == true)
							{
								continue;
							}

							if (boundary_phi_x(i, j, k) < 0)
							{
								// External force including driving pressure
								water_velocity_x(i, j, k) += R1/POW2(V0)*(gravity_dt.x + ((T)1/(density_half_x(i, j, k)*oil_density))*driving_pressure.x);
							}
							else
							{
								continue;
							}
						}
						END_GRID_ITERATION_3D;

						BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
						{
							if (water_velocity_field_mac_y->fixed(i, j, k) == true)
							{
								continue;
							}

							if (boundary_phi_y(i, j, k) < 0)
							{
								//continue;
								water_velocity_y(i, j, k) += R1/POW2(V0)*gravity_dt.y;
							}
							else
							{
								continue;
							}
						}
						END_GRID_ITERATION_3D;
	
						BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
						{
							if (water_velocity_field_mac_z->fixed(i, j, k) == true)
							{
								continue;
							}

							if (boundary_phi_z(i, j, k) < 0)
							{
								water_velocity_z(i, j, k) += R1/POW2(V0)*gravity_dt.z; 
							}
							else
							{
								continue;
							}
						}
						END_GRID_ITERATION_3D;
					}
				}
				else
				{
					BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
					{
						water_velocity_x(i, j, k) += gravity_dt.x;
					}
					END_GRID_ITERATION_3D;

					BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
					{
						water_velocity_y(i, j, k) += gravity_dt.y;
					}
					END_GRID_ITERATION_3D;

					BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
					{
						water_velocity_z(i, j, k) += gravity_dt.z;
					}
					END_GRID_ITERATION_3D;
				}	
			}
		}
	}
		
public: // Sourcing Functions
	void HeavisideFunction(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& heaviside)
	{
		T one_over_epsilon = (T)1/epsilon, one_over_pi = (T)1/PI;

		BEGIN_GRID_ITERATION_3D(phi.partial_grids[thread_id])
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
				heaviside(i, j, k) = (T)0.5 + phi(i, j, k)*one_over_epsilon*(T)0.5 + (T)0.5*one_over_pi*sin((T)PI*phi(i, j, k)*one_over_epsilon);
			}
		}
		END_GRID_ITERATION_3D;
	}

	void DetermineViscosityField(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& viscosity)
	{
		FIELD_STRUCTURE_3D<T> heaviside_phi;
		heaviside_phi.Initialize(multithreading, viscosity.grid, 2);

		HeavisideFunction(thread_id, phi, epsilon, heaviside_phi);

		multithreading->Sync(thread_id);

		GRID_ITERATION_3D(phi.partial_grids[thread_id])
		{
			if (air_water_simulation)
			{
				if (water_projection->air_bubble_rising == true)
				{
					viscosity(i, j, k) = air_viscosity + (water_viscosity - air_viscosity)*heaviside_phi(i, j, k);
				}
				if (water_projection->water_drop == true)
				{
					viscosity(i, j, k) = water_viscosity + (air_viscosity - water_viscosity)*heaviside_phi(i, j, k);
				}
			}
			if (oil_water_simulation)
			{
				if (dimensionless_form)
				{
					//viscosity(i, j, k) = oil_viscosity/water_viscosity + ((T)1 - oil_viscosity/water_viscosity)*heaviside_phi(i, j, k);
					viscosity(i, j, k) = 1 + (m - 1)*heaviside_phi(i, j, k);
				}
				else
				{
					viscosity(i, j, k) = oil_viscosity + (water_viscosity - oil_viscosity)*heaviside_phi(i, j, k);
				}
			}
		}
		
		viscosity.FillGhostCellsFrom(thread_id, viscosity.array_for_this, true);
		
		multithreading->Sync(thread_id);
	}

	void DetermineDensityField(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& density)
	{
		FIELD_STRUCTURE_3D<T> heaviside_phi;
		heaviside_phi.Initialize(multithreading, density.grid, 2);

		HeavisideFunction(thread_id, phi, epsilon, heaviside_phi);
		
		multithreading->Sync(thread_id);
		
		GRID_ITERATION_3D(phi.partial_grids[thread_id])
		{
			if (air_water_simulation)
			{
				if (water_projection->air_bubble_rising == true)
				{
					density(i, j, k) = air_density + (water_density - air_density)*heaviside_phi(i, j, k);
				}
				if (water_projection->water_drop == true)
				{
					density(i, j, k) = water_density + (air_density - water_density)*heaviside_phi(i, j, k);
				}
			}
			if (oil_water_simulation)
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

		multithreading->Sync(thread_id);
	}

	void DeltaFunction(const int& thread_id, FIELD_STRUCTURE_3D<T>& phi, const T& epsilon, FIELD_STRUCTURE_3D<T>& delta)
	{
		int i(0), j(0);
		const int i_start(phi.i_start), j_start(phi.j_start), i_end(phi.i_end), j_end(phi.j_end);

		T one_over_epsilon = (T)1/epsilon;

		GRID_ITERATION_3D(phi.partial_grids[thread_id])
		{
			if (phi(i, j, k) < -epsilon)
			{
				delta(i, j, k) = 0;
			}
			else if (phi(i, j, k) >= -epsilon && phi(i, j, k) <= epsilon)
			{
				delta(i, j, k) = (T)0.5*one_over_epsilon + (T)0.5*one_over_epsilon*cos(PI*phi(i, j, k)*one_over_epsilon);
			}
			else
			{
				delta(i, j, k) = 0;
			}
		}

		delta.FillGhostCellsFrom(thread_id, delta.array_for_this, true);

		multithreading->Sync(thread_id);
	}

	void SourcingFromObjectList(const int& thread_id)
	{
		const T large_positive_phi = (T)5*base_grid.dx;

		BEGIN_GRID_ITERATION_3D(partial_base_grids[thread_id])
		{
			object_levelset.arr(i, j, k) = large_positive_phi;
		}
		END_GRID_ITERATION_3D;

		assert(world_discretization);

		world_discretization->UpdateAABBGrids(thread_id, base_grid);

		for (int o = 0; o < world_discretization->object_list.num_of_elements; o++)
		{
			SIMULATION_OBJECT* object = world_discretization->object_list.values[o];
			
			if (object->levelset_object == water_levelset)
			{
				continue;
			}
			const FLUID_CONTROL_OPTIONS& options(object->fluid_control_options);
			bool is_sourcing_frame = false;
			if (object->fluid_control_options.source_start_frame <= frame && frame <= object->fluid_control_options.source_stop_frame)
			{
				is_sourcing_frame = true;
			}

			if(is_sourcing_frame)
			{
				BEGIN_GRID_ITERATION_3D(object->partial_aabb_grids[thread_id])
				{
					const VT cell_center = base_grid.CellCenter(i, j, k);
					const T phi_obj(object->SignedDistance(cell_center));

					if (phi_obj <= base_grid.dx*(T)3)
					{
						if (use_water_solver)
						{
							if(options.source_water_levelset == true)
							{
								water_levelset->arr(i, j, k) = MIN(phi_obj, water_levelset->arr(i, j, k));
							}
							else if (options.delete_water_levelset == true)
							{
								water_levelset->arr(i, j, k) = MAX(-phi_obj, water_levelset->arr(i, j, k));
							}
						}	

						if (phi_obj <= (T)0)		// Inside an object
						{
							if (options.source_velocity == true)
							{
								if (use_water_solver)
								{
									/*if (abs(water_levelset->arr(i, j, k)) <= base_grid.dx*(T)2)
									{
										(*water_velocity_field)(i, j, k) = options.surface_velocity_sourcing_scale*object->Velocity(cell_center);
									}
									else
									{
										(*water_velocity_field)(i, j, k) = options.velocity_sourcing_scale*object->Velocity(cell_center);
									}*/
								}
							}
						}
					}

					if (phi_obj < object_levelset.arr(i, j, k))
					{
						object_levelset.arr(i, j, k) = phi_obj;
					}
				}
				END_GRID_ITERATION_3D;

				// Need to be updated later
			}
		}
	}

	T CFLOneTimeStep()
	{
		// Time restriction for Advection
		T u_max(0), v_max(0), w_max(0);
		
		GRID_ITERATION_3D(water_velocity_field_mac_x->grid)
		{
			if (u_max <= abs(water_velocity_field_mac_x->array_for_this(i, j, k)))
			{
				u_max = abs(water_velocity_field_mac_x->array_for_this(i, j, k));
			}
		}

		GRID_ITERATION_3D(water_velocity_field_mac_y->grid)
		{
			if (v_max <= abs(water_velocity_field_mac_y->array_for_this(i, j, k)))
			{
				v_max = abs(water_velocity_field_mac_y->array_for_this(i, j, k));
			}
		}
		
		GRID_ITERATION_3D(water_velocity_field_mac_z->grid)
		{
			if (w_max <= abs(water_velocity_field_mac_z->array_for_this(i, j, k)))
			{
				w_max = abs(water_velocity_field_mac_z->array_for_this(i, j, k));
			}
		}

		if (dimensionless_form)
		{
			c_f = u_max*water_velocity_field_mac_x->grid.one_over_dx + v_max*water_velocity_field_mac_y->grid.one_over_dy + w_max*water_velocity_field_mac_z->grid.one_over_dz;
			//c_f = Re/(T)4*POW2(u_max + v_max + w_max);
		}
		else
		{
			c_f = u_max*water_velocity_field_mac_x->grid.one_over_dx + v_max*water_velocity_field_mac_y->grid.one_over_dy + w_max*water_velocity_field_mac_z->grid.one_over_dz;
		}

		// Time restriction for Gravity
		if (is_external_force_active)
		{
			if (air_water_simulation)
			{
				g_f = sqrt(abs(gravity.y)*water_velocity_field_mac_y->grid.one_over_dy);
			}
			if (oil_water_simulation)
			{
				if (dimensionless_form)
				{
					// Note that if you use nonuniform grid, then you need to change this
					//g_f = sqrt(abs(R1/POW2(V0)*(gravity.y + ((T)1/(oil_density*MIN((T)1, water_density/oil_density)))*f))*water_velocity_field_mac_y->grid.one_over_dy);
					g_f = sqrt(R1/POW2(V0)*(abs(gravity.y) + abs((T)1/MIN(oil_density, water_density)*f))*water_velocity_field_mac_y->grid.one_over_dy);
				}
				else
				{
					g_f = sqrt(abs(gravity.y)*water_velocity_field_mac_y->grid.one_over_dy);
				}
			} 
		}
		else
		{
			g_f = 0;
		}
		
		s_f = 0;

		// Time restriction for Surface Tension
		if (air_water_simulation)
		{
			s_f = sqrt((surface_tension*water_levelset->abs_max_curvature)/(MIN(water_density, air_density)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
		}
		if (oil_water_simulation)
		{
			if (dimensionless_form)
			{
				//s_f = sqrt((surface_tension*water_levelset->abs_max_curvature)/(MIN(water_density/oil_density, 1)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
				s_f = sqrt((one_over_We*water_levelset->abs_max_curvature)/(MIN(water_density/oil_density, 1)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
				// For test
				//s_f = sqrt((We*water_levelset->abs_max_curvature)/(MIN(water_density/oil_density, 1)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
				//s_f = sqrt(((T)1/We*water_levelset->grid.one_over_dx)/(MIN(water_density/oil_density, 1)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
				//s_f = sqrt(((POW2(R1)*V0)/We*water_levelset->abs_max_curvature)/(MIN(water_density, oil_density)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
			}
			else
			{
				s_f = sqrt((surface_tension*water_levelset->abs_max_curvature)/(MIN(water_density, oil_density)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
			}	
			
		}
		
		// Time restriction for Viscosity Term
		if (air_water_simulation)
		{
			v_f = max(water_viscosity/water_density, air_viscosity/air_density)*((T)2*base_grid.one_over_dx2 + (T)2*base_grid.one_over_dy2 + (T)2*base_grid.one_over_dz2);
		}
		if (oil_water_simulation)
		{
			if (viscosity_solver->semi_implicit_approach)
			{
				v_f = 0;
			}
			else
			{
				v_f = max(water_viscosity/water_density, oil_viscosity/oil_density)*((T)2*base_grid.one_over_dx2 + (T)2*base_grid.one_over_dy2 + (T)2*base_grid.one_over_dz2);
			}
		}
		
		T cfl = (T)0.5*(c_f + v_f + sqrt(POW2(c_f + v_f) + (T)4*POW2(g_f) + (T)4*POW2(s_f)));
		
		T one_over_cfl = (T)1/cfl;
		
		return one_over_cfl;
	}
	
	T CFLOneTimeStep(const int& thread_id)
	{
		// Time restriction for Advection
		T u_max(0), v_max(0), w_max(0);
		
		BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_x->partial_grids[thread_id])
		{
			if (u_max <= abs(water_velocity_field_mac_x->array_for_this(i, j, k)))
			{
				u_max = abs(water_velocity_field_mac_x->array_for_this(i, j, k));
			}
		}
		END_GRID_ITERATION_MAX_3D(u_max);

		BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_y->partial_grids[thread_id])
		{
			if (v_max <= abs(water_velocity_field_mac_y->array_for_this(i, j, k)))
			{
				v_max = abs(water_velocity_field_mac_y->array_for_this(i, j, k));
			}
		}
		END_GRID_ITERATION_MAX_3D(v_max);

		BEGIN_GRID_ITERATION_3D(water_velocity_field_mac_z->partial_grids[thread_id])
		{
			if (w_max <= abs(water_velocity_field_mac_z->array_for_this(i, j, k)))
			{
				w_max = abs(water_velocity_field_mac_z->array_for_this(i, j, k));
			}
		}
		END_GRID_ITERATION_MAX_3D(w_max);

		c_f = u_max*water_velocity_field_mac_x->grid.one_over_dx + v_max*water_velocity_field_mac_y->grid.one_over_dy + w_max*water_velocity_field_mac_z->grid.one_over_dz;

		// Time restriction for Gravity
		if (is_external_force_active)
		{
			if (air_water_simulation)
			{
				g_f = sqrt(abs(gravity.y)*water_velocity_field_mac_y->grid.one_over_dy);
			}
			if (oil_water_simulation)
			{
				if (dimensionless_form)
				{
					// Note that if you use nonuniform grid, then you need to change this
					g_f = sqrt(abs(R1/POW2(V0)*(gravity.y + ((T)1/(oil_density*MIN((T)1, water_density/oil_density)))*f))*water_velocity_field_mac_y->grid.one_over_dy);
					//g_f = sqrt((abs(gravity.y) + abs(f/oil_density))*water_velocity_field_mac_y->grid.one_over_dy);
				}
				else
				{
					g_f = sqrt(abs(gravity.y)*water_velocity_field_mac_y->grid.one_over_dy);
				}
			} 
		}
		else
		{
			g_f = 0;
		}
		
		// Time restriction for Surface Tension
		if (air_water_simulation)
		{
			s_f = sqrt((surface_tension*water_levelset->abs_max_curvature)/(MIN(water_density, air_density)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
		}
		if (oil_water_simulation)
		{
			if (dimensionless_form)
			{
				s_f = sqrt((R1/We*water_levelset->abs_max_curvature)/(MIN(water_density/oil_density, 1)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
				//s_f = sqrt(((POW2(R1)*V0)/We*water_levelset->abs_max_curvature)/(MIN(water_density, oil_density)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
			}
			else
			{
				s_f = sqrt((surface_tension*water_levelset->abs_max_curvature)/(MIN(water_density, oil_density)*POW2(MIN3(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy, water_velocity_field_mac_z->grid.dz))));
			}	
			
		}
		
		// Time restriction for Viscosity Term
		if (air_water_simulation)
		{
			v_f = max(water_viscosity/water_density, air_viscosity/air_density)*((T)2*base_grid.one_over_dx2 + (T)2*base_grid.one_over_dy2 + (T)2*base_grid.one_over_dz2);
		}
		if (oil_water_simulation)
		{
			if (viscosity_solver->semi_implicit_approach)
			{
				v_f = 0;
			}
			else
			{
				v_f = max(water_viscosity/water_density, oil_viscosity/oil_density)*((T)2*base_grid.one_over_dx2 + (T)2*base_grid.one_over_dy2 + (T)2*base_grid.one_over_dz2);
			}
		}
		
		T cfl = (T)0.5*(c_f + v_f + sqrt(POW2(c_f + v_f) + (T)4*POW2(g_f) + (T)4*POW2(s_f)));
		
		T one_over_cfl = (T)1/cfl;
		
		return one_over_cfl;
	}
		
};

class FLUID_OBJECT : public	SIMULATION_OBJECT
{
public: // Essential data
	FIELD_STRUCTURE_3D<VT>* velocity_field;
	LEVELSET_3D* water_levelset;

public: // Constructors and Destructor
	FLUID_OBJECT(void)
		: velocity_field(0), water_levelset(0)
	{}

	~FLUID_OBJECT(void)
	{}

public: // Initialization Function
	void Initialize(FIELD_STRUCTURE_3D<VT>* velocity_input, LEVELSET_3D* levelset_input)
	{
		velocity_field = velocity_input;
		water_levelset = levelset_input;
		levelset_object = levelset_object;
	}

public: // Member Functions
	const T SignedDistance(const VT& position) const
	{
		return water_levelset->SignedDistance(position);
	}

	const bool Inside(const VT& position) const
	{
		if(water_levelset->SignedDistance(position) <= 0) return true;
		else return false;
	}

	const VT ClosestPoint(const VT& position) const
	{
		return position - (water_levelset->UnitNormal(position)*water_levelset->SignedDistance(position));
	}

	const VT Velocity(const VT& position) const
	{
		return (*velocity_field)(position);
	}

	const void Normal(const VT& position, VT& normal) const
	{
		normal = water_levelset->Normal(position);
	}

	const VT Normal(const VT& position) const
	{
		return water_levelset->Normal(position);
	}

	const void UnitNormal(const VT& position, VT& normal) const
	{
		water_levelset->UnitNormal(position, normal);
	}

	const VT UnitNormal(const VT& position) const
	{
		return water_levelset->UnitNormal(position);
	}
};


