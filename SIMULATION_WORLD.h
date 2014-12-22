#pragma once

#include "stdafx.h"
#include "SIMULATION_OBJECT.h"
#include "EULERIAN_FLUID_SOLVER_3D.h"
#include "WORLD_DISCRETIZATION_3D.h"
#include "STATIC_OBJECT.h"
#include "RIGID_OBJECT.h"
#include "SPHERE.h"
#include "CYLINDER.h"
#include "LOG.h"
#include <ctype.h>
#include <boost/chrono.hpp>

class SIMULATION_WORLD
{
public:
	// Primitive Solvers which are updated
	EULERIAN_FLUID_SOLVER_3D				eulerian_solver;

public: 
	DYNAMIC_ARRAY<SIMULATION_OBJECT*>&		object_list;
	WORLD_DISCRETIZATION_3D					world_discretization;

public: // Options for Simulation
	bool									air_water_simulation;
	bool									oil_water_simulation;

public: // Options for Oil Water Simulation
	bool									is_vertical, is_parallel;

public: // Wave length and Wave number for Oil Water Simulation
	T										a, R2, A_0, alpha;

public: // Multithreading 
	MULTITHREADING*							multithreading;

public:
	// NETWORKING networking;
	std::string								script_abs_dir;
	std::string								app_abs_dir;
	std::string								script_file_path;
	std::string								script_base_name;

public:
	// Simulation properties
	T										frame_rate;
	T										dt, accu_dt, max_dt;
	T										CFL;
	int										num_substeps;
	int										num_current_frame;

	// Simulation Options
	bool									use_eulerian_solver;
	bool									use_sph_solver;
	bool									use_solid_solver;
	bool									use_sph_levelset;
	bool									use_deform_solver;
	bool									use_mpm_solver;

	// Option
	int										last_frame;
	bool									auto_run;
	bool									is_polygonize_levelset;

	// PLUGIN_FUNC plugin_initialize;
	// PLUGIN_FUNC plugin_advance;
	// PLUGIN_FUNC plugin_finalize;

	// FBX view
	bool									use_fbx_camera;
	string									fbx_file;

	T										grid_scale;
	VT										eulerian_color;

public: // Constructor and Destructor
	SIMULATION_WORLD(void)
		: use_eulerian_solver(false), use_sph_solver(false), use_solid_solver(false), use_deform_solver(false), use_sph_levelset(false),
		num_current_frame(0), last_frame(numeric_limits<int>::max()), eulerian_solver(num_current_frame), multithreading(0), object_list(world_discretization.object_list),
		dt((T)0), accu_dt((T)0), max_dt((T)0), air_water_simulation(false), oil_water_simulation(false), is_vertical(false), is_parallel(false)
	{}

	~SIMULATION_WORLD(void)
	{
		DELETE_POINTER(multithreading);
	}

public: // Member Functions
	void InitializeFormScript(const char* script)
	{
		DELETE_POINTER(multithreading);
		multithreading = new MULTITHREADING;

		SetCurrentDirectory(script_abs_dir.c_str());
		
		SCRIPT_READER script_reader(script);

		num_current_frame = 0;

		SCRIPT_BLOCK script_block_for_this = script_reader.FindBlock("SIMULATION_WORLD");

		dt = script_block_for_this.GetFloat("dt", (T)0.01);
		max_dt = script_block_for_this.GetFloat("max_dt", (T)100);
		CFL = script_block_for_this.GetFloat("CFL", (T)0);
		frame_rate = script_block_for_this.GetFloat("frame_rate", (T)24);
		
		num_substeps = script_block_for_this.GetInteger("num_substeps", (int)1);

		use_eulerian_solver = script_block_for_this.GetBoolean("use_eulerian_solver");
		use_sph_solver		= script_block_for_this.GetBoolean("use_sph_solver");
		use_solid_solver	= script_block_for_this.GetBoolean("use_solid_solver");
		use_sph_levelset	= script_block_for_this.GetBoolean("use_sph_levelset");
		use_deform_solver   = script_block_for_this.GetBoolean("use_deform_solver");
		use_mpm_solver		= script_block_for_this.GetBoolean("use_mpm_solver");

		// Simulation Options
		air_water_simulation = script_block_for_this.GetBoolean("air_water_simulation", false);
		oil_water_simulation = script_block_for_this.GetBoolean("oil_water_simulation", false);
		
		if (oil_water_simulation)
		{
			a = script_block_for_this.GetFloat("a", (T)1);
			R2 = script_block_for_this.GetFloat("R2", (T)1);
			A_0 = script_block_for_this.GetFloat("A_0", (T)1);
			A_0 = A_0/(R2/a);
			alpha = script_block_for_this.GetFloat("alpha", (T)1);
		}

		multithreading->Initialize(script_block_for_this.GetInteger("number_of_threads"));

		last_frame = script_block_for_this.GetInteger("last_frame");
		auto_run = script_block_for_this.GetBoolean("auto_run");

		// Initialize LOG
		LOG::InitializeFromScriptBlock(multithreading, script_block_for_this.FindBlock("LOG"), script_abs_dir);

		// Display given variables
		cout << "--------------SIMULATION WORLD VARIABLES--------------" << endl;
		cout << "dt : " << dt << endl;
		cout << "max dt : " << max_dt << endl;
		cout << "CFL: " << CFL << endl;
		cout << "frame rate: " << frame_rate << endl;
		cout << "number of substeps: " << num_substeps << endl;
		
		if (use_eulerian_solver == true)
		{
			cout << "use eulerian solver: " << "true" << endl;
		}
		else
		{
			cout << "use eulerian solver: " << "false" << endl;
		}
		
		if (use_sph_levelset == true)
		{
			cout << "use sph solver: " << "true" << endl;
		}
		else
		{
			cout << "use sph solver: " << "false" << endl;
		}
		
		if (use_solid_solver == true)
		{
			cout << "use solid solver: " << "true" << endl;
		}
		else
		{
			cout << "use solid solver: " << "false" << endl;
		}

		if (use_sph_levelset == true)
		{
			cout << "use sph levelset: " << "true" << endl;
		}
		else
		{
			cout << "use sph levelset: " << "false" << endl;
		}

		if (use_deform_solver == true)
		{
			cout << "use deform solver: " << "true" << endl;
		}
		else
		{
			cout << "use deform solver: " << "false" << endl;
		}
		
		if (use_mpm_solver)
		{
			cout << "use mpm solver: " << "true" << endl;
		}
		else
		{
			cout << "use mpm solver: " << "false" << endl;
		}
		
		cout << "number of threads: " << multithreading->num_threads << endl;
		cout << "last frame: " << last_frame << endl;
		cout << "auto run: " << auto_run << endl;
		
		// Initialize world discretization
		if (air_water_simulation)
		{
			world_discretization.air_water_simulation = true;
			world_discretization.oil_water_simulation = false;
			world_discretization.Initialize(multithreading, script_block_for_this.FindBlock("WORLD_DISCRETIZATION_AIR_WATER"));
		}
		if (oil_water_simulation)
		{
			world_discretization.air_water_simulation = false;
			world_discretization.oil_water_simulation = true;
			world_discretization.Initialize(multithreading, script_block_for_this.FindBlock("WORLD_DISCRETIZATION_OIL_WATER"));
			is_vertical = world_discretization.is_vertical;
			is_parallel = world_discretization.is_parallel;
		}

		InitializeDynamicsSolversFromScript(script_reader);
		InitializeObjectListFromScript(script_reader);
		
		if (object_list.num_of_elements > 0)
		{
			multithreading->RunThreads(&EULERIAN_FLUID_SOLVER_3D::SourcingFromObjectList, &eulerian_solver);
		}

		SCRIPT_BLOCK viewer_block = script_block_for_this.FindBlock("VIEWER");
		use_fbx_camera = viewer_block.GetBoolean("use_fbx_camera");
		fbx_file = viewer_block.GetString("fbx_file");

		grid_scale = viewer_block.GetFloat("grid_scale", 1);
		T darkness = 0.4f;

		eulerian_color = viewer_block.GetVector3("eulerian_color", VT((T)0.27*darkness, (T)0.38*darkness, (T)0.49*darkness));

		SetCurrentDirectory(app_abs_dir.c_str());
	}
		
	void InitializeDynamicsSolversFromScript(SCRIPT_READER& script_reader)
	{
		if (use_eulerian_solver)
		{
			eulerian_solver.world_discretization = &world_discretization;
			if (air_water_simulation)
			{
				eulerian_solver.air_water_simulation = true;
				eulerian_solver.oil_water_simulation = false;
				eulerian_solver.InitializeFromScriptBlock(multithreading, world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM_AIR_WATER"));	
			}
			if (oil_water_simulation)
			{
				eulerian_solver.air_water_simulation = false;
				eulerian_solver.oil_water_simulation = true;
				eulerian_solver.a = a;
				eulerian_solver.R2 = R2;
				eulerian_solver.InitializeFromScriptBlock(multithreading, world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM_OIL_WATER"));	
			}
		}
	}
	
	void AdvanceOneTimeStepThread(const int& thread_id, const T& dt_input)
	{
		world_discretization.Update(thread_id);

		if (use_eulerian_solver)
		{
			eulerian_solver.AdvanceOneTimeStepThread(thread_id, dt_input);
		}
	}

	void AdvanceOneFrameThread(const int& thread_id)
	{
		if(CFL == (T)0 || frame_rate == (T)0)
		{
			for(int substep = 0; substep < num_substeps; substep++)
			{
				LOG::Begin(thread_id, "Substep %d", substep);

				AdvanceOneTimeStepThread(thread_id, dt);

				LOG::End(thread_id);
			}
		}
		else
		{
			BEGIN_HEAD_THREAD_WORK
			{
				dt = DetermineTimeStep()*CFL;
				accu_dt += dt;
				cout << "Max x-velocity   : " << eulerian_solver.water_projection->max_velocity_x << endl;
				cout << "Max y-velocity   : " << eulerian_solver.water_projection->max_velocity_y << endl;
				cout << "Max z-velocity   : " << eulerian_solver.water_projection->max_velocity_z << endl;
				cout << "Time Step        : " << dt << endl;
				cout << "Accumulated Time : " << accu_dt << endl;
				cout << "-------CFL Time Step-------" << endl;
				cout << "c_f              : " << eulerian_solver.c_f << endl;
				cout << "g_f              : " << eulerian_solver.g_f << endl;
				cout << "s_f              : " << eulerian_solver.s_f << endl;
				cout << "v_f              : " << eulerian_solver.v_f << endl;
				
				if (oil_water_simulation)
				{
					if (accu_dt > max_dt)
					{
						cout << "Time step reaches at max!:)" << endl;
						exit(0);
					}
				}
			}
			END_HEAD_THREAD_WORK;

			AdvanceOneTimeStepThread(thread_id, dt);

			multithreading->Sync(thread_id);
		}
	}

	void AdvanceOneFrame()
	{
		LOG::Begin("Frame %d", num_current_frame);

		eulerian_solver.water_levelset->ComputeCurvaturesThreaded();
		multithreading->RunThreads(&SIMULATION_WORLD::AdvanceOneFrameThread, this);

		LOG::End();

		num_current_frame++;
	}

	T DetermineTimeStep()
	{
		//T dt_for_this = (T)1/frame_rate;
		T dt_for_this(0);

		if (use_eulerian_solver)
		{
			//dt_for_this = MIN(dt_for_this, eulerian_solver.CFLOneTimeStep());
			dt_for_this = eulerian_solver.CFLOneTimeStep();
		}

		return dt_for_this;
	}

	T DetermineTimeStep(const int& thread_id)
	{
		T dt_for_this = (T)1/frame_rate;

		if (use_eulerian_solver)
		{
			dt_for_this = MIN(dt_for_this, eulerian_solver.CFLOneTimeStep(thread_id));
		}

		return dt_for_this;
	}

	void Save()
	{
		SetCurrentDirectory(script_abs_dir.c_str());
		SetCurrentDirectory(app_abs_dir.c_str());
	}

	void Load()
	{
		SetCurrentDirectory(script_abs_dir.c_str());
		SetCurrentDirectory(app_abs_dir.c_str());
	}
	
	void InitializeObjectListFromScript(SCRIPT_READER& script_reader);
};

	



