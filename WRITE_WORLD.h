#pragma once

#include "WRITE_OBJECT_BASE.h"
#include "COMMON_DEFINITIONS.h"
#include "GRID_STRUCTURE_3D.h"
#include "SIMULATION_WORLD.h"
#include "WRITE_LEVELSET.h"
#include <vector>
#include <string>

#define WRITE_OPTION_FILE_NAME "write.wr"

class SIMULATION_WORLD;
class EULERIAN_FLUID_SOLVER_3D;
class MULTITHREADING;

class WRITE_WORLD
{
public: // Essential Data
	SIMULATION_WORLD*			simulation_world;
	MULTITHREADING*				multithreading;
	GRID_STRUCTURE_3D*			world_grid;

	string						app_dir;
	string						script_dir;

	vector<WRITE_OBJECT_BASE*>	all_objects;

public: // Constructor and Destructor
	WRITE_WORLD(void)
		: simulation_world(0), multithreading(0), world_grid(0), app_dir(string()), script_dir(string())
	{}

	~WRITE_WORLD(void)
	{
		DeleteAllObjects();
	}

public: // Initialization Functions
	void Initialize(SIMULATION_WORLD* simulation_world_input)
	{
		simulation_world = simulation_world_input;
		if (simulation_world == 0)
		{
			return; 
		}

		DeleteAllObjects();

		script_dir = simulation_world->script_abs_dir;
		multithreading = (simulation_world->multithreading);
		world_grid = &(simulation_world->world_discretization.world_grid);

		if (simulation_world->use_eulerian_solver)
		{
			EULERIAN_FLUID_SOLVER_3D* eulerian_solver = &(simulation_world->eulerian_solver);

			// Water levelset
			if (eulerian_solver->use_water_solver)
			{
				all_objects.push_back(new WRITE_LEVELSET("WATER_LEVELSET", eulerian_solver->water_levelset));
			}
		}
		
		// Will be updated later
		LoadWriteOptions();
	}

	void Write()
	{
		for (vector<WRITE_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
		{
			(*it)->Write(simulation_world->num_current_frame);
		}
	}
	
	void Key(unsigned char c)
	{
		switch (c)
		{
		case '1' : 
			LoadWriteOptions();
			break;
		}
	}
	
	void KeyWithAlt(unsigned char c)
	{}

	void SpecialKey(int key)
	{}

	bool LoadWriteOptions()
	{
		string file_full_path("");

		if (!script_dir.empty())
		{
			boost::filesystem::path p(script_dir);
			if (p.is_relative())
			{
				p = boost::filesystem::absolute(p);
				script_dir = p.string();
			}
			file_full_path = script_dir + WRITE_OPTION_FILE_NAME;
			if (!boost::filesystem::exists(boost::filesystem::path(file_full_path)))
			{
				return false;
			}
		}
		else
		{
			return false;
		}

		SCRIPT_READER script_reader;
		bool is_succeed = script_reader.Initialize(file_full_path.c_str());
		if (is_succeed)
		{
			SCRIPT_BLOCK script_block = script_reader.FindBlock("WRITE_WORLD");
			for (vector<WRITE_OBJECT_BASE*>::iterator it  = all_objects.begin(); it != all_objects.end() ; ++it)
			{
				(*it)->ParseWriteOptions(&script_block, script_dir.c_str());
			}
		}
		else
		{
			return false;
		}
		
		return true;
	}

	void DeleteAllObjects()
	{
		for (vector<WRITE_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
		{
			delete (*it);
		}

		vector<WRITE_OBJECT_BASE*>().swap(all_objects);
	}
};
	
		




