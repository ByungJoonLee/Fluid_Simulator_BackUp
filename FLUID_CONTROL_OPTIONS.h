#pragma once

#include "COMMON_DEFINITIONS.h"

class FLUID_CONTROL_OPTIONS
{
public: // Options
	bool	source_density;
	bool	source_temperature;
	bool	source_velocity;
	bool	source_water_levelset;
	bool	delete_water_levelset;
	bool	source_vortex_particles;
	bool	source_one_time;

	int		source_start_frame;
	int		source_stop_frame;

	bool	use_normal_as_velocity;					// Currently, implemented just for a gas simulation
	T		normal_direction_velocity_coefficient;

	T		velocity_sourcing_scale;
	T		surface_velocity_sourcing_scale;
	T		sourcing_temperature;

public: // Constructor and Destructor
	FLUID_CONTROL_OPTIONS(void)
		: source_density(false)
		, source_temperature(false)
		, source_velocity(false)
		, source_water_levelset(false)
		, delete_water_levelset(false)
		, source_vortex_particles(false)
		, source_start_frame(-numeric_limits<int>::max())
		, source_stop_frame(numeric_limits<int>::max())
		, velocity_sourcing_scale((T)1)
		, surface_velocity_sourcing_scale((T)1)
		, sourcing_temperature((T)3000)
		, use_normal_as_velocity(false)
		, normal_direction_velocity_coefficient((T)1.0)
		, source_one_time(false)
	{}

public: // Initialization Function
	void InitializeFromScript(SCRIPT_BLOCK& outer_block)
	{
		source_density = outer_block.GetBoolean("source_density", source_density);
		source_temperature = outer_block.GetBoolean("source_temperature", (bool)false);
		source_velocity = outer_block.GetBoolean("source_velocity", source_velocity);
		use_normal_as_velocity = outer_block.GetBoolean("use_normal_as_velocity", false);
		normal_direction_velocity_coefficient = outer_block.GetFloat("normal_direction_velocity_coefficient", (T)1.0);
		source_water_levelset = outer_block.GetBoolean("source_levelset", (bool)false);
		delete_water_levelset = outer_block.GetBoolean("delete_water_levelset", (bool)false);
		source_vortex_particles = outer_block.GetBoolean("source_vortex_particles", (bool)false);
		source_one_time = outer_block.GetBoolean("source_one_time", (bool)false);
		source_start_frame = outer_block.GetInteger("source_start_frame", source_start_frame);
		source_stop_frame = outer_block.GetInteger("source_stop_frame", source_stop_frame);

		velocity_sourcing_scale = outer_block.GetFloat("velocity_sourcing_scale", velocity_sourcing_scale);
		surface_velocity_sourcing_scale = outer_block.GetFloat("surface_velocity_sourcing_scale", surface_velocity_sourcing_scale);
		sourcing_temperature = outer_block.GetFloat("sourcing_temperature", (T)3000);
	}
};


