#pragma once

#include "OPENGL_SOLVER_BASE.h"
#include "EULERIAN_FLUID_SOLVER_3D.h"
#include "OPENGL_CG_VOLUME.h"
#include "OPENGL_LEVELSET.h"
#include "OPENGL_SCALARFIELD.h"
#include "OPENGL_VECTORFIELD.h"
#include <vector>

class EULERIAN_FLUID_SOLVER_3D;
class OPENGL_LEVELSET;
class OPENGL_CG_VOLUME;
class OPENGL_SCALARFIELD;

class OPENGL_EULERIAN_FLUID_SOLVER : public OPENGL_SOLVER_BASE
{
public: // Essential Data
	OPENGL_DRIVER*				driver;
	EULERIAN_FLUID_SOLVER_3D*	eulerian_solver;
	
	OPENGL_LEVELSET*			water_levelset;
	OPENGL_LEVELSET*			boundary_levelset;
	OPENGL_LEVELSET*			water_levelset_upsampled;
	
	// For Poisson Test
	OPENGL_SCALARFIELD*			pressure;

	OPENGL_VECTORFIELD*			velocity_field_x;
	OPENGL_VECTORFIELD*			velocity_field_y;
	OPENGL_VECTORFIELD*			velocity_field_z;
		
	OPENGL_CG_VOLUME*			temperature;			// Both for ordinary & upsampled
	OPENGL_CG_VOLUME*			rigid_body_levelset;
	OPENGL_CG_VOLUME*			smoke_density;			// Both for ordinary & upsampled

public: // Constructor and Destructor
	OPENGL_EULERIAN_FLUID_SOLVER(OPENGL_DRIVER* driver_input, EULERIAN_FLUID_SOLVER_3D* eulerian_object, MULTITHREADING* multithreading_input, float grid_scale)
		: driver(driver_input), eulerian_solver(eulerian_object), water_levelset(0), pressure(0), boundary_levelset(0), water_levelset_upsampled(0), smoke_density(0), velocity_field_x(0), velocity_field_y(0), velocity_field_z(0)
	{
		/*
		// water levelset
		if(eulerian_solver_->num_hybrid_water_upsampling_ == 1 && eulerian_solver_->num_water_levelset_upsampling_ == 1)
		water_levelset_ = new OPENGL_LEVELSET("EULERIAN_WATER_LEVELSET", driver_, &(eulerian_solver_->water_levelset_), multithreading, grid_scale);
		else
		water_levelset_ = new OPENGL_LEVELSET("EULERIAN_WATER_LEVELSET", driver_, &(eulerian_solver_->water_levelset_upsampled_), multithreading, grid_scale);
		*/	

		bool use_levelset_upsampling = false;
		
		if (eulerian_solver->use_water_solver)
		{
			water_levelset = new OPENGL_LEVELSET("EULERIAN_WATER_LEVELSET", driver, eulerian_solver->water_levelset, multithreading_input, grid_scale);
			AddObject(water_levelset);
			
			boundary_levelset = new OPENGL_LEVELSET("BOUNDARY_LEVELSET", driver, eulerian_solver->boundary_levelset, multithreading_input, grid_scale);
			AddObject(boundary_levelset);
			boundary_levelset->SetDrawType(OPENGL_LEVELSET::LEVELSET_DRAW_SOLID);
			boundary_levelset->is_boundary = true;
			boundary_levelset->is_levelset = false;

			pressure = new OPENGL_SCALARFIELD("PRESSURE", driver, &eulerian_solver->water_projection->pressure_field);
			AddObject(pressure);
			pressure->SetDrawType(OPENGL_SCALARFIELD::SCALARFIELD_DRAW_SHOW_Y);

			/*velocity_field_x = new OPENGL_VECTORFIELD("VELOCITY_FIELD_X", driver, eulerian_solver->water_velocity_field, eulerian_solver->water_velocity_field_mac_x);
			AddObject(velocity_field_x);
			velocity_field_x->SetDrawType(OPENGL_VECTORFIELD::VECTORFIELD_DRAW_SHOW_X);
			velocity_field_y = new OPENGL_VECTORFIELD("VELOCITY_FIELD_Y", driver, eulerian_solver->water_velocity_field, eulerian_solver->water_velocity_field_mac_y);
			AddObject(velocity_field_y);
			velocity_field_y->SetDrawType(OPENGL_VECTORFIELD::VECTORFIELD_DRAW_SHOW_Y);
			velocity_field_z = new OPENGL_VECTORFIELD("VELOCITY_FIELD_Z", driver, eulerian_solver->water_velocity_field, eulerian_solver->water_velocity_field_mac_z);
			AddObject(velocity_field_z);
			velocity_field_z->SetDrawType(OPENGL_VECTORFIELD::VECTORFIELD_DRAW_SHOW_Z);*/

			// Water levelset scalar field
			AddObject(new OPENGL_SCALARFIELD("WATER_LEVELSET_FIELD", driver, &(eulerian_solver->water_levelset->signed_distance_field), SCALARFIELD_LEVELSET_MODE));
		}
		// Need to be fixed some
	}

	~OPENGL_EULERIAN_FLUID_SOLVER(void)
	{}

public: // Member Functions
	virtual void LoadOpenGLSettings(SCRIPT_BLOCK* script_block)
	{
		// Water Levelset
		SCRIPT_BLOCK* solver_block = script_block->SearchBlock("OPENGL_EULERIAN_FLUID_SOLVER");
		if (!solver_block)
		{
			return;
		}

		SCRIPT_BLOCK* water_levelset_block = solver_block->SearchBlock("WATER_LEVELSET");
		if (water_levelset && water_levelset_block)
		{
			OPENGL_MATERIAL material;
			const char* solid_color_name = water_levelset_block->GetString("solid_material_name", 0);
			if (solid_color_name)
			{
				material.SetMaterialByName(solid_color_name);
			}
			else
			{
				VT ambient = water_levelset_block->GetVector3("solid_ambient_color");
				VT diffuse = water_levelset_block->GetVector3("solid_diffuse_color");
				VT specular = water_levelset_block->GetVector3("solid_specular_color");
				GLfloat shininess = water_levelset_block->GetFloat("solid_shininess");

				material.ambient_color.Set(ambient.x, ambient.y, ambient.z);
				material.diffuse_color.Set(diffuse.x, diffuse.y, diffuse.z);
				material.specular_color.Set(specular.x, specular.y, specular.z);
				material.shininess = shininess;
			}

			water_levelset->SetMaterial(material);

			if (water_levelset_upsampled)
			{
				water_levelset_upsampled->SetMaterial(material);
			}
			
			// Display given variables
			cout << "-------------- OPENGL EULERIAN FLUID SOLVER - LEVELSET --------------" << endl;
			cout << "solid material name : " << solid_color_name << endl;
		}

		SCRIPT_BLOCK* boundary_levelset_block = solver_block->SearchBlock("BOUNDARY_LEVELSET");
		if (boundary_levelset && boundary_levelset_block)
		{
			OPENGL_MATERIAL material;
			const char* solid_color_name = boundary_levelset_block->GetString("solid_material_name", 0);
			if (solid_color_name)
			{
				material.SetMaterialByName(solid_color_name);
			}
			else
			{
				VT ambient = boundary_levelset_block->GetVector3("solid_ambient_color");
				VT diffuse = boundary_levelset_block->GetVector3("solid_diffuse_color");
				VT specular = boundary_levelset_block->GetVector3("solid_specular_color");
				GLfloat shininess = boundary_levelset_block->GetFloat("solid_shininess");

				material.ambient_color.Set(ambient.x, ambient.y, ambient.z);
				material.diffuse_color.Set(diffuse.x, diffuse.y, diffuse.z);
				material.specular_color.Set(specular.x, specular.y, specular.z);
				material.shininess = shininess;
			}

			boundary_levelset->SetMaterial(material);
						
			// Display given variables
			cout << "-------------- OPENGL EULERIAN FLUID SOLVER - BOUNDARY LEVELSET --------------" << endl;
			cout << "solid material name : " << solid_color_name << endl;
		}
		// Smoke densith field will be added later
	}
};


		