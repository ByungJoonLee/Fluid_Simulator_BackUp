#pragma once

#include "OPENGL_LIGHT_MANAGER.h"
#include "OPENGL_SOLVER_BASE.h"
#include "SIMULATION_WORLD.h"
#include "OPENGL_EULERIAN_FLUID_SOLVER.h"
#include "OPENGL_LEVELSET.h"
#include "OPENGL_2D_TEXT.h"
#include <list>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <map>

static bool FileExists(const char* filename)
{
	if (FILE* file = fopen(filename, "r"))
	{
		fclose(file);
		return true;
	}
	return false;
}

class SIMULATION_WORLD;
class OPENGL_SIMULATION_BOX : public OPENGL_OBJECT_BASE
{
public: // Constructor and Destructor
	OPENGL_SIMULATION_BOX(OPENGL_DRIVER* driver)
		: OPENGL_OBJECT_BASE("Simulation Box", driver)
	{}

	~OPENGL_SIMULATION_BOX(void)
	{}

public: // Virtual Functions
	virtual int NextDrawType()
	{
		return 0;
	}

	virtual int PreviousDrawType()
	{
		return 0;
	}

	virtual void Render()
	{
		GetDriver()->SetDefaultRenderStatesLineMode();
		GetDriver()->DrawWireBox(GetLength().x, GetLength().y, GetLength().z);
	}
};

class OPENGL_GROUND : public OPENGL_OBJECT_BASE
{
public: // Essential Data
	float width;

public: // Constructor and Destructor
	OPENGL_GROUND(OPENGL_DRIVER* driver)
		: OPENGL_OBJECT_BASE("Ground", driver), width(1.2)
	{}

	~OPENGL_GROUND(void)
	{}

public: // Member functions
	void CalWidth(GRID_STRUCTURE_3D& grid)
	{
		float diff_x = abs(grid.x_max - grid.x_min);
		float diff_y = abs(grid.y_max - grid.y_min);
		float diff_z = abs(grid.z_max - grid.z_min);
		float diff_max = (std::max)(diff_x, diff_y);
		diff_max = (std::max)(diff_max, diff_z);
		diff_max += (diff_max / 4.0f);
		width = diff_max;
	}
	
public: // Virtual Functions
	virtual int NextDrawType()
	{
		return 0;
	}

	virtual int PreviousDrawType()
	{
		return 0;
	}

	virtual void Render()
	{
		GetDriver()->SetDefaultRenderStatesLineMode();
		GetDriver()->DrawWireSquarePlane(width, 20);
	}
};

class OPENGL_EULERIAN_FLUID_SOLVER;

class OPENGL_WORLD
{
public: // Essential Data
	SIMULATION_WORLD*				simulation_world;
	MULTITHREADING*					multithreading;
	OPENGL_DRIVER*					driver;
	float							grid_scale;

	vector<OPENGL_OBJECT_BASE*>		all_objects;
	vector<OPENGL_SOLVER_BASE*>		all_solvers;

	OPENGL_EULERIAN_FLUID_SOLVER*	opengl_fluid_solver;

	OPENGL_LIGHT_MANAGER*			light_manager;
	OPENGL_GROUND*					gl_ground;
	OPENGL_SIMULATION_BOX*			gl_simulation_box;

	OPENGL_OBJECT_BASE*				selected_object;
	int								selected_object_index;
	OPENGL_2D_MENU_TEXT*			menu_object;
	OPENGL_2D_MENU_TEXT*			menu_draw_type;
	bool							is_show_menu;
	bool							is_object_menu;

	OPENGL_2D_BLOCK_TEXT*			bottom_text;
	OPENGL_2D_BLOCK_TEXT*			right_text;

	OPENGL_COLOR					bg_upper_color;
	OPENGL_COLOR					bg_lower_color;

	string							app_abs_dir;
	string							script_abs_dir;
	bool							light_draw_for_debug;

	bool							draw_center_axis;
	bool							draw_corner_axis;
	bool							draw_simulation_box;
	bool							draw_ground;
	bool							draw_info;

	bool							draw_levelset;
	bool							draw_velocity;
	bool							draw_pressure;
	bool							draw_boundary;

	bool							draw_velocity_x;
	bool							draw_velocity_y;
	bool							draw_velocity_z;

	// For periodicity
	bool							is_periodic;
	int								periodic_num_x, periodic_num_y;

	// For report
	int								min_num_water_triangles;
	int								max_num_water_triangles;
	int								ave_num_water_triangles;
	int								ave_num_water_triangles_sum;

	int								render_count;

	string							info_text;
	string							info_text_def;

public: // Constructor and Destructor
	OPENGL_WORLD(void)
		: simulation_world(0)
		, multithreading(0)
		, grid_scale(0)
		, selected_object(0)
		, selected_object_index(0)
		, menu_object(new OPENGL_2D_MENU_TEXT(FONT_VERDANA, 13))
		, menu_draw_type(new OPENGL_2D_MENU_TEXT(FONT_VERDANA, 13))
		, is_show_menu(false)
		, is_object_menu(true)
		, driver(new OPENGL_DRIVER)
		, bg_upper_color(1.0f, 1.0f, 1.0f)
		, bg_lower_color(1.0f, 1.0f, 1.0f)
		, app_abs_dir(string())
		, script_abs_dir(string())
		, light_draw_for_debug(false)
		, draw_center_axis(false)
		, draw_corner_axis(true)
		, draw_simulation_box(true)
		, draw_ground(true)
		, draw_info(true)
		, draw_levelset(true)
		, draw_pressure(false)
		, draw_velocity(false)
		, draw_boundary(true)
		, draw_velocity_x(false)
		, draw_velocity_y(false)
		, draw_velocity_z(false)
		, is_periodic(false)
		, periodic_num_x(1)
		, periodic_num_y(1)
		, opengl_fluid_solver(0)
		, min_num_water_triangles(numeric_limits<int>::max())
		, max_num_water_triangles(0)
		, ave_num_water_triangles(0)
		, ave_num_water_triangles_sum(0)
		, render_count(0)
		, info_text("")
		, info_text_def("")
		, light_manager(new OPENGL_LIGHT_MANAGER(driver))
		, gl_ground(new OPENGL_GROUND(driver))
		, gl_simulation_box(new OPENGL_SIMULATION_BOX(driver))
		, bottom_text(new OPENGL_2D_BLOCK_TEXT(FONT_TIMES_NEW_ROMAN, 25, STYLE_REGULAR, PIXMAP_FONT, TEXT_BLOCK_BOTTOM))
		, right_text(new OPENGL_2D_BLOCK_TEXT(FONT_CONSOLAS, 14, STYLE_BOLD, PIXMAP_FONT, TEXT_BLOCK_RIGHT))
	{}

	~OPENGL_WORLD(void)
	{
		DeleteAllObjects();

		if (gl_simulation_box)
		{
			delete gl_simulation_box;
		}

		if (gl_ground)
		{
			delete gl_ground;
		}

		if (light_manager)
		{
			delete light_manager;
		}

		if (driver)
		{
			delete driver;
		}
	}

public: // Initialization Function
	void Initialize(SIMULATION_WORLD* simulation_world_input)
	{
		simulation_world = simulation_world_input;
		if (simulation_world == 0)
		{
			return;
		}

		app_abs_dir = simulation_world->app_abs_dir;
		script_abs_dir = simulation_world->script_abs_dir;
		multithreading = (simulation_world->multithreading);
		grid_scale = simulation_world->grid_scale;

		driver->InitGenericGL();
		light_manager->InitLights();

		menu_object->Init();
		menu_draw_type->Init();
		bottom_text->Init();
		right_text->Init();

		DeleteAllObjects();

		// Create 
		if (simulation_world->use_eulerian_solver)
		{
			opengl_fluid_solver = new OPENGL_EULERIAN_FLUID_SOLVER(driver, &(simulation_world->eulerian_solver), multithreading, grid_scale);
			all_solvers.push_back(opengl_fluid_solver);
		}

		// For periodic boundary condition
		is_periodic = simulation_world->eulerian_solver.water_levelset->is_periodic;
		periodic_num_x = simulation_world->eulerian_solver.water_levelset->periodic_num_x;
		periodic_num_y = simulation_world->eulerian_solver.water_levelset->periodic_num_y;

		for (vector<OPENGL_SOLVER_BASE*>::iterator it = all_solvers.begin(); it != all_solvers.end(); ++it)
		{
			all_objects.insert(all_objects.end(), (*it)->GetObjects().begin(), (*it)->GetObjects().end());
		}

		// Simulation grid min/max
		GRID_STRUCTURE_3D& base_grid = simulation_world->world_discretization.world_grid;
		
		if (is_periodic)
		{
			base_grid.Initialize(periodic_num_x*base_grid.i_res, periodic_num_y*base_grid.j_res, base_grid.k_res, base_grid.i_start, base_grid.j_start, base_grid.k_start, base_grid.x_min, base_grid.y_min, base_grid.z_min, base_grid.x_max + (periodic_num_x - 1)*(base_grid.x_max - base_grid.x_min), base_grid.y_max + (periodic_num_y - 1)*(base_grid.y_max - base_grid.y_min), base_grid.z_max); 
		}

		float half_x = (base_grid.x_max - base_grid.x_min)/2.0f;
		float half_y = (base_grid.y_max - base_grid.y_min)/2.0f;
		float half_z = (base_grid.z_max - base_grid.z_min)/2.0f;
		float cen_x = (base_grid.x_max + base_grid.x_min)/2.0f;
		float cen_y = (base_grid.y_max + base_grid.y_min)/2.0f;
		float cen_z = (base_grid.z_max + base_grid.z_min)/2.0f;

		for (vector<OPENGL_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
		{
			(*it)->SetLength(half_x, half_y, half_z);
			(*it)->SetCenter(cen_x, cen_y, cen_z);
		}

		gl_simulation_box->SetPosition(cen_x, cen_y, cen_z);
		gl_simulation_box->SetLength(abs(half_x), abs(half_y), abs(half_z));
		gl_ground->SetPosition(cen_x, base_grid.y_min, cen_z);
		gl_ground->CalWidth(base_grid);

		LoadOpenGLSettings();
	}

public: // Member Functions
	void Update()
	{
		for (vector<OPENGL_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
		{
			(*it)->Update();
		}

		UpdateDrawInfo();
	}

	void Render()
	{
		Draw();
	}

	void Draw()
	{
		// Gradient Background
		driver->SetDefaultRenderStates2DMode();
		driver->DrawGradientBackground(bg_upper_color.r, bg_upper_color.g, bg_upper_color.b, bg_lower_color.r, bg_lower_color.g, bg_lower_color.b);

		// Light
		driver->SetDefaultRenderStates3DMode();
		light_manager->UpdateLightPosition(light_draw_for_debug);

		// Center axis
		if (draw_center_axis)
		{
			driver->SetDefaultRenderStates3DMode();
			driver->DrawCenterAxis();
		}

		// ETC
		if (draw_ground)
		{
			gl_ground->Draw();
		}

		if (draw_simulation_box)
		{
			gl_simulation_box->Draw();
		}

		// Main Draw
		for (vector<OPENGL_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
		{
			glEnable(GL_BLEND);
			glBlendFunc(GL_ONE, GL_ZERO);
			if (draw_levelset)
			{
				if ((*it)->is_levelset)
				{
					(*it)->Draw();
				}
			}
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			if (draw_boundary)
			{
				if ((*it)->is_boundary)
				{
					(*it)->Draw();
				}
			}
			if (draw_velocity)
			{
				if ((*it)->is_velocity)
				{
					(*it)->Draw();
				}
			}
			if (draw_velocity_x)
			{
				if ((*it)->is_velocity_x)
				{
					(*it)->Draw();
				}
			}
			if (draw_velocity_y)
			{
				if ((*it)->is_velocity_y)
				{
					(*it)->Draw();
				}
			}
			if (draw_velocity_z)
			{
				if ((*it)->is_velocity_z)
				{
					(*it)->Draw();
				}
			}
			if (draw_pressure)
			{
				if ((*it)->is_pressure)
				{
					(*it)->Draw();
				}
			}
		}

		// Corner axis
		if (draw_corner_axis)
		{
			driver->SetDefaultRenderStates2DMode();
			driver->DrawCornerAxis();
		}

		// Menu
		DrawMenu();

		// Text
		if (draw_info)
		{
			DrawInfoText();
		}
	}

	void DrawMenu()
	{
		if (!is_show_menu)
		{
			return;
		}

		OPENGL_COLOR text_color(0.1f, 0.1f, 0.1f);
		
		OPENGL_COLOR select_text_color(0.1f, 0.1f, 0.1f);
		OPENGL_COLOR select_hg_color(1.0f, 0.1f, 0.1f);
		OPENGL_COLOR select_bg_color(0.8f, 0.9f, 0.8f);

		OPENGL_COLOR unsele_text_color(0.3f, 0.3f, 0.3f);
		OPENGL_COLOR unsele_hg_color(0.5f, 0.4f, 0.4f);
		OPENGL_COLOR unsele_bg_color(0.8f, 0.8f, 0.8f);

		// objects menu
		driver->SetDefaultRenderStates2DMode();
		if (all_objects.size() > 0)
		{
			OPENGL_COLOR select_color(1.0f, 1.0f, 1.0f);
			OPENGL_COLOR bg_color(0.8f, 0.9f, 0.8f);
			string str;

			for (vector<OPENGL_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
			{
				str += (*it)->GetName();
				str += "\n";
			}
			str = str.substr(0, str.length() - 1);

			if (is_object_menu)
			{
				menu_object->PrintMenu(str.c_str(), select_text_color, selected_object_index, select_hg_color, true, select_bg_color);
			}
			else
			{
				menu_object->PrintMenu(str.c_str(), unsele_text_color, selected_object_index, unsele_hg_color, true, unsele_bg_color);
			}
		}

		// Draw type menu
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			OPENGL_COLOR select_color(0.1f, 0.8f, 0.9f);
			OPENGL_COLOR bg_color(0.9f, 0.8f, 0.8f);
			string str;
			int selected_draw_type = selected_obj->GetDrawType();
			int draw_type_index = 0;
			int index_count = 0;
			for (map<int, string>::iterator it = selected_obj->GetDrawTypeMap().begin(); it != selected_obj->GetDrawTypeMap().end(); ++it)
			{
				if (selected_draw_type == (*it).first)
				{
					draw_type_index = index_count;
				}

				str += (*it).second;
				str += "\n";

				index_count++;
			}
			str = str.substr(0, str.length() - 1);

			if (!is_object_menu)
			{
				menu_draw_type->PrintSubMenu(menu_object->GetSubMenuPos(selected_object_index), str.c_str(), select_text_color, draw_type_index, select_hg_color, true, select_bg_color);
			}
			else
			{
				menu_draw_type->PrintSubMenu(menu_object->GetSubMenuPos(selected_object_index), str.c_str(), unsele_text_color, draw_type_index, unsele_hg_color, true, unsele_bg_color);
			}
		}
	}

	void DrawInfoText()
	{
		if (info_text.empty())
		{
			if (info_text_def.empty())
			{
				return;
			}
			
			driver->SetDefaultRenderStates2DMode();
			right_text->PrintBlockText(OPENGL_COLOR(0.2, 0.1, 0.1), true, OPENGL_COLOR(0.9, 0.9, 1.0), info_text_def.c_str());
		}
		else
		{
			driver->SetDefaultRenderStates2DMode();
			right_text->PrintBlockText(OPENGL_COLOR(0.2, 0.1, 0.1), true, OPENGL_COLOR(0.9, 0.9, 1.0), info_text.c_str());
		}
	}

	void DrawWithName()
	{
		for (vector<OPENGL_OBJECT_BASE*>::iterator it = all_objects.begin(); it !=  all_objects.end(); ++it)
		{
			(*it)->DrawWithName();
		}
	}

	void ToggleMenu()
	{
		ResetPickObjectID();
		is_show_menu = !is_show_menu;
		is_object_menu = true;
		selected_object_index = 0;
	}

	void Key(unsigned char c)
	{
		for (vector<OPENGL_SOLVER_BASE*>::iterator it = all_solvers.begin(); it != all_solvers.end(); ++it)
		{
			(*it)->Key(c);
		}

		switch (c)
		{
		// Enter Key
		case 13:
			ToggleMenu();
			break;
		case '0':
			LoadOpenGLSettings();
			break;
		case '=':
		case '+': 
			IncrementValue();
			break;
		case '-':
			DecrementValue();
			break;
		case '[':
		case '/': 
			LeftShiftValue();
			break;
		case ']':
		case '*':
			RightShiftValue();
		case ',':
			draw_center_axis = !draw_center_axis;
			break;
		case '.':
			draw_corner_axis = !draw_corner_axis;
			break;
		case ';':
			draw_ground = !draw_ground;
			break;
		case '\'':
			draw_simulation_box = !draw_simulation_box;
			break;
		case 'i':
			draw_info = !draw_info;
			break;
		case 'l':
			draw_levelset = !draw_levelset;
			break;
		case 'b':
			draw_boundary = !draw_boundary;
			break;
		case 'v':
			draw_velocity = !draw_velocity;
			break;
		case 'o':
			draw_velocity_x = !draw_velocity_x;
			break;
		case 'j':
			draw_velocity_y = !draw_velocity_y;
			break;
		case 'k':
			draw_velocity_z = !draw_velocity_z;
			break;
		case 'u':
			draw_pressure = !draw_pressure;
			break;
		}
	}

	void KeyWithAlt(unsigned char c)
	{
		for (vector<OPENGL_SOLVER_BASE*>::iterator it = all_solvers.begin();  it != all_solvers.end(); ++it)
		{
			(*it)->KeyWithAlt(c);
		}

		// Lighting on/off for debugging
		switch (c)
		{
		case '1':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '2':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '3':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '4':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '5':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '6':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '7':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		case '8':
			light_manager->ToggleEnableLight((int) (c - '1')); 
			break;
		}
	}

	void SpecialKey(int key)
	{
		switch (key)
		{
		case GLUT_KEY_UP:
		case GLUT_KEY_DOWN:
		case GLUT_KEY_LEFT:
		case GLUT_KEY_RIGHT:
			// ToggleMenu();
			break;
		}
	}

	void Report()
	{
		++render_count;

		if (opengl_fluid_solver)
		{
			MARCHING_CUBES_ALGORITHM* marching_cube = 0;

			if (opengl_fluid_solver->water_levelset_upsampled)
			{
				marching_cube = opengl_fluid_solver->water_levelset_upsampled->MarchingCubeAlgorithm();
			}
			else
			{
				marching_cube = opengl_fluid_solver->water_levelset->MarchingCubeAlgorithm();
			}

			int num_tri = marching_cube->triangles.num_elements;

			if (min_num_water_triangles > num_tri)
			{
				min_num_water_triangles = num_tri;
			}
			if (max_num_water_triangles < num_tri)
			{
				max_num_water_triangles = num_tri;
			}

			ave_num_water_triangles += num_tri/500;

			cout << "min Number of water triangles = " << min_num_water_triangles << endl;
			cout << "max Number of water triangles = " << max_num_water_triangles << endl;
			cout << "ave Number of water triangles = " << ave_num_water_triangles << endl;
		}
	}

	void PrepareGL()
	{
		driver->InitGenericGL();
		LoadOpenGLSettings(true);

		gl_simulation_box->SetPosition(0, 0, 0);
		gl_simulation_box->SetLength(0.5f, 0.5f, 0.5f);
		gl_ground->SetPosition(0, -0.5f, 0);
	}

	void InitGL()
	{
		// driver
		driver->InitGenericGL();
		driver->SetDefaultRenderStates3DMode();

		// lights
		light_manager->InitLights();
		light_manager->UpdateLightPosition(false, true);
	}

	void SetSelectedObject(OPENGL_OBJECT_BASE* obj)
	{
		selected_object = obj;
	}

	OPENGL_OBJECT_BASE* GetSelectedObject()
	{
		return selected_object;
	}

	OPENGL_DRIVER* GetDriver()
	{
		return driver;
	}

	bool IsDrawCenterAxis()
	{
		return draw_center_axis;
	}

	bool IsDrawCornerAxis()
	{
		return draw_corner_axis;
	}

	bool IsDrawSimulationBox()
	{
		return draw_simulation_box;
	}

	bool IsDrawGround()
	{
		return draw_ground;
	}

	void SetPickObjectID(int pick_id)
	{
		driver->SetSelectedNameId(pick_id);
	}

	void ResetPickObjectID()
	{
		driver->SetSelectedNameId(-1);
	}

	void SetInfoText(string& info)
	{
		info_text = info;
	}

	void DeleteAllObjects()
	{
		all_objects.swap(vector<OPENGL_OBJECT_BASE*>());

		for (vector<OPENGL_SOLVER_BASE*>::iterator it = all_solvers.begin(); it != all_solvers.end(); ++it)
		{
			delete (*it);
		}
		
		all_solvers.swap(vector<OPENGL_SOLVER_BASE*>());
	}

	void IncrementValue()
	{
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			selected_obj->UserAction(OPENGL_OBJECT_BASE::ACTION_INCREMENT);
		}
	}

	void DecrementValue()
	{
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			selected_obj->UserAction(OPENGL_OBJECT_BASE::ACTION_DECREMENT);
		}
	}

	void LeftShiftValue()
	{
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			selected_obj->UserAction(OPENGL_OBJECT_BASE::ACTION_LEFT);
		}
	}

	void RightShiftValue()
	{
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			selected_obj->UserAction(OPENGL_OBJECT_BASE::ACTION_RIGHT);
		}
	}

	void NextObject()
	{
		SelectNextDrawType();
	}

	void PrevObject()
	{
		SelectPrevDrawType();
	}

	void SelectNextObject()
	{
		ResetPickObjectID();
		selected_object_index++;
		if (selected_object_index >= all_objects.size())
		{
			selected_object_index = 0;
		}
	}

	void SelectPrevObject()
	{
		ResetPickObjectID();
		selected_object_index--;
		if (selected_object_index < 0)
		{
			selected_object_index = (int)(all_objects.size() - 1);
		}
	}

	void SelectNextDrawType()
	{
		ResetPickObjectID();
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			selected_obj->NextDrawType();
		}
	}

	void SelectPrevDrawType()
	{
		ResetPickObjectID();
		OPENGL_OBJECT_BASE* selected_obj = all_objects[selected_object_index];
		if (selected_obj)
		{
			selected_obj->PreviousDrawType();
		}
	}

	void UpdateDrawInfo()
	{
		std::ostringstream string_stream;

		if (simulation_world)
		{
			string_stream << "Current Frame: " << simulation_world->num_current_frame;
		}

		// Need to be updated when you need the information
		info_text_def = string_stream.str();
	}

	void LoadOpenGLSettings(bool bFirst = false)
	{
		char buffer[MAX_PATH];
		GetCurrentDirectory(MAX_PATH, buffer);

		if (!script_abs_dir.empty())
		{
			SetCurrentDirectory(script_abs_dir.c_str());
			if (!FileExists("opengl.gl"))
			{
				SetCurrentDirectory(app_abs_dir.c_str());
			}
		}
		else
		{
			SetCurrentDirectory(app_abs_dir.c_str());
		}

		SCRIPT_READER script_reader;
		bool is_succeed = script_reader.Initialize("opengl.gl");
		
		if (is_succeed)
		{
			SCRIPT_BLOCK script_block = script_reader.FindBlock("OPENGL_WORLD");
			light_draw_for_debug = script_block.GetBoolean("light_draw_for_debug");
			bool add_default_point_lights = script_block.GetBoolean("add_default_point_lights");

			// BackGround
			SCRIPT_BLOCK* bg_property = script_block.SearchBlock("BACKGROUND");
			VT bg_upper = bg_property->GetVector3("upper_color");
			bg_upper_color.Set(bg_upper.x, bg_upper.y, bg_upper.z);
			VT bg_lower = bg_property->GetVector3("lower_color");
			bg_lower_color.Set(bg_lower.x, bg_lower.y, bg_lower.z);
			
			// Anti-Alising
			SCRIPT_BLOCK* anti_property = script_block.SearchBlock("ANTI_ALISING");
			bool use_multisample = anti_property->GetBoolean("use_multisample", false);
			bool use_line_smooth = anti_property->GetBoolean("use_line_smooth", false);
			bool use_point_smooth = anti_property->GetBoolean("use_point_smooth", false);

			ANTIALISING_MODE anti_mode = ANTIALISING_OFF;
			if (use_multisample)
			{
				anti_mode = (ANTIALISING_MODE) (anti_mode | ANTIALISING_MULTISAMPLE);
			}
			if (use_line_smooth)
			{
				anti_mode = (ANTIALISING_MODE) (anti_mode | ANTIALISING_LINE_SMOOTH);
			}
			if (use_point_smooth)
			{
				anti_mode = (ANTIALISING_MODE) (anti_mode | ANTIALISING_POINT_SMOOTH);
			}
			driver->SetAntialising(anti_mode);

			// OpenGL Light
			if (bFirst)
			{
				light_manager->LoadOpenGLSettings(&script_block, simulation_world->world_discretization.world_grid, false);
			}
			else
			{
				light_manager->LoadOpenGLSettings(&script_block, simulation_world->world_discretization.world_grid, add_default_point_lights);
			}

			// Solver List
			for (vector<OPENGL_SOLVER_BASE*>::iterator it = all_solvers.begin(); it != all_solvers.end(); ++it)
			{
				(*it)->LoadOpenGLSettings(&script_block);
			}
		}
		SetCurrentDirectory(buffer);
	}

	OPENGL_2D_BLOCK_TEXT* Get2DText()
	{
		return bottom_text;
	}

	OPENGL_2D_BLOCK_TEXT* GetInfoText()
	{
		return right_text;
	}
};

