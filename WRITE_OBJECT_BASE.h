#pragma once

#include "TRIANGULAR_SURFACE.h"
#include <string>
#include <Windows.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

enum FILE_FORMAT
{
	FORMAT_NONE = 0,

	PARTICLE_FORMAT_PTC,				// arnold
	PARTICLE_FORMAT_BIN,				// realflow bin
	PARTICLE_FORMAT_OBJ,
	PARTICLE_FORMAT_PBRT,				// pbrt
	PARTICLE_FORMAT_MB_PBRT,			// pbrt (motion blur)
	PARTICLE_FORMAL_LXO,				// luxrender

	PARTICLE_DIFF_FORMAT_PTC,			// arnold
	PARTICLE_DIFF_FORMAT_BIN,			// realflow bin

	RIGIDOBJECT_FORMAT_OBJ,
	RIGIDOBJECT_FORMAT_PBRT,
	RIGIDOBJECT_FORMAT_LXO,				// luxrender

	LEVELSET_FORMAT_MES,				// arnold
	LEVELSET_FORMAT_OBJ,			
	LEVELSET_FORMAT_PBRT,
	LEVELSET_FORMAT_LXO,				// luxrender

	SCALARFIELD_DENSITY_FORMAT_PBRT,	// density PBRT
	// Add the other one when you need
};

struct WRITE_OPTIONS
{
public: // Essential Data
	bool			write_file_on;
	string			write_root_folder;
	string			write_folder_name;
	string			write_file_name_prefix;
	FILE_FORMAT		write_file_format;

	VI				write_grid;
	float			write_velocity_scale_factor;

	bool			write_proxy_file;
	int				write_proxy_reduction_ratio;

public: // Constructor
	WRITE_OPTIONS(void)
		: write_file_on(false)
		, write_folder_name("data")
		, write_file_name_prefix("data")
		, write_file_format(FORMAT_NONE)
		, write_grid(100, 100, 100)
		, write_velocity_scale_factor(0.015f)
		, write_proxy_file(false)
		, write_proxy_reduction_ratio(1)
	{}
};

class WRITE_OBJECT_BASE
{
public: // Essential Data
	WRITE_OPTIONS	write_options;
	string			name;

public: // Constructor and Destructor
	WRITE_OBJECT_BASE(const char* script_name)
		: name(script_name)
	{}

	virtual ~WRITE_OBJECT_BASE(void)
	{}

public: // Member Functions
	WRITE_OPTIONS& GetOptions()
	{
		return write_options;
	}

	void SetOptions(const WRITE_OPTIONS& opt)
	{
		write_options = opt;
	}

	void GenFileFullPath(string& file_full_path, int current_frame, const char* extension, const char* postfix = 0)
	{
		char file_path[MAX_PATH] = {0, };
		char folder_path[MAX_PATH] = {0, };

		if (postfix)
		{
			_snprintf(folder_path, MAX_PATH, "%s/%s_%s", GetOptions().write_root_folder.c_str(), GetOptions().write_folder_name.c_str(), postfix);
			_snprintf(file_path, MAX_PATH, "%s/%s_%s_%04d.%s", folder_path, GetOptions().write_file_name_prefix.c_str(), postfix, current_frame, extension);
		}
		else
		{
			_snprintf(folder_path, MAX_PATH, "%s/%s", GetOptions().write_root_folder.c_str(), GetOptions().write_folder_name.c_str());
			_snprintf(file_path, MAX_PATH, "%s/%s_%04d.%s", folder_path, GetOptions().write_file_name_prefix.c_str(),current_frame, extension);
		}
		CreateDir(folder_path);

		file_full_path = file_path;
	}

	void GenFilePath(string& file_path, int current_frame, const char* extension)
	{
		char file_path_t[MAX_PATH] = {0, };
		_snprintf(file_path_t, MAX_PATH, "%s_%04d.%s", GetOptions().write_file_name_prefix.c_str(), current_frame, extension);

		file_path = file_path_t;
	}

	void GenFileAbsPath(string& file_path, int current_frame, const char* extension)
	{
		char file_path_t[MAX_PATH] = {0, };
		_snprintf(file_path_t, MAX_PATH, "%s/%s_%04d.%s", GetOptions().write_folder_name.c_str(), GetOptions().write_file_name_prefix.c_str(), current_frame, extension);

		file_path = file_path_t;
	}

	void CreateDir(const char* dir_path)
	{
		boost::filesystem::path p(dir_path);
		if ((!boost::filesystem::exists(p)) || (boost::filesystem::exists(p) && boost::filesystem::is_directory(p)))
		{
			boost::filesystem::create_directories(p);
		}
	}

	void ParseWriteOptions(SCRIPT_BLOCK* root_block, const char* write_root_folder)
	{
		if (!root_block)
		{
			return;
		}

		SCRIPT_BLOCK* block = root_block->SearchBlock(name.c_str());
		if (!block)
		{
			return;
		}

		write_options.write_file_on = block->GetBoolean("wirte_file_on", false);
		write_options.write_root_folder = write_root_folder;
		const char* str = block->GetString("write_folder_name", 0);
		if (str)
		{
			write_options.write_folder_name = str;
		}
		else
		{
			write_options.write_folder_name = "data";	// Default
		}

		str = block->GetString("write_file_format", 0);
		if (str)
		{
			if(boost::iequals(str, "particle_ptc")) write_options.write_file_format = PARTICLE_FORMAT_PTC;
			else if(boost::iequals(str, "particle_bin")) write_options.write_file_format = PARTICLE_FORMAT_BIN;
			else if(boost::iequals(str, "particle_obj")) write_options.write_file_format = PARTICLE_FORMAT_OBJ;
			else if(boost::iequals(str, "particle_pbrt")) write_options.write_file_format = PARTICLE_FORMAT_PBRT;
			else if(boost::iequals(str, "particle_mb_pbrt")) write_options.write_file_format = PARTICLE_FORMAT_MB_PBRT;
			
			else if(boost::iequals(str, "particle_diff_bin")) write_options.write_file_format = PARTICLE_DIFF_FORMAT_BIN;
			else if(boost::iequals(str, "particle_diff_obj")) write_options.write_file_format = PARTICLE_DIFF_FORMAT_BIN;

			else if(boost::iequals(str, "rigidobject_obj")) write_options.write_file_format = RIGIDOBJECT_FORMAT_OBJ;
			else if(boost::iequals(str, "rigidobject_pbrt")) write_options.write_file_format = RIGIDOBJECT_FORMAT_PBRT;
			else if(boost::iequals(str, "rigidobject_lxo")) write_options.write_file_format = RIGIDOBJECT_FORMAT_LXO;

			else if(boost::iequals(str, "levelset_mes")) write_options.write_file_format = LEVELSET_FORMAT_MES;
			else if(boost::iequals(str, "levelset_obj")) write_options.write_file_format = LEVELSET_FORMAT_OBJ;
			else if(boost::iequals(str, "levelset_pbrt")) write_options.write_file_format = LEVELSET_FORMAT_PBRT;
			else if(boost::iequals(str, "levelset_lxo")) write_options.write_file_format = LEVELSET_FORMAT_LXO;

			else if(boost::iequals(str, "scalarfield_density_pbrt")) write_options.write_file_format = SCALARFIELD_DENSITY_FORMAT_PBRT;
		}
		else
		{
			write_options.write_file_format = FORMAT_NONE;
		}

		str = block->GetString("write_file_name_prefix", 0);
		if (str)
		{
			write_options.write_file_name_prefix = str;
		}
		else
		{
			write_options.write_file_name_prefix = "data";		// Default
		}

		write_options.write_grid = block->GetInt3("write_grid", VI(100, 100, 100));

		write_options.write_velocity_scale_factor = block->GetFloat("write_velocity_scale_factor", 0.015f);

		write_options.write_proxy_file = block->GetBoolean("write_proxy_file", false);
		write_options.write_proxy_reduction_ratio = block->GetInteger("write_proxy_reduction_ration", 1);

		if (write_options.write_proxy_reduction_ratio < 1)
		{
			write_options.write_proxy_reduction_ratio = 1;
		}
	}

public: // Virtual Functions
	virtual void Write(int current_frame) = 0;
};





