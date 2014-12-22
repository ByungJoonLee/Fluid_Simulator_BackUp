#pragma once

#include "WRITE_OBJECT_BASE.h"
#include "LEVELSET_3D.h"
#include "MARCHING_CUBES_ALGORITHM.h"

class WRITE_LEVELSET : public WRITE_OBJECT_BASE
{
public: // Essnetial Data
	LEVELSET_3D* levelset_object;

public: // Constructor and Destructor
	WRITE_LEVELSET(const char* script_name, LEVELSET_3D* levelset_object_input)
		: WRITE_OBJECT_BASE(script_name)
		, levelset_object(levelset_object_input)
	{}

	~WRITE_LEVELSET(void)
	{}

public: // Member Functions
	void WriteOBJ(int current_frame)
	{
		string file_path("");
		GenFileFullPath(file_path, current_frame, "obj");

		GRID_STRUCTURE_3D& grid(levelset_object->grid);
		MARCHING_CUBES_ALGORITHM levelset_polygonizer;
		levelset_polygonizer.Initialize(levelset_object->multithreading, levelset_object->grid);
		levelset_polygonizer.Polygonize(*levelset_object);

		// Write OBJ
		ofstream file(file_path);
		if (!file.is_open())
		{
			cout << "[ERROR] Fail to open file: " << file_path << endl;
			return;
		}

		ARRAY<TRIANGULAR_SURFACE*>& arr = levelset_polygonizer.surfaces;
		ARRAY<int> vertex_count(arr.num_elements);

		// Write vertex
		for (int i = 0; i < arr.num_elements; ++i)
		{
			TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
			int vertices_size = (int)triangular_surface->vertices.size();
			for (int index = 0; index < vertices_size; ++index)
			{
				file << "v " << triangular_surface->vertices[index]->x[0];
				file << " " << triangular_surface->vertices[index]->x[1];
				file << " " << triangular_surface->vertices[index]->x[2];
			}
		}

		// Write normal
		for (int i = 0; i < arr.num_elements; ++i)
		{
			TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
			int vertices_size = (int)triangular_surface->vertices.size();
			for (int index = 0; index < vertices_size; ++index)
			{
				file << "vn " << triangular_surface->vertices[index]->n[0];
				file << " " << triangular_surface->vertices[index]->n[1];
				file << " " << triangular_surface->vertices[index]->n[2];
			}
		}

		// Write index
		for (int i = 0; i < arr.num_elements; ++i)
		{
			TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
			for (list<TRIANGLE*>::iterator itr_triangle = triangular_surface->triangles.begin(); itr_triangle != triangular_surface->triangles.end(); ++itr_triangle)
			{
				int index[3] = {(*itr_triangle)->vertices[0]->index, (*itr_triangle)->vertices[1]->index, (*itr_triangle)->vertices[2]->index};
				file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << endl;
			}
		}

		file.close();
	}

	// Need to be updated when you need more
public: // Virtual Functions
	virtual void Write(int current_frame)
	{
		if (!GetOptions().write_file_on)
		{
			return;
		}

		if (GetOptions().write_file_format == LEVELSET_FORMAT_OBJ)
		{
			WriteOBJ(current_frame);
		}
		// Need to be updated when you need more
	}
};
