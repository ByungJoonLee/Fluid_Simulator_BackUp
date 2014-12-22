#pragma once

#include "GRID_STRUCTURE_3D.h"

class MAC_GRID_STRUCTURE_3D
{
public: // Essential Data
	GRID_STRUCTURE_3D		x_grid;
	GRID_STRUCTURE_3D		y_grid;
	GRID_STRUCTURE_3D		z_grid;

public: // Speed up variable
	union              // Base Grid Resolution
	{
		struct{int i_res, j_res, k_res;};
		int res[3];
	};
	
	union             // Base Grid Domain
	{
		struct{T x_min, y_min, z_min, x_max, y_max, z_max;};
		struct{T min[3], max[3];};
	};

	union             // Base Grid Start Indices
	{
		struct{int i_start, j_start, k_start;};
		int start_ix[3];
	};
	union              // Grid Resolution x
	{
		struct{int i_res_x, j_res_x, k_res_x;};
		int res_x[3];
	};
	
	union              // Grid Resolution y
	{
		struct{int i_res_y, j_res_y, k_res_y;};
		int res_y[3];
	};

	union              // Grid Resolution z
	{
		struct{int i_res_z, j_res_z, k_res_z;};
		int res_z[3];
	};
	

	union             // Grid Domain x
	{
		struct{T x_min_x, y_min_x, z_min_x, x_max_x, y_max_x, z_max_x;};
		struct{T min_x[3], max_x[3];};
	};

	union             // Grid Domain y
	{
		struct{T x_min_y, y_min_y, z_min_y, x_max_y, y_max_y, z_max_y;};
		struct{T min_y[3], max_y[3];};
	};

	union             // Grid Domain z
	{
		struct{T x_min_z, y_min_z, z_min_z, x_max_z, y_max_z, z_max_z;};
		struct{T min_z[3], max_z[3];};
	};

	union            // Grid Spacing -- Assume that given problem is solved on the uniform grid space
	{
		struct{T dx, dy, dz;};
		T dw[3];
	};

public: // Constructors and Destructor
	MAC_GRID_STRUCTURE_3D(void)
	{}

	MAC_GRID_STRUCTURE_3D(const int& i_res_input, const int& j_res_input, const int& k_res_input, const int& i_start_input, const int& j_start_input, const int& k_start_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input)
	{
		Initialize(i_res_input, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input);
	}


	MAC_GRID_STRUCTURE_3D(const GRID_STRUCTURE_3D& grid_3d_input)
	{
		Initialize(grid_3d_input);
	}

	~MAC_GRID_STRUCTURE_3D(void)
	{}

public: // Initialization Functions
	void Initialize(const int& i_res_input, const int& j_res_input, const int& k_res_input, const int& i_start_input, const int& j_start_input, const int& k_start_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input)
	{
		i_res = i_res_input;
		j_res = j_res_input;
		k_res = k_res_input;

		x_min = x_min_input;
		y_min = y_min_input;
		z_min = z_min_input;

		x_max = x_max_input;
		y_max = y_max_input;
		z_max = z_max_input;
		
		dx = (x_max_input - x_min_input)/i_res_input;
		dy = (y_max_input - y_min_input)/j_res_input;
		dz = (z_max_input - z_min_input)/k_res_input;

		x_grid.Initialize(i_res_input + 1, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input - (T)0.5*dx, y_min_input, z_min_input, x_max_input + (T)0.5*dx, y_max_input, z_max_input);
		y_grid.Initialize(i_res_input, j_res_input + 1, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input - (T)0.5*dy, z_min_input, x_max_input, y_max_input + (T)0.5*dy, z_max_input);
		z_grid.Initialize(i_res_input, j_res_input, k_res_input + 1, i_start_input, j_start_input, k_start_input, x_min_input - (T)0.5*dx, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input + (T)0.5*dz);
		
		i_res_x = x_grid.i_res;
		j_res_x = x_grid.j_res;
		k_res_x = x_grid.k_res;

		i_res_y = y_grid.i_res;
		j_res_y = y_grid.j_res;
		k_res_y = y_grid.k_res;

		i_res_z = z_grid.i_res;
		j_res_z = z_grid.j_res;
		k_res_z = z_grid.k_res;

		x_min_x = x_grid.x_min;
		y_min_x = x_grid.y_min;
		z_min_x = x_grid.z_min;

		x_min_y = y_grid.x_min;
		y_min_y = y_grid.y_min;
		z_min_y = y_grid.z_min;

		x_min_z = z_grid.x_min;
		y_min_z = z_grid.y_min;
		z_min_z = z_grid.z_min;

		x_max_x = x_grid.x_max;
		y_max_x = x_grid.y_max;
		z_max_x = x_grid.z_max;

		x_max_y = y_grid.x_max;
		y_max_y = y_grid.y_max;
		z_max_y = y_grid.z_max;

		x_max_z = z_grid.x_max;
		y_max_z = z_grid.y_max;
		z_max_z = z_grid.z_max;
	}

	void Initialize(const GRID_STRUCTURE_3D& grid_3d_input)
	{
		Initialize(grid_3d_input.i_res, grid_3d_input.j_res, grid_3d_input.k_res, grid_3d_input.i_start, grid_3d_input.j_start, grid_3d_input.k_start, grid_3d_input.x_min, grid_3d_input.y_min, grid_3d_input.z_min, grid_3d_input.x_max, grid_3d_input.y_max, grid_3d_input.z_max);
	}

public: // Member Functions
	MAC_GRID_STRUCTURE_3D Enlarged(const int& width) const
	{
		MAC_GRID_STRUCTURE_3D(i_res + 2*width, j_res + 2*width, k_res + 2*width, i_start - width, j_start - width, k_start - width, x_min - (T)width*dx, y_min - (T)width*dy, z_min - (T)width*dz, x_max + (T)width*dx, y_max + (T)width*dy, z_max + (T)width*dz);
	}
}; 

// As you proceeding the given code, modify this again




	





