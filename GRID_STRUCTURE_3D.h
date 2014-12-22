#pragma once

#include "COMMON_DEFINITIONS.h"
#include "SCRIPT_READER.h"
#include "ARRAY.h"
#include "BOX.h"

class GRID_STRUCTURE_3D
{
public: // Essential Data
	union		// Grid Resolution
	{
		struct{int i_res, j_res, k_res;};
		int res[3];
	};

	union		// Start Indices
	{
		struct{int i_start, j_start, k_start;};
		int start_ix[3];
	};

	union		// End indices
	{
		struct{int i_end, j_end, k_end;};
		int end_ix[3];
	};

	union		// Grid Domain (Excluding the ghost cells)
	{
		struct{T x_min, y_min, z_min, x_max, y_max, z_max;};
		struct{T min[3], max[3];};
	};

	union		// Grid Spacing
	{
		struct{T dx, dy, dz;};
		T dxdydz[3];
	};

	union		// Square of Grid Spacing
	{
		struct{T dx2, dy2, dz2;};
		T dx2dy2dz2[3];
	};

	union		// Inverse of grid spacing
	{
		struct{T one_over_dx, one_over_dy, one_over_dz;};
		T one_over_dw[3];
	};

	union		// Inverse of grid spacing
	{
		struct{T one_over_2dx, one_over_2dy, one_over_2dz;};
		T one_over_2dw[3];
	};

	union		// Inverse of square of grid spacing
	{
		struct{T one_over_dx2, one_over_dy2, one_over_dz2;};
		T one_over_dw2[3];
	};

	// Speed up variables
	int ij_res, ijk_res;

public: // Constructors and Destructor
	GRID_STRUCTURE_3D(void)
	{}

	GRID_STRUCTURE_3D(const int& i_res_input, const int& j_res_input, const int& k_res_input, const int& i_start_input, const int& j_start_input, const int& k_start_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input)
	{
		Initialize(i_res_input, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input);
	}

	GRID_STRUCTURE_3D(const VI& res_input, const VI& start_input, const VT& min_input, const VT& max_input)
	{
		Initialize(res_input, start_input, min_input, max_input);
	}

	GRID_STRUCTURE_3D(const GRID_STRUCTURE_3D& grid_input)
	{
		Initialize(grid_input);
	}

	~GRID_STRUCTURE_3D(void)
	{}

public: // Initialization Function
	void Initialize(const int& i_res_input, const int& j_res_input, const int& k_res_input, const int& i_start_input, const int& j_start_input, const int& k_start_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input)
	{
		i_res = i_res_input;
		j_res = j_res_input;
		k_res = k_res_input;

		i_start = i_start_input;
		j_start = j_start_input;
		k_start = k_start_input;

		i_end = i_start + i_res - 1;
		j_end = j_start + j_res - 1;
		k_end = k_start + k_res - 1;

		x_min = x_min_input;
		y_min = y_min_input;
		z_min = z_min_input;

		x_max = x_max_input;
		y_max = y_max_input;
		z_max = z_max_input;

		T diff_x = x_max - x_min, diff_y = y_max - y_min, diff_z = z_max - z_min;
		
		dx = diff_x/(i_res - 1);
		dy = diff_y/(j_res - 1);
		dz = diff_z/(k_res - 1);

		dx2 = dx*dx;
		dy2 = dy*dy;
		dz2 = dz*dz;

		one_over_dx = (T)1/dx;
		one_over_dy = (T)1/dy;
		one_over_dz = (T)1/dz;

		one_over_2dx = (T)0.5/dx;
		one_over_2dy = (T)0.5/dy;
		one_over_2dz = (T)0.5/dz;

		one_over_dx2 = (T)1/dx2;
		one_over_dy2 = (T)1/dy2;
		one_over_dz2 = (T)1/dz2;

		ij_res = i_res*j_res;
		ijk_res = ij_res*k_res;
	}

	void Initialize(const VI& res_input, const VI& start_input, const VT& min_input, const VT& max_input)
	{
		Initialize(res_input.i, res_input.j, res_input.k, start_input.i, start_input.j, start_input.k, min_input.x, min_input.y, min_input.z, max_input.x, max_input.y, max_input.z);
	}
	
	void Initialize(const GRID_STRUCTURE_3D& grid_input)
	{
		Initialize(grid_input.i_res, grid_input.j_res, grid_input.k_res, grid_input.i_start, grid_input.j_start, grid_input.k_start, grid_input.x_min, grid_input.y_min, grid_input.z_min, grid_input.x_max, grid_input.y_max, grid_input.z_max);
	}
	
	void Initialize(const GRID_STRUCTURE_3D& grid_input, const T& resolution_scale)
	{
		Initialize((int)((T)grid_input.i_res*resolution_scale), (int)((T)grid_input.j_res*resolution_scale), (int)((T)grid_input.k_res*resolution_scale), grid_input.i_start, grid_input.j_start, grid_input.k_start, grid_input.x_min, grid_input.y_min, grid_input.z_min, grid_input.x_max, grid_input.y_max, grid_input.z_max);
	}

	void Initialize(const BOX& box_input, const int& max_res)
	{
		const VT edge_length = box_input.max - box_input.min;

		if (edge_length.z >= edge_length.x && edge_length.z >= edge_length.y)
		{
			const int z_res = max_res;
			const T dz = edge_length.z/(T)z_res;
			const int y_res = (int)ceil(edge_length.y/dz);
			const int x_res = (int)ceil(edge_length.x/dz);
			const VT corrected_max(box_input.min.x + (T)x_res*dz, box_input.min.y + (T)y_res*dz, box_input.min.z + (T)z_res*dz);	// To make dx = dy = dz

			Initialize(x_res, y_res, z_res, 0, 0, 0, box_input.min.x, box_input.min.y, box_input.min.z, corrected_max.x, corrected_max.y, corrected_max.z);
		}
		else if (edge_length.y >= edge_length.x)
		{
			const int y_res = max_res;
			const T dy = edge_length.y/(T)y_res;
			const int z_res = (int)ceil(edge_length.z/dy);
			const int x_res = (int)ceil(edge_length.x/dy);
			const VT corrected_max(box_input.min.x + (T)x_res*dy, box_input.min.y + (T)y_res*dy, box_input.min.z + (T)z_res*dy);

			Initialize(x_res, y_res, z_res, 0, 0, 0, box_input.min.x, box_input.min.y, box_input.min.z, corrected_max.x, corrected_max.y, corrected_max.z);
		}
		else
		{
			const int x_res = max_res;
			const T dx = edge_length.x/(T)x_res;
			const int z_res = (int)ceil(edge_length.z/dx);
			const int y_res = (int)ceil(edge_length.y/dx);
			const VT corrected_max(box_input.min.x + (T)x_res*dx, box_input.min.y + (T)y_res*dx, box_input.min.z + (T)z_res*dx);

			Initialize(x_res, y_res, z_res, 0, 0, 0, box_input.min.x, box_input.min.y, box_input.min.z, corrected_max.x, corrected_max.y, corrected_max.z);
		}
	}

	void Initialize(const BOX& box_input, const T& dx_input)
	{
		const VT edge_length = box_input.max - box_input.min;
		const int x_res = (int)ceil(edge_length.x/dx_input);
		const int y_res = (int)ceil(edge_length.y/dx_input);
		const int z_res = (int)ceil(edge_length.z/dx_input);
		const VT corrected_max(box_input.min.x + (T)x_res*dx_input, box_input.min.y + (T)y_res*dx_input, box_input.min.z + (T)z_res*dx_input);

		Initialize(x_res, y_res, z_res, 0, 0, 0, box_input.min.x, box_input.min.y, box_input.min.z, corrected_max.x, corrected_max.y, corrected_max.z);
	}
	
	void InitializeFromBlock(const SCRIPT_BLOCK& block)
	{
		const T resolution_scale = block.GetFloat("resolution_scale", (T)1);
		const VI start = block.GetInt3("start_indices");
		const VI res = block.GetInt3("base_grid_resolution");
		const VI res_scaled((int)((T)res.x*resolution_scale), (int)((T)res.y*resolution_scale), (int)((T)res.z*resolution_scale));
		const VT min = block.GetVector3("base_grid_min");
		const VT max = block.GetVector3("base_grid_max");

		// Display
		cout << "resolution scale: " << resolution_scale << endl;
		cout << "start indices: (" << start.i << " ," << start.j << " ," << start.k << ") " << endl;
		cout << "resolution: (" << res.i << " ," << res.j << " ," << res.k << ") " << endl;
		cout << "scaled resolution: (" << res_scaled.i << " ," << res_scaled.j << " ," << res_scaled.k << ") " << endl;
		cout << "base grid min: (" << min.i << " ," << min.j << " ," << min.k << ") " << endl;
		cout << "base grid max: (" << max.i << " ," << max.j << " ," << max.k << ") " << endl;

		Initialize(res_scaled.i, res_scaled.j, res_scaled.k, start.i, start.j, start.k, min.i, min.j, min.k, max.i, max.j, max.k);
	}

public: // Operator Overloading
	void operator = (const GRID_STRUCTURE_3D& grid_input)
	{
		Initialize(grid_input);
	}

public: // Indexing Functions
	inline const VI ClampedIndex(const int& i, const int& j, const int& k) const
	{
		return VI(CLAMP(i, i_start, i_end), CLAMP(j, j_start, j_end), CLAMP(k, k_start, k_end));
	}

	inline const VI ClampedIndex(const VI& ix) const
	{
		return VI(CLAMP(ix.i, i_start, i_end), CLAMP(ix.j, j_start, j_end), CLAMP(ix.k, k_start, k_end));
	}

	inline const int Index1D(const int& i, const int& j, const int& k) const
	{
		return (i - i_start) + (j - j_start)*i_res + (k - k_start)*ij_res;
	}

	inline const VI Index3D(const int& index_1d) const
	{
		const int n = index_1d/ij_res;
		const int n_ij_res = n*ij_res;
		const int m = (index_1d - n_ij_res)/i_res;
		const int l = index_1d - i_res*m - n_ij_res;

		return VI(l,m,n);
	}

	inline const int ClampedIndexI(const int& i) const
	{
		return CLAMP(i, i_start, i_end);
	}

	inline const int ClampedIndexJ(const int& j) const
	{
		return CLAMP(j, j_start, j_end);
	}

	inline const int ClampedIndexK(const int& k) const
	{
		return CLAMP(k, k_start, k_end);
	}

public: // Member Functions
	inline const VT CellCenter(const int& i, const int& j, const int& k) const
	{
		return VT(x_min + ((T)0.5 + (T)(i - i_start))*dx, y_min + ((T)0.5 + (T)(j - j_start))*dy, z_min + ((T)0.5 + (T)(k - k_start))*dz);
	}

	inline const VT CellCenter(const VI& ix) const
	{
		return VT(x_min + ((T)0.5 + (T)(ix.i - i_start))*dx, y_min + ((T)0.5 + (T)(ix.j - j_start))*dy, z_min + ((T)0.5 + (T)(ix.k - k_start))*dz);
	}

	inline void Center(const int& i, const int& j, const int& k, VT& position) const
	{
		position.x = x_min + ((T)0.5 + (T)(i - i_start))*dx;
		position.y = y_min + ((T)0.5 + (T)(j - j_start))*dy;
		position.z = z_min + ((T)0.5 + (T)(k - k_start))*dz;
	}

	inline void Center(const VI& ix, VT& position) const 
	{
		position.x = x_min + ((T)0.5 + (T)(ix.i - i_start))*dx;
		position.y = y_min + ((T)0.5 + (T)(ix.j - j_start))*dy;
		position.z = z_min + ((T)0.5 + (T)(ix.k - k_start))*dz;
	}

	// Note that this one returns index containing the cell
	inline const VI Cell(const VT& position) const		 
	{
		return VI(i_start + (int)(floor((position.x - x_min)*one_over_dx)), j_start + (int)(floor((position.y - y_min)*one_over_dy)), k_start + (int)(floor((position.z - z_min)*one_over_dz)));
	}

	inline const void Cell(const VT& position, int& i, int& j, int& k) const
	{
		i = i_start + (int)(floor((position.x - x_min)*one_over_dx));
		j = j_start + (int)(floor((position.y - y_min)*one_over_dy));
		k = k_start + (int)(floor((position.z - z_min)*one_over_dz));

		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);
		assert(k >= k_start && k <= k_end);
	}

	inline const VI ClampedCell(const VT& position) const
	{
		VI index(i_start + (int)(floor((position.x - x_min)*one_over_dx)), j_start + (int)(floor((position.y - y_min)*one_over_dy)), k_start + (int)(floor((position.z - z_min)*one_over_dz)));

		index.i = CLAMP(index.i, i_start, i_end);
		index.j = CLAMP(index.j, j_start, j_end);
		index.k = CLAMP(index.k, k_start, k_end);
	
		return index;
	}

	inline const VI ClampedCell(const VT& position, VI& index) const
	{
		index.Assign(i_start + (int)(floor((position.x - x_min)*one_over_dx)), j_start + (int)(floor((position.y - y_min)*one_over_dy)), k_start + (int)(floor((position.z - z_min)*one_over_dz)));

		index.i = CLAMP(index.i, i_start, i_end);
		index.j = CLAMP(index.j, j_start, j_end);
		index.k = CLAMP(index.k, k_start, k_end);

	  	return index;
	}

	inline const void ClampedCell(const VT& position, int& i, int& j, int& k) const
	{
		i = i_start + (int)(floor((position.x - x_min)*one_over_dx));
		j = j_start + (int)(floor((position.y - y_min)*one_over_dy));
		k = k_start + (int)(floor((position.z - z_min)*one_over_dz));

		i = CLAMP(i, i_start, i_end);
		j = CLAMP(j, j_start, j_end);
		k = CLAMP(k, k_start, k_end);
	}

	inline const VI LeftBottomCell(const VT& position) const
	{
		return VI(i_start + (int)((position.x - x_min)*one_over_dx - (T)0.5), j_start + (int)((position.y - y_min)*one_over_dy - (T)0.5), k_start + (int)((position.z - z_min)*one_over_dz));
	}

	inline void LeftBottomCell(const VT& position, int& i, int& j, int& k) const
	{
		i = i_start + (int)floor((position.x - x_min)*one_over_dx - (T)0.5);
		j = j_start + (int)floor((position.y - y_min)*one_over_dy - (T)0.5);
		k = k_start + (int)floor((position.z - z_min)*one_over_dz - (T)0.5);
	}

	inline void LeftBottomCell(const VT& position, VI& ix) const
	{
		return LeftBottomCell(position, ix.i, ix.j, ix.k);
	}

	inline bool Inside(const VT& position) const
	{
		if(position.x <= x_min) return false;
		else if(position.x >= x_max) return false;
		else if(position.y <= y_min) return false;
		else if(position.y >= y_max) return false;
		else if(position.z <= z_min) return false;
		else if(position.z >= z_max) return false;
		return true;
	}

	inline bool Inside(const VT& position, const T& width) const
	{
		if(position.x <= x_min + width) return false;
		else if(position.x >= x_max - width) return false;
		else if(position.y <= y_min + width) return false;
		else if(position.y >= y_max - width) return false;
		else if(position.z <= z_min + width) return false;
		else if(position.z >= z_max - width) return false;
		return true;
	}

	inline bool Inside(const int& i, const int& j, const int& k) const
	{
		if(i < i_start) return false;
		else if(i > i_end) return false;
		else if(j < j_start) return false;
		else if(j > j_end) return false;
		else if(k < k_start) return false;
		else if(k > k_end) return false;
		return true;
	}

	inline bool Inside(const VI& ix) const
	{
		if(ix.i < i_start) return false;
		else if(ix.i > i_end) return false;
		else if(ix.j < j_start) return false;
		else if(ix.j > j_end) return false;
		else if(ix.k < k_start) return false;
		else if(ix.k > k_end) return false;
		return true;
	}

	inline bool Inside(const VI& ix, const int& inner_width) const
	{
		if(ix.i < i_start + inner_width) return false;
		else if(ix.i > i_end - inner_width) return false;
		else if(ix.j < j_start + inner_width) return false;
		else if(ix.j > j_end - inner_width) return false;
		else if(ix.k < k_start + inner_width) return false;
		else if(ix.k > k_end - inner_width) return false;
		return true;
	}

	inline bool InsideI(const int& i) const
	{
		if(i < i_start) return false;
		else if(i > i_end) return false;
		return true;
	}

	inline bool InsideJ(const int& j) const
	{
		if(j < j_start) return false;
		else if(j > j_end) return false;
		return true;
	}

	inline bool InsideK(const int& k) const
	{
		if(k < k_start) return false;
		else if(k > k_end) return false;
		return true;
	}

	GRID_STRUCTURE_3D Enlarged(const int& width) const
	{
		return GRID_STRUCTURE_3D(i_res + 2*width, j_res + 2*width, k_res + 2*width, i_start - width, j_start - width, k_start - width, x_min - (T)width*dx, y_min - (T)width*dy, z_min - (T)width*dz, x_max + (T)width*dx, y_max + (T)width*dy, z_max + (T)width*dz);
	}

	void Enlarge(const int& width) 
	{
		Initialize(i_res + 2*width, j_res + 2*width, k_res + 2*width, i_start - width, j_start - width, k_start - width, x_min - (T)width*dx, y_min - (T)width*dy, z_min - (T)width*dz, x_max + (T)width*dx, y_max + (T)width*dy, z_max + (T)width*dz);
	}
	
	void Translate(const VT& deviation)
	{
		x_min += deviation.x;
		y_min += deviation.y;
		z_min += deviation.z;

		x_max += deviation.x;
		y_max += deviation.y;
		z_max += deviation.z;
	}
	
public: // Function for multithreading - defined on the cpp file
	void SplitInXDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids);
	void SplitInYDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids);
	void SplitInZDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids);
	void SplitInMaxDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids);

	void SplitInZDierctionPartially(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids);

	int RecommendMaxMultigridLevel(const int& lowest_level_res);
};

// Ostream object overloading
//inline std::ostream& 
//operator<<(std::ostream& output, const GRID_STRUCTURE_3D& grid)
//{
//	return output << "GRID_STRUCTURE_3D ["
//		          << "Resolution = " << grid.i_res << " " << grid.j_res	<< " " << grid.k_res		
//				  << " Index range =(" << grid.i_start << " " << grid.j_start << " " << grid.k_start << ") to (" << grid.i_end << " " << grid.j_end << " " << grid.k_end << ") "
//				  << " Range = (" << grid.x_min << " " << grid.y_min << " " << grid.z_min << ") to (" << grid.x_max << " " << grid.y_max << " " << grid.z_max << ") "
//				  << " DX = " << grid.dx << " " << grid.dy << " " << grid.dz << "]";
//}
