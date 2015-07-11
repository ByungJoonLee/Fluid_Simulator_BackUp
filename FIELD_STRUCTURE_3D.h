#pragma once

#include "MAC_GRID_STRUCTURE_3D.h"
#include "GRID_STRUCTURE_3D.h"
#include "ARRAY_3D.h"
#include "MULTITHREADING.h"
#include "LEVELSET_OBJECT.h"
#include "GS_COEFFICIENTS.h"

template<class TT>
class FIELD_STRUCTURE_3D
{
public: // Essential Data
	GRID_STRUCTURE_3D				grid;
	GRID_STRUCTURE_3D				grid_ghost;							// Grid including ghost cells. Same as grid when ghost_width = 0
	ARRAY<GRID_STRUCTURE_3D>		partial_grids;						// Partial grid for multithreading
	ARRAY<GRID_STRUCTURE_3D>		partial_grids_ghost;				// Partial grid for grid_ghost
	MAC_GRID_STRUCTURE_3D			mac_grid;
	MAC_GRID_STRUCTURE_3D			mac_grid_ghost;
	ARRAY<MAC_GRID_STRUCTURE_3D>	partial_grids_mac;
	ARRAY<MAC_GRID_STRUCTURE_3D>	partial_grids_ghost_mac;
	ARRAY_3D<TT>					array_for_this;
	ARRAY_3D<bool>					fixed;

public: // Properties
	int								ghost_width;
	
public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Speed Up Constants and Variables
	int								i_start, j_start, k_start, i_start_g, j_start_g, k_start_g;
	int								i_end, j_end, k_end, i_end_g, j_end_g, k_end_g;
	int								i_res_g, ij_res_g;
	T								x_min, y_min, z_min;
	T								dx, dy, dz;
	T								one_over_dx, one_over_dy, one_over_dz;
	T								one_over_2dx, one_over_2dy, one_over_2dz;
	T								one_over_dx2, one_over_dy2, one_over_dz2;
	TT*								values;

public: // Choosing the type
	bool							is_scalar;
	bool							is_vector;

public: // Choosing which component is
	bool							is_x_component;
	bool							is_y_component;
	bool							is_z_component;

public: // Constructors and Destructor
	FIELD_STRUCTURE_3D(void)
		: multithreading(0), is_scalar(true), is_vector(false), is_x_component(false), is_y_component(false), is_z_component(false)
	{}

	FIELD_STRUCTURE_3D(MULTITHREADING* multithreading_input, const VI& res_input, const VI& start_input, const VT& min_input, const VT& max_input, const int& ghost_width_input = 0)
		: partial_grids(0), is_scalar(true), is_vector(false), is_x_component(false), is_y_component(false), is_z_component(false)
	{
		Initialize(multithreading_input, res_input, start_input, min_input, max_input, ghost_width);
	}

	~FIELD_STRUCTURE_3D(void)
	{}

public: // Initialization Functions
	void Initialize(MULTITHREADING* multithreading_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const int& i_start_input, const int& j_start_input, const int& k_start_input, 
					const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input, const int& ghost_width_input = 0, const bool& is_scalar_input = true, const bool& is_vector_input = false)
	{
		// For general grid
		grid.Initialize(i_res_input, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input);
		
		ghost_width = ghost_width_input;

		grid_ghost.Initialize(grid.Enlarged(ghost_width));
		
		mac_grid.Initialize(grid);

		i_start = grid.i_start;
		j_start = grid.j_start;
		k_start = grid.k_start;

		i_end = grid.i_end;
		j_end = grid.j_end;
		k_end = grid.k_end;

		i_start_g = grid_ghost.i_start;
		j_start_g = grid_ghost.j_start;
		k_start_g = grid_ghost.k_start;

		i_end_g = grid_ghost.i_end;
		j_end_g = grid_ghost.j_end;
		k_end_g = grid_ghost.k_end;

		x_min = grid.x_min;
		y_min = grid.y_min;
		z_min = grid.z_min;

		dx = grid.dx;
		dy = grid.dy;
		dz = grid.dz;

		one_over_dx = grid.one_over_dx;
		one_over_dy = grid.one_over_dy;
		one_over_dz = grid.one_over_dz;

		one_over_2dx = grid.one_over_2dx;
		one_over_2dy = grid.one_over_2dy;
		one_over_2dz = grid.one_over_2dz;

		one_over_dx2 = grid.one_over_dx2;
		one_over_dy2 = grid.one_over_dy2;
		one_over_dz2 = grid.one_over_dz2;

		// Let array_for_this have the information for grid_ghost
		array_for_this.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.k_start - ghost_width, 
								  grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, grid.k_res + 2*ghost_width, true);
		
		fixed.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.k_start - ghost_width, 
								  grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, grid.k_res + 2*ghost_width, true);

		fixed.AssignAllValues((bool)false);

		i_start_g = array_for_this.i_start;
		j_start_g = array_for_this.j_start;
		k_start_g = array_for_this.k_start;

		i_end_g = array_for_this.i_end;
		j_end_g = array_for_this.j_end;
		k_end_g = array_for_this.k_end;
		
		values = array_for_this.values;

		i_res_g = array_for_this.i_res;
		ij_res_g = array_for_this.ij_res;

		is_scalar = is_scalar_input;
		is_vector = is_vector_input;

		// Multithreading
		assert(multithreading_input != 0);
		multithreading = multithreading_input;
		grid.SplitInZDirection(multithreading->num_threads, partial_grids);
		grid_ghost.SplitInZDirection(multithreading->num_threads, partial_grids_ghost);
	}
	
	void Initialize(MULTITHREADING* multithreading_input, const VI& res_input, const VI& start_input, const VT& min_input, const VT& max_input, const int& ghost_width_input = 0)
	{
		Initialize(multithreading_input, res_input.i, res_input.j, res_input.k, start_input.i, start_input.j, start_input.k, min_input.i, min_input.j, min_input.k, max_input.i, max_input.j, max_input.k, ghost_width_input);
	}

	void Initialize(MULTITHREADING* multithreading_input, const GRID_STRUCTURE_3D& grid_input, const int& ghost_width = 0, const bool& is_scalar_input = true, const bool& is_vector_input = false)
	{
		Initialize(multithreading_input, grid_input.i_res, grid_input.j_res, grid_input.k_res, grid_input.i_start, grid_input.j_start, grid_input.k_start, grid_input.x_min, grid_input.y_min, grid_input.z_min, grid_input.x_max, grid_input.y_max, grid_input.z_max, ghost_width, is_scalar_input, is_vector_input);
	}

public: // 	Operator Overloading
	inline TT& operator()(const VI& ix) const
	{
		return array_for_this(ix);
	}

	inline TT& operator()(const int& ix) const
	{
		return array_for_this(ix);
	}

	inline TT& operator()(const int& i, const int& j, const int& k) const
	{
		return array_for_this(i, j, k);
	}

	// Performance Check
	inline TT operator()(const VT& position) const
	{
		return TriLinearInterpolation(position);
	}

public: // Indexing Functions
	inline int ClampI(const int& i) const
	{
		if(i < i_start_g)
			return i_start_g;
		else if(i > i_end_g)
			return i_end_g;

		return i;
	}

	inline int ClampJ(const int& j) const
	{
		if(j < j_start_g)
			return j_start_g;
		else if(j > j_end_g)
			return j_end_g;

		return j;
	}

	inline int ClampK(const int& k) const
	{
		if(k < k_start_g)
			return k_start_g;
		else if(k > k_end_g)
			return k_end_g;

		return k;
	}

	inline const int Index1D(const int& i, const int& j, const int& k) const
	{
		return array_for_this.Index1D(i, j, k);
	}

	inline const int Index1D(const VI& index) const
	{
		return array_for_this.Index1D(index);
	}

public: // Member Functions
	inline VT CellCenter(const int& i, const int& j, const int& k) const
	{
		return grid.CellCenter(i, j, k);
	}

	inline TT TriLinearInterpolation(const VT& pos) const
	{
		// TODO: Clamp posx, posy, posz instead of clamping indices
		const T &posx(pos.x), &posy(pos.y), &posz(pos.z);
		
		const int i0 = (int)floor(((posx - x_min)*one_over_dx - (T)0.5));
		const int j0 = (int)floor(((posy - y_min)*one_over_dy - (T)0.5));
		const int k0 = (int)floor(((posz - z_min)*one_over_dz - (T)0.5));

		const T a = (posx - (x_min + ((T)0.5 + (T)(i0 - i_start))*dx))*one_over_dx;
		const T b = (posy - (y_min + ((T)0.5 + (T)(j0 - j_start))*dy))*one_over_dy;
		const T c = (posz - (z_min + ((T)0.5 + (T)(k0 - k_start))*dz))*one_over_dz;

		const int ci0(ClampI(i0)), cj0(ClampJ(j0)), ck0(ClampK(k0));
		const int ci0p1(ClampI(i0+1)), cj0p1(ClampJ(j0+1)), ck0p1(ClampK(k0+1));

		const TT v000(array_for_this(ci0, cj0, ck0)), v100(array_for_this(ci0p1, cj0, ck0)), v010(array_for_this(ci0, cj0p1, ck0)),
				 v001(array_for_this(ci0, cj0, ck0p1)), v110(array_for_this(ci0p1, cj0p1, ck0)), v101(array_for_this(ci0p1, cj0, ck0p1)),
				 v011(array_for_this(ci0, cj0p1, ck0p1)), v111(array_for_this(ci0p1, cj0p1, ck0p1));

		const T mamb = ((T)1 - a)*((T)1 - b);
		const T bc = b*c;

		return (v000*(mamb*((T)1-c)) + (a*((T)1-b)*((T)1-c))*v100 + (((T)1-a)*b*((T)1-c))*v010 + (mamb*c)*v001
			   + (a*b*((T)1-c))*v110 + (a*((T)1-b)*c)*v101 + (((T)1-a)*bc)*v011 + (a*bc)*v111);
	}

	inline T ComponentwiseMin(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7)
	{
		return MIN8(v0, v1, v2, v3, v4, v5, v6, v7);
	}

	inline T ComponentwiseMax(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7)
	{
		return MAX8(v0, v1, v2, v3, v4, v5, v6, v7);
	}

	inline VT ComponentwiseMin(const VT& v0, const VT& v1, const VT& v2, const VT& v3, const VT& v4, const VT& v5, const VT& v6, const VT& v7)
	{
		return VT(MIN8(v0.x, v1.x, v2.x, v3.x, v4.x, v5.x, v6.x, v7.x), MIN8(v0.y, v1.y, v2.y, v3.y, v4.y, v5.y, v6.y, v7.y), MIN8(v0.z, v1.z, v2.z, v3.z, v4.z, v5.z, v6.z, v7.z));
	}

	inline VT ComponentwiseMax(const VT& v0, const VT& v1, const VT& v2, const VT& v3, const VT& v4, const VT& v5, const VT& v6, const VT& v7)
	{
		return VT(MAX8(v0.x, v1.x, v2.x, v3.x, v4.x, v5.x, v6.x, v7.x), MAX8(v0.y, v1.y, v2.y, v3.y, v4.y, v5.y, v6.y, v7.y), MAX8(v0.z, v1.z, v2.z, v3.z, v4.z, v5.z, v6.z, v7.z));
	}

	inline void PrepareTriLinearInterpolation(const VT& position, bool& inside, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1,
											  T& w000, T& w001, T& w010, T& w011, T& w100, T& w101, T& w110, T& w111) 
	{
		T posx(position.x), posy(position.y), posz(position.z);

		i0 = (int)floor((posx - x_min)*one_over_dx - (T)0.5) + i_start;
		j0 = (int)floor((posy - y_min)*one_over_dy - (T)0.5) + j_start;
		k0 = (int)floor((posz - z_min)*one_over_dx - (T)0.5) + k_start;

		i1 = i0 + 1;
		j1 = j0 + 1;
		k1 = k0 + 1;

		if(grid.Inside(i0,j0,k0) == false || grid.Inside(i1,k1,k1) == false)
		{
			inside = false;
			return;
		}

		const T x0 = x_min + ((T)0.5 + (T)(i0 - i_start))*dx;
		const T y0 = y_min + ((T)0.5 + (T)(j0 - j_start))*dy;
		const T z0 = z_min + ((T)0.5 + (T)(k0 - k_start))*dz;

		const T a = (posx - x0)*one_over_dx;
		const T b = (posy - y0)*one_over_dy;
		const T c = (posz - z0)*one_over_dz;

		w000 = ((T)1 - a)*((T) - b)*((T)1 - c);
		w100 = a*((T)1 - b)*((T)1 - c);
		w010 = ((T)1 - a)*b*((T)1 - c);
		w001 = ((T)1 - a)*((T)1 - b)*c;
		w110 = a*b*((T)1 - c);
		w101 = a*((T)1 - b)*c;
		w011 = ((T)a - a)*b*c;
		w111 = a*b*c;

		assert(w000 >= (T)0 && w000 <= (T)1 && w100 >= (T)0 && w100 <= (T)1 && w010 >= (T)0 && w010 <= (T)1 && w001 >= (T)0 && w001 <= (T)1 &&
			   w110 >= (T)0 && w110 <= (T)1 && w101 >= (T)0 && w101 <= (T)1	&& w011 >= (T)0 && w011 <= (T)1 && w111 >= (T)0 && w111 <= (T)1);

		inside = true;
	}
	
	inline TT TriLinearInterpolation(const VT& position, const TT& componentwise_min, const TT& componentwise_max)
	{
		const T &posx(pos.x), &posy(pos.y), &posz(pos.z);
		
		const int i0 = (int)floor(((posx - x_min)*one_over_dx - (T)0.5));
		const int j0 = (int)floor(((posy - y_min)*one_over_dy - (T)0.5));
		const int k0 = (int)floor(((posz - z_min)*one_over_dz - (T)0.5));

		const T a = (posx - (x_min + ((T)0.5 + (T)(i0 - i_start))*dx))*one_over_dx;
		const T b = (posy - (y_min + ((T)0.5 + (T)(j0 - j_start))*dy))*one_over_dy;
		const T c = (posz - (z_min + ((T)0.5 + (T)(k0 - k_start))*dz))*one_over_dz;

		const TT v000(ClampedArrayValue(i0, j0, k0)), v100(ClampedArrayValue(i0+1, j0, k0)), v010(ClampedArrayValue(i0, j0+1, k0)),
				 v001(ClampedArrayValue(i0, j0, k0+1)), v110(ClampedArrayValue(i0+1, j0+1, k0)), v101(ClampedArrayValue(i0+1, j0, k0+1)),
				 v011(ClampedArrayValue(i0, j0+1, k0+1)), v111(ClampedArrayValue(i0+1, j0+1, k0+1));

		const T mamb = ((T)1 - a)*((T)1 - b);
		const T bc = b*c;

		return (v000*(mamb*((T)1-c)) + (a*((T)1-b)*((T)1-c))*v100 + (((T)1-a)*b*((T)1-c))*v010 + (mamb*c)*v001
			   + (a*b*((T)1-c))*v110 + (a*((T)1-b)*c)*v101 + (((T)1-a)*bc)*v011 + (a*bc)*v111);
	}

	void MakePartialGrids(const int& num_threads)
	{
		grid.SplitInMaxDirection(num_threads, partial_grids);
	}

	void AssignValuesInsideObject(const LEVELSET_OBJECT* obj, const TT& value)
	{
		GRID_ITERATION_3D(grid)
		{
			if((*obj)(CellCenter(i, j, k)) <= (T)0) 
			{
				array_for_this(i, j, k) = value;
			}
		}
	}

	void AssignValuesInsideObject(const int& thread_id, const LEVELSET_OBJECT* obj, const TT& value)
	{
		GRID_ITERATION_3D(grid)
		{
			if((*obj)(CellCenter(i, j, k)) <= (T)0) 
			{
				array_for_this(i, j, k) = value;
			}
		}
	}
	
	void AssignAllValues(const TT& value)
	{
		GRID_ITERATION_3D(grid)
		{
			array_for_this(i, j, k) = value;
		}
	}

	void AssignValues(const int& i, const int& j, const int& k, const TT& value) 
	{
		array_for_this(i, j, k) = value;
	}

	void AssignAllValues(const int& thread_id, const TT& value)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			array_for_this(i, j, k) = value;
		}
		END_GRID_ITERATION_3D;
	}

	void AssignAllValuesZeroGhost()
	{
		memset(array_for_this.values, 0, (array_for_this.ijk_res)*sizeof(TT));
	}

	void AssignAllValuesZeroGhost(const int& thread_id)
	{
		PREPARE_FOR_1D_ITERATION(array_for_this.ijk_res);
		std::memset(array_for_this.values + multithreading->start_ix_1D[thread_id], 0, (multithreading->end_ix_1D[thread_id] - multithreading->start_ix_1D[thread_id] + 1)*sizeof(TT));
	}

	void AssignAllValuesGhost(const int& thread_id, const TT& value)
	{
		PREPARE_FOR_1D_ITERATION(array_for_this.ijk_res);

		BEGIN_1D_ITERATION
		{
			array_for_this.values[p] = value;
		}
		END_1D_ITERATION
	}

	void AddAllValues(const int& thread_id, const TT& value)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			array_for_this(i, j, k) += value;
		}
		END_GRID_ITERATION_3D;
	}

	template<class TTT>
	void AssignAllValues(const int& thread_id, const ARRAY_3D<TTT>& arr, const TTT& value)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			arr(i, j, k) = value;
		}
		END_GRID_ITERATION_3D;
	}

	void CopyAllValuesFrom(const int& thread_id, const FIELD_STRUCTURE_3D& field_input)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			array_for_this(i, j, k) = field_input.array_for_this(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void CopyAllValuesGhostFrom(const int& thread_id, const FIELD_STRUCTURE_3D& field_input)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids_ghost[thread_id])
		{
			array_for_this(i, j, k) = field_input.array_for_this(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	int FindPartialGridK(const int& k)	// Use "k" only since we are generating partial grids by depth splitting
	{
		const int num_partial_grids(partial_grids.num_elements);
		for(int i = 0; i < num_partial_grids; i++)
		{
			if(partial_grids[i].InsideK(k)) return i;
		}
		
		assert(false);

		return 0;
	}

	int FindPartialGridGhostK(const int& k)
	{
		const int num_partial_grids(partial_grids_ghost.num_elements);
		for(int i = 0; i < num_partial_grids; i++)
		{
			if(partial_grids_ghost[i].InsideK(k)) return i;
		}

		assert(false);
	
		return 0;
	}

	static void Sampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_3D<T>& from_input, FIELD_STRUCTURE_3D<T>& to_input)
	{
		// Note : Only levelset upsampling handles ghost cells, whereas FIELD_STRUCTURE_3D doesn't
		ARRAY_3D<T>& to_array(to_input.array_for_this);

		GRID_ITERATION_3D(to_input.partial_grids[thread_id])
		{
			to_array(i,j,k) = from_input.TriLinearInterpolation(to_input.CellCenter(i,j,k));
		}
		multithreading->Sync(thread_id);
	}
	
	static void Sampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_3D<VT>& from_input, FIELD_STRUCTURE_3D<VT>& to_input)
	{
		ARRAY_3D<VT>& to_array(to_input.array_for_this);

		GRID_ITERATION_3D(to_input.partial_grids[thread_id])
		{
			to_array(i,j,k) = from_input.TriLinearInterpolation(to_input.CellCenter(i,j,k));
		}
		multithreading->Sync(thread_id);
	}
	
	void Read(const char* file_name, MULTITHREADING* multithreading_input);
	void Write(const char* file_name);

public: // Speedup Functions
	inline const TT ArrayValue(const int& i, const int& j, const int& k) const
	{
		return *(values + (i - i_start_g) + (j - j_start_g)*i_res_g + (k - k_start_g)*ij_res_g);
	}

	inline const VI ClampedArrayIndex(const int& i, const int& j, const int& k) const
	{
		return VI(CLAMP(i,i_start_g,i_end_g), CLAMP(j,j_start_g,j_end_g), CLAMP(k,k_start_g,k_end_g));
	}

	inline const VI ClampedArrayIndex(const VI& ix) const
	{
		return VI(CLAMP(ix.i,i_start_g,i_end_g), CLAMP(ix.j,j_start_g,j_end_g), CLAMP(ix.k,k_start_g,k_end_g));
	}

	inline const TT ClampedArrayValue(const VI& ix) const
	{
		return values[(CLAMP(ix.i,i_start_g,i_end_g) - i_start_g) + (CLAMP(ix.j,j_start_g,j_end_g) - j_start_g)*i_res_g + (CLAMP(ix.k,k_start_g,k_end_g) - k_start_g)*ij_res_g];
	}

	inline const TT ClampedArrayValue(const int& i, const int& j, const int& k) const
	{
		return *(values + ((CLAMP(ix.i,i_start_g,i_end_g) - i_start_g) + (CLAMP(ix.j,j_start_g,j_end_g) - j_start_g)*i_res_g + (CLAMP(ix.k,k_start_g,k_end_g) - k_start_g)*ij_res_g));
	}

	inline void LeftBottomCell(const VT& position, int& i, int& j, int& k) const
	{
		i = (int)((position.x - x_min)*one_over_dx - (T)0.5) + i_start;
		j = (int)((position.y - y_min)*one_over_dy - (T)0.5) + j_start;
		k = (int)((position.z - z_min)*one_over_dz - (T)0.5) + k_start;
	}

	inline void LeftBottomCell(const VT& position, VI& ix) const
	{
		ix.i = (int)((position.x - x_min)*one_over_dx - (T)0.5) + i_start;
		ix.j = (int)((position.y - y_min)*one_over_dy - (T)0.5) + j_start;
		ix.k = (int)((position.z - z_min)*one_over_dz - (T)0.5) + k_start;
	}

	// Need to fill this out according to the number of your system's cores
	void FillGhostCellsZero(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
					  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);

		// Face left
		if(0 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = (TT)0;
					}
				}
			}
		}

		// Face right
		if(1 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = (TT)0;
					}
				}
			}
		}

		// Face bottom
		if(2 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (TT)0;
					}
				}
			}
		}

		// Face up
		if(3 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (TT)0;
					}
				}
			}
		}

		// Face front
		if(4 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (TT)0;
					}
				}
			}
		}

		// Face back
		if(5 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (TT)0;
					}
				}
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFrom(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
					  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);
		
		// Face left
		if(0 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j, k);
					}
				}
			}
		}

		// Face right
		if(1 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j, k);
					}
				}
			}
		}

		// Face bottom
		if(2 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j_start, k);
					}
				}
			}
		}

		// Face up
		if(3 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j_end, k);
					}
				}
			}
		}

		// Face front
		if(4 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k_start);
					}
				}
			}
		}

		// Face back
		if(5 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k_end);
					}
				}
			}
		}

		// Edge 1 normal to z plane
		if(6 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j_start, k);
					}
				}
			}
		}

		// Edge 2 normal to z plane
		if(7 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j_end, k);
					}
				}
			}
		}

		// Edge 3 normal to z plane
		if(8 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j_end, k);
					}
				}
			}
		}

		// Edge 4 normal to z plane
		if(9 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j_start, k);
					}
				}
			}
		}

		// Edge 1 normal to y plane
		if(10 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j, k_start);
					}
				}
			}
		}

		// Edge 2 normal to y plane
		if(11 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j, k_end);
					}
				}
			}
		}

		// edge 3 normal to y plane
		if(12 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j, k_start);
					}
				}
			}
		}

		// edge 4 normal to y plane
		if(13 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j, k_end);
					}
				}
			}
		}

		// edge 1 normal to x plane
		if(14 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j_start, k_start);
					}
				}
			}
		}

		// edge 2 normal to x plane
		if(15 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j_end, k_end);
					}
				}
			}
		}

		// edge 3 normal to x plane
		if(16 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j_start, k_end);
					}
				}
			}
		}

		// edge 4 normal to x plane
		if(17 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j_end, k_start);
					}
				}
			}
		}
	
		// vertex 1
		if(18 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j_start, k_start);
					}
				}
			}
		}

		// vertex 2
		if(19 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j_end, k_end);
					}
				}
			}
		}
		
		// vertex 3
		if(20 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j_start, k_start);
					}
				}
			}
		}

		// vertex 4
		if(21 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i,j,k) = phi_real(i_start, j_end, k_start);
					}
				}
			}
		}

		// vertex 5
		if(22 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j_start, k_end);
					}
				}
			}
		}

		// vertex 6
		if(23 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = phi_real(i_end, j_end, k_start);
					}
				}
			}
		}

		// vertex 7
		if(24 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i,j,k) = phi_real(i_end, j_start, k_end);
					}
				}
			}
		}

		// vertex 8
		if(25 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = phi_real(i_start, j_end, k_end);
					}
				}
			}
		}

		multithreading->Sync(thread_id);
	}
	
	void FillGhostCellsPeriodicInXDirection(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
					  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);
		
		if (is_x_component)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_start; i >= i_start - ghost_width; i--)
						{
							array_for_this(i, j, k) = phi_real(i_end + (i - (i_start + 1)), j, k);
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_end; i <= i_end + ghost_width; i++)
						{
							array_for_this(i, j, k) = phi_real(i_start + (i - (i_end - 1)), j, k);
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
		if (is_y_component || is_z_component || is_scalar)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_start; i >= i_start - ghost_width; i--)
						{
							array_for_this(i, j, k) = phi_real(i_end + (i - i_start), j, k);
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_end + 1; i <= i_end + ghost_width; i++)
						{
							array_for_this(i, j, k) = phi_real(i_start + (i - i_end), j, k); 
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
	}

	void FillGhostCellsPeriodicInXDirection(const int& thread_id, ARRAY_3D<TT>& phi_real, ARRAY_3D<TT>& boundary_levelset, const bool& copy_real_data)
	{
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
				  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						if (boundary_levelset(i, j, k) <= 0)
						{
							array_for_this(i, j, k) = phi_real(i, j, k);
						}
						else
						{
							continue;
						}
					}
				}
			}
			multithreading->Sync(thread_id);
		}
		
		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);
		
		if (is_x_component)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_start; i >= i_start - ghost_width; i--)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i_end + (i - (i_start + 1)), j, k);
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_end; i <= i_end + ghost_width; i++)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i_start + (i - (i_end - 1)), j, k);
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
		if (is_y_component || is_z_component || is_scalar)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_start; i >= i_start - ghost_width; i--)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i_end + (i - i_start), j, k);
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(j = j_start; j <= j_end; j++)
					{
						for(i = i_end + 1; i <= i_end + ghost_width; i++)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i_start + (i - i_end), j, k); 
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
	}

	void FillGhostCellsPeriodicInYDirection(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		int i, j, k;
				
		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
				  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
			multithreading->Sync(thread_id);
		}
		
		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int num_threads(partial_grids.num_elements);
		
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);

		if (is_y_component)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_start; j >= j_start - ghost_width; j--)
						{
							array_for_this(i, j, k) = phi_real(i, j_end + (j - (j_start + 1)), k);
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_end; j <= j_end + ghost_width; j++)
						{
							array_for_this(i, j, k) = phi_real(i, j_start + (j - (j_end - 1)), k);
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
		if (is_x_component || is_z_component || is_scalar)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_start - 1; j >= j_start - ghost_width; j--)
						{
							array_for_this(i, j, k) = phi_real(i, j_end + (j - j_start), k);
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_end + 1; j <= j_end + ghost_width; j++)
						{
							array_for_this(i, j, k) = phi_real(i, j_start + (j - j_end), k); 
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
	}

	void FillGhostCellsPeriodicInYDirection(const int& thread_id, ARRAY_3D<TT>& phi_real, ARRAY_3D<TT>& boundary_levelset, const bool& copy_real_data)
	{
		int i, j, k;
				
		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
				  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						if (boundary_levelset(i, j, k) < (double)1e-8)
						{
							array_for_this(i, j, k) = phi_real(i, j, k);
						}
						else
						{
							array_for_this(i, j, k) = 0;
						}
					}
				}
			}
			multithreading->Sync(thread_id);
		}
		
		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int num_threads(partial_grids.num_elements);
		
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);

		if (is_y_component)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_start; j >= j_start - ghost_width; j--)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i, j_end + (j - (j_start + 1)), k);
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_end; j <= j_end + ghost_width; j++)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i, j_start + (j - (j_end - 1)), k);
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}
		if (is_x_component || is_z_component || is_scalar)
		{
			// i_start
			if(0 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_start - 1; j >= j_start - ghost_width; j--)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i, j_end + (j - j_start), k);
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			// i_end
			if(1 % num_threads == thread_id)
			{
				for(k = k_start; k <= k_end; k++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						for(j = j_end + 1; j <= j_end + ghost_width; j++)
						{
							if (boundary_levelset(i, j, k) < (double)1e-8)
							{
								array_for_this(i, j, k) = phi_real(i, j_start + (j - j_end), k); 
							}
							else
							{
								continue;
							}
						}
					}
				}
			}

			multithreading->Sync(thread_id); 
		}

		if(2 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						if (boundary_levelset(i - 1, j, k) < 0 && boundary_levelset(i, j, k) > 0)
						{
							for (int t = i; t <= i_end + ghost_width; t++)
							{
								array_for_this(t, j, k) = (T)2*phi_real(t - 1, j, k) - phi_real(t - 2, j, k);
							}
						}
						else
						{
							continue;	
						}
					}
				}
			}
		}

		if(3 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						if (boundary_levelset(i + 1, j, k) < 0 && boundary_levelset(i, j, k) > 0)
						{
							for (int t = i; t >= i_start - ghost_width; t--)
							{
								array_for_this(t, j, k) = (T)2*phi_real(t + 1, j, k) - phi_real(t + 2, j, k);
							}
						}
						else
						{
							continue;	
						}
					}
				}
			}
		}

		if(4 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						if (boundary_levelset(i, j, k - 1) < 0 && boundary_levelset(i, j, k) > 0)
						{
							for (int t = k; t <= k_end + ghost_width; t++)
							{
								array_for_this(i, j, t) = (T)2*phi_real(i, j, t - 1) - phi_real(i, j, t - 2);
							}
						}
						else
						{
							continue;	
						}
					}
				}
			}
		}

		if(4 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						if (boundary_levelset(i, j, k + 1) < 0 && boundary_levelset(i, j, k) > 0)
						{
							for (int t = k; t >= k_start - ghost_width; t--)
							{
								array_for_this(i, j, t) = (T)2*phi_real(i, j, t + 1) - phi_real(i, j, t + 2);
							}
						}
						else
						{
							continue;	
						}
					}
				}
			}
		}
		//// Face right
		//if(1 % num_threads == thread_id)
		//{
		//	for(k = k_start; k <= k_end; k++)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k);
		//			}
		//		}
		//	}
		//}
		//
		//// Face front
		//if(2 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_start; i <= i_end; i++)
		//			{
		//				array_for_this(i, j, k) = (T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2);
		//			}
		//		}
		//	}
		//}

		//// Face back
		//if(3 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_start; i <= i_end; i++)
		//			{
		//				array_for_this(i, j, k) = (T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2);
		//			}
		//		}
		//	}
		//}
	}

	void FillGhostCellsPeriodicInZDirection(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
					  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);
		
		// i_start
		if(0 % num_threads == thread_id)
		{
			for(j = j_start; j <= j_end; j++)
			{
				for(i = i_start; i <= i_end; i++)
				{
					for(k = k_start - 1; k >= k_start - ghost_width; k--)
					{
						array_for_this(i, j, k) = phi_real(i, j, k_end + (k - k_start));
					}
				}
			}
		}

		// i_end
		if(1 % num_threads == thread_id)
		{
			for(j = j_start; j <= j_end; j++)
			{
				for(i = i_start; i <= i_end; i++)
				{
					for(k = k_end + 1; k <= k_end + ghost_width; k++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k_start + (k - k_end));
					}
				}
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsContinuousDerivativesFrom(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
					  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);
		
		// Face left
		if(0 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - 1; i >= i_start - ghost_width; i--)
					{
						array_for_this(i, j, k) = (T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k);
					}
				}
			}
		}

		// Face right
		if(1 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = (T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k);
					}
				}
			}
		}

		// Face bottom
		if(2 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - 1; j >= j_start - ghost_width; j--)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k);
					}
				}
			}
		}

		// Face up
		if(3 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k);
					}
				}
			}
		}

		// Face front
		if(4 % num_threads == thread_id)
		{
			for(k = k_start - 1; k >= k_start - ghost_width; k--)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2);
					}
				}
			}
		}

		// Face back
		if(5 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = (T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2);
					}
				}
			}
		}

		// Edge 1 normal to z plane
		if(6 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - 1; j >= j_start - ghost_width; j--)
				{
					for(i = i_start - 1; i >= i_start - ghost_width; i--)
					{
						array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)));
					}
				}
			}
		}

		// Edge 2 normal to z plane
		if(7 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)));
					}
				}
			}
		}

		// Edge 3 normal to z plane
		if(8 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - 1; i >= i_start - ghost_width; i--)
					{
						array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)));
					}
				}
			}
		}

		//// Edge 4 normal to z plane
		//if(9 % num_threads == thread_id)
		//{
		//	for(k = k_start; k <= k_end; k++)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)));
		//			}
		//		}
		//	}
		//}

		//// Edge 1 normal to y plane
		//if(10 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_start - 1; i >= i_start - ghost_width; i--)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// Edge 2 normal to y plane
		//if(11 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		//// edge 3 normal to y plane
		//if(12 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// edge 4 normal to y plane
		//if(13 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_start; j <= j_end; j++)
		//		{
		//			for(i = i_start - 1; i >= i_start - ghost_width; i--)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		//// edge 1 normal to x plane
		//if(14 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_start; i <= i_end; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// edge 2 normal to x plane
		//if(15 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_end + 1; j <= j_end + ghost_width; j++)
		//		{
		//			for(i = i_start; i <= i_end; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// edge 3 normal to x plane
		//if(16 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_start; i <= i_end; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		//// edge 4 normal to x plane
		//if(17 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_end + 1; j <= j_end + ghost_width; j++)
		//		{
		//			for(i = i_start; i <= i_end; i++)
		//			{
		//				array_for_this(i, j, k) = (T)0.5*(((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}
	
		//// vertex 1
		//if(18 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_start - 1; i >= i_start - ghost_width; i--)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// vertex 2
		//if(19 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_end + 1; j <= j_end + ghost_width; j++)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}
		//
		//// vertex 3
		//if(20 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// vertex 4
		//if(21 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_end + 1; j <= j_end + ghost_width; j++)
		//		{
		//			for(i = i_start - 1; i >= i_start - ghost_width; i--)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)) + ((T)2*phi_real(i, j, k + 1) - phi_real(i, j, k + 2)));
		//			}
		//		}
		//	}
		//}

		//// vertex 5
		//if(22 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_start - 1; i >= i_start - ghost_width; i--)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		//// vertex 6
		//if(23 % num_threads == thread_id)
		//{
		//	for(k = k_start - 1; k >= k_start - ghost_width; k--)
		//	{
		//		for(j = j_end + 1; j <= j_end + ghost_width; j++)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		//// vertex 7
		//if(24 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_start - 1; j >= j_start - ghost_width; j--)
		//		{
		//			for(i = i_end + 1; i <= i_end + ghost_width; i++)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i - 1, j, k) - phi_real(i - 2, j, k)) + ((T)2*phi_real(i, j + 1, k) - phi_real(i, j + 2, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		//// vertex 8
		//if(25 % num_threads == thread_id)
		//{
		//	for(k = k_end + 1; k <= k_end + ghost_width; k++)
		//	{
		//		for(j = j_end + 1; j <= j_end + ghost_width; j++)
		//		{
		//			for(i = i_start - 1; i >= i_start - ghost_width; i--)
		//			{
		//				array_for_this(i, j, k) = (T)1/(T)3*(((T)2*phi_real(i + 1, j, k) - phi_real(i + 2, j, k)) + ((T)2*phi_real(i, j - 1, k) - phi_real(i, j - 2, k)) + ((T)2*phi_real(i, j, k - 1) - phi_real(i, j, k - 2)));
		//			}
		//		}
		//	}
		//}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFromDirichlet(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{
		// This function fills Ghost Cells with Zero values
		int i, j, k;

		// Fill real region
		if(copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), k_start(partial_grids[thread_id].k_start),
					  i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end), k_end(partial_grids[thread_id].k_end);

			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = phi_real(i, j, k);
					}
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill ghost region
		const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);
		const int num_threads(partial_grids.num_elements);
		
		// Face left
		if(0 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Face right
		if(1 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Face bottom
		if(2 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Face up
		if(3 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Face front
		if(4 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Face back
		if(5 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Edge 1 normal to z plane
		if(6 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Edge 2 normal to z plane
		if(7 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Edge 3 normal to z plane
		if(8 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Edge 4 normal to z plane
		if(9 % num_threads == thread_id)
		{
			for(k = k_start; k <= k_end; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Edge 1 normal to y plane
		if(10 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// Edge 2 normal to y plane
		if(11 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// edge 3 normal to y plane
		if(12 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// edge 4 normal to y plane
		if(13 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start; j <= j_end; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// edge 1 normal to x plane
		if(14 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// edge 2 normal to x plane
		if(15 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// edge 3 normal to x plane
		if(16 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// edge 4 normal to x plane
		if(17 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start; i <= i_end; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}
	
		// vertex 1
		if(18 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// vertex 2
		if(19 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}
		
		// vertex 3
		if(20 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// vertex 4
		if(21 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i,j,k) = TT();
					}
				}
			}
		}

		// vertex 5
		if(22 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// vertex 6
		if(23 % num_threads == thread_id)
		{
			for(k = k_start - ghost_width; k <= k_start - 1; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		// vertex 7
		if(24 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					for(i = i_end + 1; i <= i_end + ghost_width; i++)
					{
						array_for_this(i,j,k) = TT();
					}
				}
			}
		}

		// vertex 8
		if(25 % num_threads == thread_id)
		{
			for(k = k_end + 1; k <= k_end + ghost_width; k++)
			{
				for(j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					for(i = i_start - ghost_width; i <= i_start - 1; i++)
					{
						array_for_this(i, j, k) = TT();
					}
				}
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFromFreeSlip(const int& thread_id, ARRAY_3D<TT>& phi_real, const bool& copy_real_data)
	{}

	void PrintKPlaneValues(const int& k, const bool& print_ghost);

	void Translate(const VT& deviation)
	{
		grid.Translate(deviation);
		grid_ghost.Translate(deviation);

		for(int i = 0; i < multithreading->num_threads; i++)
		{
			partial_grids[i].Translate(deviation);
			partial_grids_ghost[i].Translate(deviation);
		}

		// Translate SpeedUp variables
		x_min += deviation.x;
		y_min += deviation.y;
		z_min += deviation.z;
	}

	/*void GetMaxValue(const int& thread_id)
	{
		max_abs_value = 0;

		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			if (max_abs_value <= abs(array_for_this(i, j, k)))
			{
				max_abs_value = abs(array_for_this(i, j, k));
			}
		}
		END_GRID_ITERATION_3D;
	}*/
};

static void Sampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_3D<T>& from_input, FIELD_STRUCTURE_3D<T>& to_input)
{
	ARRAY_3D<T>& to_arr(to_input.array_for_this);

	GRID_ITERATION_3D(to_input.partial_grids[thread_id])
	{
		to_arr(i, j, k) = from_input.TriLinearInterpolation(to_input.CellCenter(i, j, k));
	}
	multithreading->Sync(thread_id);
}

static void Sampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_3D<VT>& from_input, FIELD_STRUCTURE_3D<VT>& to_input)
{
	ARRAY_3D<VT>& to_arr(to_input.array_for_this);

	GRID_ITERATION_3D(to_input.partial_grids[thread_id])
	{
		to_arr(i, j, k) = from_input.TriLinearInterpolation(to_input.CellCenter(i, j, k));
	}
	multithreading->Sync(thread_id);
}

static void DownSampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_3D<T>& from_input, FIELD_STRUCTURE_3D<T>& to_input)
{
	const ARRAY_3D<T>& from_arr(from_input.array_for_this);
	ARRAY_3D<T>& to_arr(to_input.array_for_this);

	GRID_ITERATION_3D(to_input.partial_grids[thread_id])
	{
		to_arr(i,j,k) = ( from_arr(i*2,   j*2,   k*2  ) 
						+ from_arr(i*2,   j*2,   k*2+1) 
						+ from_arr(i*2,   j*2+1, k*2  )
						+ from_arr(i*2,   j*2+1, k*2+1)
						+ from_arr(i*2+1, j*2,   k*2  )
						+ from_arr(i*2+1, j*2,   k*2+1)
						+ from_arr(i*2+1, j*2+1, k*2  )
						+ from_arr(i*2+1, j*2+1, k*2+1) ) / (T)8;		
	}	
	multithreading->Sync(thread_id);
}

static void DownSampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_3D<VT>& from_input, FIELD_STRUCTURE_3D<VT>& to_input)
{
	const ARRAY_3D<VT>& from_arr(from_input.array_for_this);
	ARRAY_3D<VT>& to_arr(to_input.array_for_this);

	GRID_ITERATION_3D(to_input.partial_grids[thread_id])
	{
		to_arr(i,j,k) = ( from_arr(i*2,   j*2,   k*2  ) 
						+ from_arr(i*2,   j*2,   k*2+1) 
						+ from_arr(i*2,   j*2+1, k*2  )
						+ from_arr(i*2,   j*2+1, k*2+1)
						+ from_arr(i*2+1, j*2,   k*2  )
						+ from_arr(i*2+1, j*2,   k*2+1)
						+ from_arr(i*2+1, j*2+1, k*2  )
						+ from_arr(i*2+1, j*2+1, k*2+1) ) / (T)8;		
	}	
	multithreading->Sync(thread_id);
}






		










		


