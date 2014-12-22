#pragma once

#include "GRID_STRUCTURE_2D.h"
#include "ARRAY_2D.h"
#include "MULTITHREADING.h"
#include "LEVELSET_OBJECT.h"
#include "GS_COEFFICIENTS.h"

template<class TT>
class FIELD_STRUCTURE_2D
{
public: // Essential Data
	GRID_STRUCTURE_2D grid;
	GRID_STRUCTURE_2D grid_ghost;								// Grid including ghost cells. Same as grid when ghost_width = 0
	ARRAY<GRID_STRUCTURE_2D> partial_grids;						// Partial grid for multithreading
	ARRAY<GRID_STRUCTURE_2D> partial_grids_ghost;				// Partial grid for grid_ghost
	ARRAY_2D<TT> array_for_this;

public: // Properties
	int ghost_width;
	
public: // Multithreading
	MULTITHREADING* multithreading;

public: // Speed Up Constants and Variables
	int i_start, j_start, i_start_g, j_start_g;
	int i_end, j_end, i_end_g, j_end_g;
	int i_res_g, ij_res_g;
	T x_min, y_min;
	T dx, dy;
	T one_over_dx, one_over_dy;
	T one_over_2dx, one_over_2dy;
	T one_over_dx2, one_over_dy2;
	TT* values;

public: // Constructors and Destructor
	FIELD_STRUCTURE_2D(void)
		: multithreading(0)
	{}

	FIELD_STRUCTURE_2D(MULTITHREADING* multithreading_input, const VI& res_input, const VI& start_input, const VT& min_input, const VT& max_input, const int& ghost_width_input = 0)
		: partial_grids(0)
	{
		Initialize(multithreading_input, res_input, start_input, min_input, max_input, ghost_width);
	}

	~FIELD_STRUCTURE_2D(void)
	{}

public: // Initialization Functions
	void Initialize(MULTITHREADING* multithreading_input, const int& i_res_input, const int& j_res_input, const int& i_start_input, const int& j_start_input, 
					const T& x_min_input, const T& y_min_input, const T& x_max_input, const T& y_max_input, const int& ghost_width_input = 0)
	{
		grid.Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input);

		ghost_width = ghost_width_input;

		grid_ghost.Initialize(grid.Enlarged(ghost_width));

		i_start = grid.i_start;
		j_start = grid.j_start;
		
		i_end = grid.i_end;
		j_end = grid.j_end;
		
		i_start_g = grid_ghost.i_start;
		j_start_g = grid_ghost.j_start;
		
		i_end_g = grid_ghost.i_end;
		j_end_g = grid_ghost.j_end;
		
		x_min = grid.x_min;
		y_min = grid.y_min;
		
		dx = grid.dx;
		dy = grid.dy;
		
		one_over_dx = grid.one_over_dx;
		one_over_dy = grid.one_over_dy;
		
		one_over_2dx = grid.one_over_2dx;
		one_over_2dy = grid.one_over_2dy;
		
		one_over_dx2 = grid.one_over_dx2;
		one_over_dy2 = grid.one_over_dy2;
		
		// Let array_for_this have the information for grid_ghost
		array_for_this.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, true);

		values = array_for_this.values;

		i_res_g = array_for_this.i_res;
		ij_res_g = array_for_this.ij_res;

		// Multithreading
		assert(multithreading_input != 0);
		multithreading = multithreading_input;
	}
	
	void Initialize(MULTITHREADING* multithreading_input, const VI& res_input, const VI& start_input, const VT& min_input, const VT& max_input, const int& ghost_width_input = 0)
	{
		Initialize(multithreading_input, res_input.i, res_input.j, start_input.i, start_input.j, min_input.i, min_input.j, max_input.i, max_input.j, ghost_width_input);
	}

	void Initialize(MULTITHREADING* multithreading_input, const GRID_STRUCTURE_2D& grid_input, const int& ghost_width = 0)
	{
		Initialize(multithreading_input, grid_input.i_res, grid_input.j_res, grid_input.i_start, grid_input.j_start, grid_input.x_min, grid_input.y_min, grid_input.x_max, grid_input.y_max, ghost_width);
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

	inline TT& operator()(const int& i, const int& j) const
	{
		return array_for_this(i, j);
	}

	// Performance Check
	inline TT operator()(const VT& position) const
	{
		return BiLinearInterpolation(position);
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

	inline const int Index1D(const int& i, const int& j) const
	{
		return array_for_this.Index1D(i, j);
	}

	inline const int Index1D(const VI& index) const
	{
		return array_for_this.Index1D(index);
	}

public: // Member Functions
	inline VT CellCenter(const int& i, const int& j) const
	{
		return grid.CellCenter(i, j);
	}

	inline TT BiLinearInterpolation(const VT& pos) const
	{
		// TODO: Clamp posx, posy, posz instead of clamping indices
		const T &posx(pos.x), &posy(pos.y);
		
		const int i0 = (int)floor(((posx - x_min)*one_over_dx - (T)0.5));
		const int j0 = (int)floor(((posy - y_min)*one_over_dy - (T)0.5));
		
		const T a = (posx - (x_min + ((T)0.5 + (T)(i0 - i_start))*dx))*one_over_dx;
		const T b = (posy - (y_min + ((T)0.5 + (T)(j0 - j_start))*dy))*one_over_dy;
		
		const int ci0(ClampI(i0)), cj0(ClampJ(j0));
		const int ci0p1(ClampI(i0+1)), cj0p1(ClampJ(j0+1));

		const TT v00(array_for_this(ci0,cj0)), v10(array_for_this(ci1,cj0)), v01(array_for_this(ci0,cj1)), v11(array_for_this(ci1,cj1));
		
		const T mamb = ((T)1 - a)*((T)1 - b);
		
		return mamb*v00 + a*((T)1 - b)*v10 + b*((T)1 - a)*v01 + b*a*v11;
	}

	inline void PrepareBiLinearInterpolation(const VT& position, bool& inside, int& i0, int& j0, int& i1, int& j1, T& w00, T& w01, T& w10, T& w11) 
	{
		T posx(position.x), posy(position.y);

		i0 = (int)floor((posx - x_min)*one_over_dx - (T)0.5) + i_start;
		j0 = (int)floor((posy - y_min)*one_over_dy - (T)0.5) + j_start;
		
		i1 = i0 + 1;
		j1 = j0 + 1;
		
		if(grid.Inside(i0,j0) == false || grid.Inside(i1,k1) == false)
		{
			inside = false;
			return;
		}

		const T x0 = x_min + ((T)0.5 + (T)(i0 - i_start))*dx;
		const T y0 = y_min + ((T)0.5 + (T)(j0 - j_start))*dy;
		
		const T a = (posx - x0)*one_over_dx;
		const T b = (posy - y0)*one_over_dy;
		
		w00 = ((T)1 - a)*((T) - b);
		w01 = b*((T)1 - a);
		w10 = ((T)1 - b)*a;
		w11 = a*b;

		assert(w00 >= 0 && w00 <= 1 && w01 >= 0 && w01 <= 1 && w10 >= 0 && w10 <= 1 && w11 >= 0 && w11 <= 1);

		inside = true;
	}
	
	inline TT BiLinearInterpolation(const VT& position, const TT& componentwise_min, const TT& componentwise_max)
	{
		const T &posx(pos.x), &posy(pos.y);
		
		const int i0 = (int)floor(((posx - x_min)*one_over_dx - (T)0.5));
		const int j0 = (int)floor(((posy - y_min)*one_over_dy - (T)0.5));
		
		const T a = (posx - (x_min + ((T)0.5 + (T)(i0 - i_start))*dx))*one_over_dx;
		const T b = (posy - (y_min + ((T)0.5 + (T)(j0 - j_start))*dy))*one_over_dy;
		
		const TT v00(ClampedArrayValue(i0,j0)), v01(ClampedArrayValue(i0,j0+1)), v10(ClampedArrayValue(i0+1,j0)), v11(ClampedArrayValue(i0+1,j0+1)); 

		const T mamb = ((T)1 - a)*((T)1 - b);
		
		return mamb*v00 + a*((T)1 - b)*v10 + b*((T)1 - a)*v01 + b*a*v11;
	}

	void MakePartialGrids(const int& num_threads)
	{
		grid.SplitInMaxDirection(num_threads, partial_grids);
	}

	void AssignValuesInsideObject(const LEVELSET_OBJECT* obj, const TT& value)
	{
		GRID_ITERATION_2D(grid)
		{
			if((*obj)(CellCenter(i, j)) <= (T)0) 
			{
				array_for_this(i, j) = value;
			}
		}
	}

	void AssignValuesInsideObject(const int& thread_id, const LEVELSET_OBJECT* obj, const TT& value)
	{
		GRID_ITERATION_2D(grid)
		{
			if((*obj)(CellCenter(i, j)) <= (T)0) 
			{
				array_for_this(i, j) = value;
			}
		}
	}
	
	void AssignAllValues(const TT& value)
	{
		GRID_ITERATION_2D(grid)
		{
			array_for_this(i, j) = value;
		}
	}

	void AssignAllValues(const int& thread_id, const TT& value)
	{
		BEGIN_GRID_ITERATION_2D(partial_grids[thread_id])
		{
			array_for_this(i, j) = value;
		}
		END_GRID_ITERATION_2D;
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
		BEGIN_GRID_ITERATION_2D(partial_grids[thread_id])
		{
			array_for_this(i, j, k) += value;
		}
		END_GRID_ITERATION_2D;
	}

	template<class TTT>
	void AssignAllValues(const int& thread_id, const ARRAY_2D<TTT>& arr, const TTT& value)
	{
		BEGIN_GRID_ITERATION_2D(partial_grids[thread_id])
		{
			arr(i, j) = value;
		}
		END_GRID_ITERATION_2D;
	}

	void CopyAllValuesFrom(const int& thread_id, const FIELD_STRUCTURE_2D& field_input)
	{
		BEGIN_GRID_ITERATION_2D(partial_grids[thread_id])
		{
			array_for_this(i, j) = field_input.array_for_this(i, j);
		}
		END_GRID_ITERATION_2D;
	}

	void CopyAllValuesGhostFrom(const int& thread_id, const FIELD_STRUCTURE_2D& field_input)
	{
		BEGIN_GRID_ITERATION_2D(partial_grid_ghost[thread_id])
		{
			array_for_this(i, j) = field_input.array_for_this(i, j);
		}
		END_GRID_ITERATION_2D;
	}

	static void Sampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_2D<T>& from_input, FIELD_STRUCTURE_2D<T>& to_input)
	{
		// Note : Only levelset upsampling handles ghost cells, whereas FIELD_STRUCTURE_2D doesn't
		ARRAY_2D<T>& to_array(to_input.array_for_this);

		GRID_ITERATION_2D(to_input.partial_grids[thread_id])
		{
			to_array(i,j) = from_input.TriLinearInterpolation(to_input.CellCenter(i,j));
		}
		multithreading->Sync(thread_id);
	}
	
	static void Sampling(MULTITHREADING* multithreading, const int& thread_id, const FIELD_STRUCTURE_2D<VT>& from_input, FIELD_STRUCTURE_2D<VT>& to_input)
	{
		ARRAY_2D<VT>& to_array(to_input.array_for_this);

		GRID_ITERATION_2D(to_input.partial_grids[thread_id])
		{
			to_array(i,j) = from_input.TriLinearInterpolation(to_input.CellCenter(i,j));
		}
		multithreading->Sync(thread_id);
	}
	
	void Read(const char* file_name, MULTITHREADING* multithreading_input);
	void Write(const char* file_name);

public: // Speedup Functions
	inline const TT ArrayValue(const int& i, const int& j) const
	{
		return *(values + (i - i_start_g) + (j - j_start_g)*i_res_g);
	}

	inline const VI ClampedArrayIndex(const int& i, const int& j) const
	{
		return VI(CLAMP(i,i_start_g,i_end_g), CLAMP(j,j_start_g,j_end_g), 0);
	}

	inline const VI ClampedArrayIndex(const VI& ix) const
	{
		return VI(CLAMP(ix.i,i_start_g,i_end_g), CLAMP(ix.j,j_start_g,j_end_g), 0);
	}

	inline const TT ClampedArrayValue(const VI& ix) const
	{
		return values[(CLAMP(ix.i,i_start_g,i_end_g) - i_start_g) + (CLAMP(ix.j,j_start_g,j_end_g) - j_start_g)*i_res_g];
	}

	inline const TT ClampedArrayValue(const int& i, const int& j) const
	{
		return *(values + ((CLAMP(ix.i,i_start_g,i_end_g) - i_start_g) + (CLAMP(ix.j,j_start_g,j_end_g) - j_start_g)*i_res_g));
	}

	inline void LeftBottomCell(const VT& position, int& i, int& j) const
	{
		i = (int)((position.x - x_min)*one_over_dx - (T)0.5) + i_start;
		j = (int)((position.y - y_min)*one_over_dy - (T)0.5) + j_start;
	}	
	
	inline void LeftBottomCell(const VT& position, VI& ix) const
	{
		ix.i = (int)((position.x - x_min)*one_over_dx - (T)0.5) + i_start;
		ix.j = (int)((position.y - y_min)*one_over_dy - (T)0.5) + j_start;
	}

	void FillGhostCellsFrom(const int& thread_id, ARRAY_2D<TT>& phi_real, const bool& copy_real_data);
	void FillGhostCellsFromDirichlet(const int& thread_id, ARRAY_2D<TT>& phi_real, const bool& copy_real_data);
	void FillGhostCellsFromFreeSlip(const int& thread_id, ARRAY_2D<TT>& phi_real, const bool& copy_real_data);

	void PrintKPlaneValues(const int& k, const bool& print_ghost);

	void Translate(const VT& deviation)
	{
		grid.Translate(deviation);
		grid_ghost.Translate(deviaiton);

		for(int i = 0; i < multithreading->num_threads; i++)
		{
			partial_grids[i].Translate(deviation);
			partial_grids_ghost[i].Translate(deviation);
		}

		// Translate SpeedUp variables
		x_min += deviation.x;
		y_min += deviation.y;
	}
};








		










		


		