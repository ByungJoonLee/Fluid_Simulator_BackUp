#pragma once

#include "COMMON_DEFINITIONS.h"
#include "FIELD_STRUCTURE_3D.h"
#include "LEVELSET_OBJECT.h"
#include "MULTITHREADING.h"
#include "MOVING_LEAST_SQUARES.h"
#include "TRIANGULAR_SURFACE.h"
#include "DYNAMIC_ARRAY.h"
#include "LOG.h"

class LEVELSET_3D : public LEVELSET_OBJECT
{
public: // Essential Data
	FIELD_STRUCTURE_3D<T>		signed_distance_field;
	FIELD_STRUCTURE_3D<T>		scalar_field_ghost;				// For interfacial smoothing

	GRID_STRUCTURE_3D&			grid;							
	ARRAY<GRID_STRUCTURE_3D>&	partial_grids;					// For the multithreading and related to the signed_distance_field
	ARRAY<GRID_STRUCTURE_3D>&	partial_grids_ghost;			// For the multithreading and related to the signed_distance_field
	int&						ghost_width;
	ARRAY_3D<T>					&phi, &arr;						// Related to the signed_distance_field

	MULTITHREADING*				multithreading;

	FIELD_STRUCTURE_3D<VT>		normal;
	FIELD_STRUCTURE_3D<T>		curvature;

	int							sweep_direction;
	
	T							max_abs_value;					// For debug
	T							abs_max_curvature;

public: // The way of calculating curvature vector
	bool						curvature_by_normal_vector;
	bool						curvature_by_levelset;

public: // Subfunctions for calculating curvature vector;
	FIELD_STRUCTURE_3D<T>		phi_x, phi_y, phi_z, phi_xx, phi_yy, phi_zz, phi_xy, phi_xz, phi_yz;

public: // Speedup Varibles
	ARRAY_3D<bool>				fixed;

public: // For drawing
	bool						is_periodic;
	int							periodic_num_x, periodic_num_y;

public: // Constructors and Destructor
	LEVELSET_3D(void)
		: grid(signed_distance_field.grid), partial_grids(signed_distance_field.partial_grids), partial_grids_ghost(signed_distance_field.partial_grids_ghost), ghost_width(signed_distance_field.ghost_width),
		 phi(signed_distance_field.array_for_this), arr(signed_distance_field.array_for_this), sweep_direction(0), max_abs_value(0), abs_max_curvature(0), is_periodic(false), periodic_num_x((int)1), periodic_num_y((int)1)
	{}

	~LEVELSET_3D(void)
	{}

public: // Operator Overloading
	inline T& operator()(const VI& ix) const
	{
		return phi(ix);
	}

	inline T& operator()(const int& i, const int& j, const int& k) const
	{
		return phi(i, j, k);
	}

	inline T operator()(const VT& position) const
	{
		return TriLinearInterpolation(position);
	}

public: // Initialization Functions
	void Initialize(MULTITHREADING* multithreading_input, const GRID_STRUCTURE_3D& grid_input, const int& ghost_cell_width_input)
	{
		assert(ghost_cell_width_input >= 2);		// When we calculate the curvature, we need to make ghost_width 2, at least

		multithreading = multithreading_input;
		signed_distance_field.Initialize(multithreading_input, grid_input, ghost_cell_width_input);
		scalar_field_ghost.Initialize(multithreading_input, grid_input, ghost_cell_width_input);

		fixed.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.k_start - ghost_width, 
						 grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, grid.k_res + 2*ghost_width, false);

		normal.Initialize(multithreading_input, grid_input, 1);
		curvature.Initialize(multithreading_input, grid_input, 0);

		// Initialize subfunctions
		phi_x.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_y.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_z.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_xx.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_yy.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_zz.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_xy.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_xz.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_yz.Initialize(multithreading, signed_distance_field.grid, 2);
		
		phi.AssignAllValues(grid.dx);

		sweep_direction = 0;

		max_abs_value = 0;

		curvature_by_normal_vector = false;
		curvature_by_levelset = false;
	}

	void Initialize(MULTITHREADING* multithreading_input, const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input, const int& ghost_width_input)
	{
		assert(ghost_width_input >= 2);				// When we calculate the curvature, we need to make ghost_width 2, at least

		multithreading = multithreading_input;
		signed_distance_field.Initialize(multithreading_input, i_res_input, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input, ghost_width);
		scalar_field_ghost.Initialize(multithreading_input, i_res_input, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input, ghost_width);

		fixed.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.k_start - ghost_width, 
						 grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, grid.k_res + 2*ghost_width, false);


		normal.Initialize(multithreading_input, i_res_input+2, j_res_input+2, k_res_input+2, i_start_input-1, j_start_input-1, k_start_input-1, x_min_input-grid.dx, y_min_input-grid.dy, z_min_input-grid.dz, x_max_input+grid.dx, y_max_input+grid.dy, z_max_input+grid.dz, 0);
		curvature.Initialize(multithreading_input, i_res_input, j_res_input, k_res_input, i_start_input, j_start_input, k_start_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input, 0);

		// Initialize subfunctions
		phi_x.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_y.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_z.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_xx.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_yy.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_zz.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_xy.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_xz.Initialize(multithreading, signed_distance_field.grid, 2);
		phi_yz.Initialize(multithreading, signed_distance_field.grid, 2);

		phi.AssignAllValues(grid.dx);

		sweep_direction = 0;

		max_abs_value = 0;

		curvature_by_normal_vector = false;
		curvature_by_levelset = false;
	}

public: // Member Functions
	inline VT CellCenter(const int& i, const int& j, const int& k) const
	{
		return grid.CellCenter(i, j, k);
	}

	void FillGhostCellsFrom(const int& thread_id, ARRAY_3D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(thread_id, phi_real, copy_real_data);
	}

	void FillGhostCellsFromPointer(const int& thread_id, ARRAY_3D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(thread_id, *phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivatesFrom(const int& thread_id, ARRAY_3D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(thread_id, phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivatesFromPointer(const int& thread_id, ARRAY_3D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(thread_id, *phi_real, copy_real_data);
	}

	void FillGhostCellsPeriodicInYDirection(const int& thread_id, ARRAY_3D<T>& phi_real, ARRAY_3D<T>& boundary_levelset, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsPeriodicInYDirection(thread_id, phi_real, boundary_levelset, copy_real_data);	
	}

	void FillGhostCellsFromThreaded(ARRAY_3D<T>* phi_real, const bool& copy_real_data)
	{
		multithreading->RunThreads(&LEVELSET_3D::FillGhostCellsFromPointer, this, phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivatesFromThreaded(ARRAY_3D<T>* phi_real, const bool& copy_real_data)
	{
		multithreading->RunThreads(&LEVELSET_3D::FillGhostCellsContinuousDerivatesFromPointer, this, phi_real, copy_real_data);
	}

	inline T TriLinearInterpolation(const VT& position) const
	{
		return signed_distance_field.TriLinearInterpolation(position);
	}

	void MinMaxIndexLevelSet(VI& min_ix, VI& max_ix)
	{
		min_ix = VI();
		max_ix = VI();
	
		T inner_range = grid.dx*(T)0.5;

		GRID_ITERATION_3D(grid)
		{
			if(phi(i, j, k) < inner_range)
			{
				min_ix.i = MIN(i, min_ix.i);
				min_ix.j = MIN(j, min_ix.j);
				min_ix.k = MIN(k, min_ix.k);

				max_ix.i = MAX(i, max_ix.i);
				max_ix.j = MAX(j, max_ix.j);
				max_ix.k = MAX(k, max_ix.k);
			}
		}
	}

	void MinMaxPositionLevelset(VT& min_pos, VT& max_pos)
	{
		VI min_ix, max_ix;
		MinMaxIndexLevelSet(min_ix, max_ix);

		min_pos = grid.CellCenter(min_ix);
		max_pos = grid.CellCenter(max_ix);

		min_pos -= VT(grid.dx, grid.dy, grid.dz);
		max_pos += VT(grid.dx, grid.dy, grid.dz);
	}

	void MinMaxCellTriangle(const VT& v0, const VT& v1, const VT& v2, VI& min, VI& max)
	{
		VT min_pos(MIN3(v0.x, v1.x, v2.x), MIN3(v0.y, v1.y, v2.y), MIN3(v0.z, v1.z, v2.z));
		VT max_pos(MAX3(v0.x, v1.x, v2.x), MAX3(v0.y, v1.y, v2.y), MAX3(v0.z, v1.z, v2.z));

		min_pos = min_pos - VT(grid.dx*(T)3, grid.dy*(T)3, grid.dz*(T)3);
		max_pos = max_pos + VT(grid.dx*(T)3, grid.dy*(T)3, grid.dz*(T)3);

		min = grid.ClampedCell(min_pos);
		max = grid.ClampedCell(max_pos);
	}

	void MinMaxCellTriangle(const TRIANGLE& triangle, VI& min, VI& max)
	{
		VT v0(triangle.vertices[0]->x);
		VT v1(triangle.vertices[1]->x);
		VT v2(triangle.vertices[2]->x);

		VT min_pos(MIN3(v0.x, v1.x, v2.x), MIN3(v0.y, v1.y, v2.y), MIN3(v0.z, v1.z, v2.z));
		VT max_pos(MAX3(v0.x, v1.x, v2.x), MAX3(v0.y, v1.y, v2.y), MAX3(v0.z, v1.z, v2.z));

		min_pos = min_pos - VT(grid.dx*(T)3, grid.dy*(T)3, grid.dz*(T)3);
		max_pos = max_pos + VT(grid.dx*(T)3, grid.dy*(T)3, grid.dz*(T)3);

		min = grid.ClampedCell(min_pos);
		max = grid.ClampedCell(max_pos);
	}

	void MinMaxCellTriangle(const GRID_STRUCTURE_3D& grid_input, const TRIANGLE& triangle, VI& min, VI& max)
	{
		VT v0(triangle.vertices[0]->x);
		VT v1(triangle.vertices[1]->x);
		VT v2(triangle.vertices[2]->x);

		VT min_pos(MIN3(v0.x, v1.x, v2.x), MIN3(v0.y, v1.y, v2.y), MIN3(v0.z, v1.z, v2.z));
		VT max_pos(MAX3(v0.x, v1.x, v2.x), MAX3(v0.y, v1.y, v2.y), MAX3(v0.z, v1.z, v2.z));

		min_pos = min_pos - VT(grid_input.dx*(T)3, grid_input.dy*(T)3, grid_input.dz*(T)3);
		max_pos = max_pos + VT(grid_input.dx*(T)3, grid_input.dy*(T)3, grid_input.dz*(T)3);

		min = grid_input.ClampedCell(min_pos);
		max = grid_input.ClampedCell(max_pos);
	}

	void MinMaxCellEdge(const EDGE& edge, VI& min, VI& max)
	{
		VT vertex_position = VT(edge.vertices[0]->x[0], edge.vertices[0]->x[1], edge.vertices[0]->x[2]);

		VI min_index = grid.ClampedCell(vertex_position);
		VI max_index = min_index;

		const int edge_num = (int)edge.vertices.size();
		for(int i = 1; i < edge_num; i++)
		{
			vertex_position = VT(edge.vertices[i]->x[0], edge.vertices[i]->x[1], edge.vertices[i]->x[2]);

			VI vertex_index = grid.ClampedCell(vertex_position);

			if(min_index.i > vertex_index.i) min_index.i = vertex_index.i;
			if(min_index.j > vertex_index.j) min_index.j = vertex_index.j;
			if(min_index.k > vertex_index.k) min_index.k = vertex_index.k;

			if(max_index.i < vertex_index.i) max_index.i = vertex_index.i;
			if(max_index.j < vertex_index.j) max_index.j = vertex_index.j;
			if(max_index.k < vertex_index.k) max_index.k = vertex_index.k;
		}

		min = min_index;
		max = max_index;
	}

	void Translate(const VT& deviation)
	{
		signed_distance_field.Translate(deviation);
		normal.Translate(deviation);
		curvature.Translate(deviation);
	}

	void SweepFill(const VI& start_idx, const VI& end_idx, const T& target_value, const T& assign_value)
	{
		int i = start_idx.i;
		int j = start_idx.j;
		int k = start_idx.k;

		if(start_idx.k < end_idx.k)
		{
			for(; k <= end_idx.k; k++)
			{
				T& phi_ = phi(i, j, k);
				if(phi_ == target_value) phi_ = assign_value;
				else break;
			}
		}
		else
		{
			for(; k >= end_idx.k; k--)
			{
				T& phi_ = phi(i, j, k);
				if(phi_ == target_value) phi_ = assign_value;
				else break;
			}
		}
	}

	void FloodFillNonRecursive(const int& i, const int& j, const int& k, const T& old_value, const T& new_value)
	{
		if(!grid.Inside(i, j, k)) return;

		DYNAMIC_ARRAY<VI> link_index_arr;
		link_index_arr.Initialize(phi.ijk_res, phi.ijk_res);
		link_index_arr.Push(VI(i, j, k));
		
		while(link_index_arr.num_of_elements)
		{
			VI idx;
			link_index_arr.Pop(idx);

			phi(idx.i, idx.j, idx.k) = new_value;

			if(grid.Inside(idx.i + 1, idx.j, idx.k) && phi(idx.i + 1, idx.j, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i + 1, idx.j, idx.k));
			if(grid.Inside(idx.i - 1, idx.j, idx.k) && phi(idx.i - 1, idx.j, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i - 1, idx.j, idx.k));
			if(grid.Inside(idx.i, idx.j + 1, idx.k) && phi(idx.i, idx.j + 1, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j + 1, idx.k));
			if(grid.Inside(idx.i, idx.j - 1, idx.k) && phi(idx.i, idx.j - 1, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j - 1, idx.k));
			if(grid.Inside(idx.i, idx.j, idx.k + 1) && phi(idx.i, idx.j, idx.k + 1) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j, idx.k + 1));
			if(grid.Inside(idx.i, idx.j, idx.k - 1) && phi(idx.i, idx.j, idx.k - 1) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j, idx.k - 1));
		}
	}

	void FloodFillNonRecursive(const int& i, const int& j, const int& k)
	{
		if(!grid.Inside(i, j, k)) return;

		T old_value = phi(i, j, k);
		DYNAMIC_ARRAY<VI> link_index_arr;
		link_index_arr.Initialize(phi.ijk_res, phi.ijk_res);

		while(link_index_arr.num_of_elements)
		{
			VI idx;
			link_index_arr.Pop(idx);

			phi(idx.i, idx.j, idx.k) = old_value*(T)-1;

			if(grid.Inside(idx.i + 1, idx.j, idx.k) && phi(idx.i + 1, idx.j, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i + 1, idx.j, idx.k));
			if(grid.Inside(idx.i - 1, idx.j, idx.k) && phi(idx.i - 1, idx.j, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i - 1, idx.j, idx.k));
			if(grid.Inside(idx.i, idx.j + 1, idx.k) && phi(idx.i, idx.j + 1, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j + 1, idx.k));
			if(grid.Inside(idx.i, idx.j - 1, idx.k) && phi(idx.i, idx.j - 1, idx.k) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j - 1, idx.k));
			if(grid.Inside(idx.i, idx.j, idx.k + 1) && phi(idx.i, idx.j, idx.k + 1) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j, idx.k + 1));
			if(grid.Inside(idx.i, idx.j, idx.k - 1) && phi(idx.i, idx.j, idx.k - 1) == (old_value))
				link_index_arr.Push(VI(idx.i, idx.j, idx.k - 1));
		}
	}

	inline bool InsideOBB(const VT& position) const
	{
		return grid.Inside(position);
	}

	inline const T SignedDistance(const VT& position) const
	{
		return TriLinearInterpolation(position);
	}

	void AssignAllValuesLevelsetThreaded(const T& value)
	{
		multithreading->RunThreads(&LEVELSET_3D::AssignAllValuesLevelset, this, value);
	}

	void AssignAllValuesLevelset(const int& thread_id, const T& value)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			GRID_STRUCTURE_3D& grid_input = partial_grids[thread_id];
			phi(i, j, k) = value;
		}
		END_GRID_ITERATION_3D;
	}

	void AssignMinPositiveValuesLevelset(const int& thread_id)
	{
		GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = grid.dx*(T)3;
		}
	}

	void AssignMinPositiveValuesLevelsetThreaded()
	{
		multithreading->RunThreads(&LEVELSET_3D::AssignMinPositiveValuesLevelset, this);
	}

	void AssignMinNegativeValuesLevelset(const int& thread_id)
	{
		GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = -grid.dx*(T)3;
		}
	}

	void AssignMinNegativeValuesLevelsetThreaded()
	{
		multithreading->RunThreads(&LEVELSET_3D::AssignMinNegativeValuesLevelset, this);
	}

	void AssignTriangleLevelsetBarycentric(const TRIANGLE& triangle)
	{
		VI min, max;
		MinMaxCellTriangle(triangle, min, max);

		VERTEX* v0 = triangle.vertices[0];
		VERTEX* v1 = triangle.vertices[1];
		VERTEX* v2 = triangle.vertices[2];

		const VT p0 = VT(v0->x);
		const VT p1 = VT(v1->x);
		const VT p2 = VT(v2->x);

		const VT n0 = VT(v0->n);
		const VT n1 = VT(v1->n);
		const VT n2 = VT(v2->n);

		int i, j, k;
		LOOPS_3D(i, j, k, min.i, min.j, min.k, max.i, max.j, max.k)
		{
			T& cell_phi = phi(i, j, k);
			VT cell_point = grid.CellCenter(i, j, k);
			
			const T phi0 = DotProduct(n0, cell_point - p0);
			const T phi1 = DotProduct(n1, cell_point - p1);
			const T phi2 = DotProduct(n2, cell_point - p2);

			VT w = triangle.BaryCentricCoordinates(cell_point, p0, p1, p2);
			
			w.x = MAX(w.x, 0);
			w.y = MAX(w.y, 0);
			w.z = MAX(w.z, 0);

			cell_phi = MIN_ABS(phi0*w.x + phi1*w.y + phi2*w.z, cell_phi);
		}
	}

	void AssignTriangularSurfaceLevelSet(TRIANGULAR_SURFACE& triangular_surface)
	{
		list<TRIANGLE*>::iterator itr_triangle = triangular_surface.triangles.begin();
		for(; itr_triangle != triangular_surface.triangles.end(); itr_triangle++)
		{
			TRIANGLE& triangle = **itr_triangle;

			VI min, max;
			MinMaxCellTriangle(triangle, min, max);

			int i, j, k;
			LOOPS_3D(i, j, k, min.i, min.j, min.k, max.i, max.j, max.k)
			{
				VT cell_point = grid.CellCenter(i, j, k);
				T signed_distance = triangle.SignedDistance(cell_point);

				phi(i, j, k) = MIN_ABS(phi(i, j, k), signed_distance);
			}
		}
	}

	void AssignMLSLevelset(MOVING_LEAST_SQUARES* mls)
	{
		GRID_ITERATION_3D(grid)
		{
			phi(i, j, k) = mls->GetScalar(grid.CellCenter(i, j, k));
		}
	}

	void AssignMLSLevelset(const int& thread_id, MOVING_LEAST_SQUARES* mls)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = mls->GetScalar(grid.CellCenter(i, j, k));
		}
		END_GRID_ITERATION_3D;
	}

	void AssignObjectLevelset(const LEVELSET_OBJECT* obj, const T& value)
	{
		GRID_ITERATION_3D(grid)
		{
			phi(i, j, k) = (*obj)(grid.CellCenter(i, j, k));
		}
	}

	void UnionObjectLevelset(const int& thread_id, const LEVELSET_OBJECT* obj)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = MIN(phi(i, j, k), (*obj)(grid.CellCenter(i, j, k)));
		}
		END_GRID_ITERATION_3D;
	}

	void UnionObjectLevelset(const LEVELSET_OBJECT* obj)
	{
		GRID_ITERATION_3D(grid)
		{
			phi(i, j, k) = MIN(phi(i, j, k), (*obj)(grid.CellCenter(i, j, k)));
		}
	}

	void InterfacialSmoothing(const int& thread_id)
	{
		cout << "InterfacialSmoothing()" << endl;

		VT position, closest_position;
		T  closest_distance;
		VT unit_normal;
		
		// Ghost filled levelset copy
		scalar_field_ghost.FillGhostCellsFrom(thread_id, phi, true);
		ARRAY_3D<T>& phi_copy = scalar_field_ghost.array_for_this;

		bool is_interfacial, is_phi_positive;
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			// Check interfacial cell
			is_interfacial = false;
			
			if(phi_copy(i, j, k) > (T)0)
			{
				is_phi_positive = true;
				if(phi_copy(i+1, j, k) <= (T)0)			is_interfacial = true;
				else if(phi_copy(i-1, j, k) <= (T)0)	is_interfacial = true;
				else if(phi_copy(i, j+1, k) <= (T)0)	is_interfacial = true;
				else if(phi_copy(i, j-1, k) <= (T)0)	is_interfacial = true;
				else if(phi_copy(i, j, k+1) <= (T)0)	is_interfacial = true;
				else if(phi_copy(i, j, k-1) <= (T)0)	is_interfacial = true;
			}
			else
			{
				is_phi_positive = false;
				if(phi_copy(i+1, j, k) > (T)0)			is_interfacial = true;
				else if(phi_copy(i-1, j, k) > (T)0)		is_interfacial = true;
				else if(phi_copy(i, j+1, k) > (T)0)		is_interfacial = true;
				else if(phi_copy(i, j-1, k) > (T)0)		is_interfacial = true;
				else if(phi_copy(i, j, k+1) > (T)0)		is_interfacial = true;
				else if(phi_copy(i, j, k-1) > (T)0)		is_interfacial = true;
			}

			if(is_interfacial == false) continue;

			// For the interfacial cell
			// Position
			position = grid.CellCenter(i, j, k);
			// Closest Distance
			closest_distance = scalar_field_ghost(i, j, k);
			// Unit Normal
			unit_normal.x = (phi_copy(i+1, j, k) - phi_copy(i-1, j, k))*signed_distance_field.one_over_2dx;
			unit_normal.y = (phi_copy(i, j+1, k) - phi_copy(i, j-1, k))*signed_distance_field.one_over_2dy;
			unit_normal.z = (phi_copy(i, j, k+1) - phi_copy(i, j, k-1))*signed_distance_field.one_over_2dz;
			unit_normal.Normalize();

			closest_position = position - unit_normal*closest_distance;

			if(is_phi_positive) phi(i, j, k) =  FastDistance(position, closest_position);
			else				phi(i, j, k) = -FastDistance(position, closest_position);
		}
		END_GRID_ITERATION_3D;
	}

	void ReinitializeByFastSweeping(const int& thread_id)
	{
		LOG::Begin(thread_id, "Reinitialization by Fast Sweeping");

		FillGhostCellsFrom(thread_id, phi, false);
		FastSweepingMethod(thread_id);
		FillGhostCellsFrom(thread_id, phi, false);
		Update(thread_id);

		LOG::End(thread_id);
	}

	void ReintializeByFastSweepingThreaded()
	{
		multithreading->RunThreads(&LEVELSET_3D::ReinitializeByFastSweeping, this);
	}

	void CopyAllValuesFrom(const int& thread_id, const LEVELSET_3D& from_levelset)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = from_levelset.phi(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void CopyAllValuesGhostFrom(const int& thread_id, const LEVELSET_3D& from_levelset)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids_ghost[thread_id])
		{
			phi(i, j, k) = from_levelset.phi(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void UnionWith(const int& thread_id, const LEVELSET_3D& from_levelset)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = MIN(from_levelset.phi(i, j, k), phi(i, j, k));
		}
		END_GRID_ITERATION_3D;
	}

	void AssignTriangle(const VT& v0, const VT& v1, const VT& v2)
	{
		VI min, max;
		MinMaxCellTriangle(v0, v1, v2, min, max);

		int i, j, k;
		LOOPS_3D(i, j, k, min.i, min.j, min.k, max.i, max.j, max.k)
		{
			VT cell_point = grid.CellCenter(i, j, k);

			T  distance = TRIANGLE::SignedDistanceFromTriangle(cell_point, v0, v1, v2);
			T& phi_ = phi(i, j, k);

			phi_ = MIN_ABS(distance, phi_);
		}
	}

	void AssignTriangularMeshLevelset(const int& thread_id, const DYNAMIC_ARRAY<VT*>* vertices, const DYNAMIC_ARRAY<VI>* indices, DYNAMIC_ARRAY<VI>* sweep_list)
	{
		T target_value = (T) grid.dx*(T)3;
		T assign_value = (T)-grid.dx*(T)3;

		AssignAllValuesLevelset(thread_id, target_value);

		PREPARE_FOR_1D_ITERATION(indices->num_of_elements);
		BEGIN_1D_ITERATION
		{
			const DYNAMIC_ARRAY<VT*>& v = *(vertices);

			VI& face = indices->values[p];
			VT& v0 = *(v.values[face.i]);
			VT& v1 = *(v.values[face.j]);
			VT& v2 = *(v.values[face.k]);

			BEGIN_SCOPED_LOCK
				AssignTriangle(v0, v1, v2);
			END_SCOPED_LOCK
		}
		END_1D_ITERATION;

		// Fine cell starting sweep cell
		BEGIN_HEAD_THREAD_WORK
			sweep_list->Emptify();
		END_HEAD_THREAD_WORK

		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			if(phi(i, j, k) < (T)0)
			{
				int k_1 = grid.ClampedIndexK(k+1);
				if(phi(i, j, k_1) == target_value)
				{
					for(int end_k = k_1; end_k <= grid.k_end; end_k++)
					{
						T& phi_ = phi(i, j, end_k);
						if(phi_ != target_value)
						{
							BEGIN_SCOPED_LOCK
							{
								sweep_list->Push(VI(i, j, k_1));
								sweep_list->Push(VI(i, j, end_k));
							}
							END_SCOPED_LOCK
							break;
						}
					}
				}
			}
		}
		END_GRID_ITERATION_3D;

		PREPARE_FOR_1D_ITERATION(sweep_list->num_of_elements/2);
		BEGIN_1D_ITERATION
		{
			int n = p*2;
			VI sta_idx = sweep_list->values[n+0];
			VI end_idx = sweep_list->values[n+1];

			SweepFill(sta_idx, end_idx, target_value, assign_value);
		}
		END_1D_ITERATION;
	}

	void AssignTriangularMeshLevelsetThreaded(const DYNAMIC_ARRAY<VT*>* vertices, const DYNAMIC_ARRAY<VI>* indices, DYNAMIC_ARRAY<VI>* sweep_list)
	{
		multithreading->RunThreads(&LEVELSET_3D::AssignTriangularMeshLevelset, this, vertices, indices, sweep_list);
	}

	void FastSweepingMethod(const int& thread_id)
	{
		const T nb_width = (T)20*grid.dx;
		
		signed_distance_field.AssignAllValues<bool>(thread_id, fixed, false);

		GRID_ITERATION_3D(partial_grids[thread_id])
		{
			bool& fix_center(fixed(i, j, k));

			if(fix_center == true) continue;

			T& phi_center(phi(i, j, k));

			// Fix values of the cells too far from interface
			if(ABS(phi_center) > nb_width)
			{
				fix_center = true;
				continue;
			}

			// Fix (do not update) interfacial cells
			if(phi_center > (T)0)
			{
				if(phi(i+1, j, k) <= (T)0) {fix_center = true; fixed(i+1, j, k) = true;}
				else if(phi(i-1, j, k) <= (T)0) {fix_center = true; fixed(i-1, j, k) = true;}
				else if(phi(i, j+1, k) <= (T)0) {fix_center = true; fixed(i, j+1, k) = true;}
				else if(phi(i, j-1, k) <= (T)0) {fix_center = true; fixed(i, j-1, k) = true;}
				else if(phi(i, j, k+1) <= (T)0) {fix_center = true; fixed(i, j, k+1) = true;}
				else if(phi(i, j, k-1) <= (T)0) {fix_center = true; fixed(i, j, k-1) = true;}
			}
			else
			{
				if(phi(i+1, j, k) > (T)0) {fix_center = true; fixed(i+1, j, k) = true;}
				else if(phi(i-1, j, k) > (T)0) {fix_center = true; fixed(i-1, j, k) = true;}
				else if(phi(i, j+1, k) > (T)0) {fix_center = true; fixed(i, j+1, k) = true;}
				else if(phi(i, j-1, k) > (T)0) {fix_center = true; fixed(i, j-1, k) = true;}
				else if(phi(i, j, k+1) > (T)0) {fix_center = true; fixed(i, j, k+1) = true;}
				else if(phi(i, j, k-1) > (T)0) {fix_center = true; fixed(i, j, k-1) = true;}
			}

			// Assign large magnitude values to non-interfacial cells
			if(fix_center == false)
			{
				if(phi_center > (T)0) phi_center = nb_width;
				else phi_center = -nb_width;
			}
		}

		multithreading->Sync(thread_id);

		FillGhostCellsFrom(thread_id, phi, false);
		
		INIT_GRIDRANGE_3D(partial_grids[thread_id], i_start, j_start, k_start, i_end, j_end, k_end);

		switch(sweep_direction % 4)
		{
		case 0:
			Sweep(i_start, j_start, k_start, i_end, j_end, k_end); 
			multithreading->Sync(thread_id);
			Sweep(i_end, j_end, k_end, i_start, j_start, k_start);
			multithreading->Sync(thread_id);
			break;

		case 1:
			Sweep(i_end, j_start, k_end, i_start, j_end, k_start);
			multithreading->Sync(thread_id);
			Sweep(i_start, j_end, k_start, i_end, j_start, k_end);
			multithreading->Sync(thread_id);
			break;

		case 2:
			Sweep(i_start, j_end, k_end, i_end, j_start, k_start);
			multithreading->Sync(thread_id);
			Sweep(i_end, j_start, k_start, i_start, j_end, k_end);
			multithreading->Sync(thread_id);
			break;

		case 3:
			Sweep(i_end, j_end, k_start, i_start, j_start, k_end);
			multithreading->Sync(thread_id);
			Sweep(i_start, j_start, k_end, i_end, j_end, k_start);
			multithreading->Sync(thread_id);
			break;
		}

		HEAD_THREAD_WORK(sweep_direction++);
	}

	void FastSweepingMethodThreaded()
	{
		multithreading->RunThreads(&LEVELSET_3D::FastSweepingMethod, this);
	}

	void FullFastSweepingMethodThreaded()
	{
		for(int d = 0; d < 4; d++)
		{
			sweep_direction = d;
			multithreading->RunThreads(&LEVELSET_3D::FastSweepingMethod, this);
		}
	}

	void Sweep(const int& i_start, const int& j_start, const int& k_start, const int& i_end, const int& j_end, const int& k_end)
	{
		const T one_over_three((T)1/(T)3), one_over_two((T)1/(T)2);
		const T h(grid.dx), hh(h*h), hh2(hh*(T)2), hh3(hh*(T)3);

		T a, b, c;
		T a1, a2, a3;
		T update;

		// Traverse orders
		int ip(1), jp(1), kp(1);		// Increasing order
		if(i_start > i_end) ip = -1;	// Decreasing order
		if(j_start > j_end) jp = -1;	// Decreasing order
		if(k_start > k_end) kp = -1;	// Decreasing order

		int i, j, k;
		k = k_start - kp;
		while(true)
		{
			k += kp;
			j = j_start - jp;
			while(true)
			{
				j += jp;
				i = i_start - ip;
				while(true)
				{
					i += ip;

					T& phi_center(phi(i, j, k));
					if(fixed(i, j, k) == false)
					{
						if(phi_center > (T)0)
						{
							// See definitions of u^h_{x_min}, u^h_{y_min} in page 605
							a = MIN(phi(i-1, j, k), phi(i+1, j, k));
							b = MIN(phi(i, j-1, k), phi(i, j+1, k));
							c = MIN(phi(i, j, k-1), phi(i, j, k+1));

							// Follows the n dimension implementation explained in 606
							// See "First we order the a_k's in ~"
							INCREASING_SORT3(a, b, c, a1, a2, a3);
							
							update = a1 + h;
							if(update > a2)
							{
								update = one_over_two*(a1 + a2 + sqrt(hh2 - POW2(a1 - a2)));

								if(update > a3)
									update = one_over_three*(a1 + a2 + a3 + sqrt(POW2(a1 + a2 + a3) - (T)3*(a1*a1 + a2*a2 + a3*a3 - hh)));
							}

							phi_center = MIN(update, phi_center);
						}
						else
						{
							// See Equation (2.8) in page 607 for the implementation of negative parts
							a = MAX(phi(i-1, j, k), phi(i+1, j, k));
							b = MAX(phi(i, j-1, k), phi(i, j+1, k));
							c = MAX(phi(i, j, k-1), phi(i, j, k+1));
							
							INCREASING_SORT3(a, b, c, a1, a2, a3);
							
							update = a3 - h;
							if(update < a2)
							{
								update = one_over_two*(a3 + a2 - sqrt(hh2 - POW2(a3 - a2)));

								if(update < a1)
									update = one_over_three*(a1 + a2 + a3 - sqrt(POW2(a1 + a2 + a3) - (T)3*(a1*a1 + a2*a2 + a3*a3 - hh)));
							}

							phi_center = MAX(update, phi_center);
						}
					}
					if(i == i_end) break;
				}
				if(j == j_end) break;
			}
			if(k == k_end) break;
		}
	}
	
	void AssignObjectLevelset(const int& thread_id, const LEVELSET_OBJECT* obj)
	{
		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			phi(i, j, k) = (*obj)(CellCenter(i, j, k));
		}
		END_GRID_ITERATION_3D;
	}

	void ComputeNormals(const int& thread_id)
	{
		signed_distance_field.FillGhostCellsFrom(thread_id, phi, false);

		BEGIN_GRID_ITERATION_3D(normal.partial_grids[thread_id])
		{
			VT& nor(normal(i, j, k));

			T nor_x = (phi(i+1, j, k) - phi(i-1, j, k))*signed_distance_field.one_over_2dx, nor_y = (phi(i, j+1, k) - phi(i, j-1, k))*signed_distance_field.one_over_2dy, nor_z = (phi(i, j, k+1) - phi(i, j, k-1))*signed_distance_field.one_over_2dz;
			T mag_of_normal = sqrt(POW2(nor_x) + POW2(nor_y) + POW2(nor_z));

			if (mag_of_normal != 0)
			{
				nor.x = nor_x;
				nor.y = nor_y;
				nor.z = nor_z;
			}
			else
			{
				nor.x = (phi(i + 1, j,  k) - phi(i, j, k))*signed_distance_field.one_over_dx;
				nor.y = (phi(i, j + 1,  k) - phi(i, j, k))*signed_distance_field.one_over_dy;
				nor.z = (phi(i, j,  k + 1) - phi(i, j, k))*signed_distance_field.one_over_dz;

			}
			nor.Normalize();
		}
		END_GRID_ITERATION_3D;
	}

	void ComputeNormalsThreaded()
	{
		multithreading->RunThreads(&LEVELSET_3D::ComputeNormals, this);
	}

	void ComputeCurvatures(const int& thread_id)
	{
		T tolerance = (T)1/MIN3(grid.dx, grid.dy, grid.dz);
		
		BEGIN_HEAD_THREAD_WORK
		{
			abs_max_curvature = (T)0;
		}
		END_HEAD_THREAD_WORK;

		// Curvature calculated by Normal
		if (curvature_by_normal_vector)
		{
			normal.FillGhostCellsFrom(thread_id, normal.array_for_this, false);
			
			BEGIN_GRID_ITERATION_3D(curvature.partial_grids[thread_id])
			{
				T curv((normal(i+1, j, k).x - normal(i-1, j, k).x)*signed_distance_field.one_over_2dx);
				curv += (normal(i, j+1, k).y - normal(i, j-1, k).y)*signed_distance_field.one_over_2dy;
				curv += (normal(i, j, k+1).z - normal(i, j, k-1).z)*signed_distance_field.one_over_2dz;

				if (abs(curv) <= tolerance)
				{
					curvature(i, j, k) = curv;
				}
				else
				{
					curvature(i, j, k) = -tolerance;
				}
			
				abs_max_curvature = MAX(abs(curv), abs_max_curvature);
			}
			END_GRID_ITERATION_MAX_3D(abs_max_curvature);
		}
		
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(thread_id, signed_distance_field.array_for_this, true);

		// Curvature calculated by Levelset
		if (curvature_by_levelset)
		{
			BEGIN_GRID_ITERATION_3D(signed_distance_field.partial_grids[thread_id])
			{
				phi_x(i, j, k) = (phi(i + 1, j, k) - phi(i - 1, j, k))*signed_distance_field.one_over_2dx;
				phi_y(i, j, k) = (phi(i, j + 1, k) - phi(i, j - 1, k))*signed_distance_field.one_over_2dy;
				phi_z(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k - 1))*signed_distance_field.one_over_2dz;
				phi_xx(i, j, k) = (phi(i + 1, j, k) - 2*phi(i, j, k) + phi(i - 1, j, k))*signed_distance_field.one_over_dx2;
				phi_yy(i, j, k) = (phi(i, j + 1, k) - 2*phi(i, j, k) + phi(i, j - 1, k))*signed_distance_field.one_over_dy2;
				phi_zz(i, j, k) = (phi(i, j, k + 1) - 2*phi(i, j, k) + phi(i, j, k - 1))*signed_distance_field.one_over_dz2;
				phi_xy(i, j, k) = (phi(i + 1, j + 1, k) - phi(i + 1, j - 1, k) - phi(i - 1, j + 1, k) + phi(i - 1, j - 1, k))*signed_distance_field.one_over_2dx*signed_distance_field.one_over_2dy;
				phi_xz(i, j, k) = (phi(i + 1, j, k + 1) - phi(i + 1, j, k - 1) - phi(i - 1, j, k + 1) + phi(i - 1, j, k - 1))*signed_distance_field.one_over_2dx*signed_distance_field.one_over_2dz;
				phi_yz(i, j, k) = (phi(i, j + 1, k + 1) - phi(i, j + 1, k - 1) - phi(i, j - 1, k + 1) + phi(i, j - 1, k - 1))*signed_distance_field.one_over_2dy*signed_distance_field.one_over_2dz;
			}
			END_GRID_ITERATION_3D;

			T px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz;

			BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
			{
				px = phi_x(i, j, k), py = phi_y(i, j, k), pz = phi_z(i, j, k);
				pxx = phi_xx(i, j, k), pyy = phi_yy(i, j, k), pzz = phi_zz(i, j, k);
				pxy = phi_xy(i, j, k), pxz = phi_xz(i, j, k), pyz = phi_yz(i, j, k);

				T magnitude = sqrt(POW2(px) + POW2(py) + POW2(pz));
				T deno = POW3(magnitude);

				T curv(0);

				if (deno != 0)
				{
					T one_over_deno = (T)1/deno;
					curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx + POW2(px)*pzz - (T)2*px*pz*pxz + POW2(pz)*pxx + POW2(py)*pzz - (T)py*pz*pyz + POW2(pz)*pyy)*one_over_deno;
				}
				else
				{
					px = (phi(i + 1, j, k) - phi(i, j, k))*grid.one_over_dx, py = (phi(i, j + 1, k) - phi(i, j, k))*grid.one_over_dy, pz = (phi(i, j, k + 1) - phi(i, j, k))*grid.dz;
					magnitude = sqrt(POW2(px) + POW2(py) + POW2(pz));
					deno = POW3(magnitude);
					if (deno == 0)
					{
						cout << "Denominator cannot be zero!!" << endl;
						exit(0);
					}
					T one_over_deno = (T)1/deno;
					curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx + POW2(px)*pzz - (T)2*px*pz*pxz + POW2(pz)*pxx + POW2(py)*pzz - (T)py*pz*pyz + POW2(pz)*pyy)*one_over_deno;
				}
				
				/*if (abs(curv) <= tolerance)
				{
					curvature(i, j, k) = curv;
				}
				else if (curv > tolerance)
				{
					curvature(i, j, k) = tolerance;
				}
				else if (curv < -tolerance)
				{
					curvature(i, j, k) = -tolerance;
				}*/
				if (curv < -tolerance)
				{
					curvature(i, j, k) = -tolerance;
				}
				else
				{
					curvature(i, j, k) = curv;
				}
				
				if (abs(signed_distance_field(i, j, k)) < 1e-1)
				{
					abs_max_curvature = MAX(abs(curvature(i, j, k)), abs_max_curvature);
				}
				else
				{
					continue;
				}
				
			}
			END_GRID_ITERATION_MAX_3D(abs_max_curvature);
		}
	}

	void ComputeCurvaturesThreaded()
	{
		multithreading->RunThreads(&LEVELSET_3D::ComputeCurvatures, this);
	}

	void CurvatureFlow(const int& thread_id, const T& dt)
	{
		BEGIN_GRID_ITERATION_3D(curvature.partial_grids[thread_id])
		{
			phi(i, j, k) += curvature(i, j, k)*dt;
		}
		END_GRID_ITERATION_3D;
	}

	inline const VT Normal(const VT& position) const
	{
		// Do not fill ghost cells here to avoid duplicated operations
		return normal.TriLinearInterpolation(position);
	}

	inline void UnitNormal(const VT& position, VT& normal_output) const
	{
		normal_output = normal.TriLinearInterpolation(position);
		normal_output.Normalize();
	}

	inline const VT UnitNormal(const VT& position) const
	{
		return (normal.TriLinearInterpolation(position)).Normalized();
	}

	inline T Curvature(const VT& position) const
	{
		return curvature.TriLinearInterpolation(position);
	}

	void Update(const int& thread_id)
	{
		ComputeNormals(thread_id);
		ComputeCurvatures(thread_id);
	}

	void UpdateThreaded()
	{
		ComputeNormalsThreaded();
		ComputeCurvaturesThreaded();
	}

	void Write(const char* filename)
	{
		ofstream ost(filename, ios::binary);
		if(!ost.is_open()) cout << "LEVELSET_3D::Failed to write" << filename << endl;

		// To record grid data
		ost.write((char*)&grid.i_start, sizeof(grid.i_start));
		ost.write((char*)&grid.j_start, sizeof(grid.j_start));
		ost.write((char*)&grid.k_start, sizeof(grid.k_start));
		ost.write((char*)&grid.i_res, sizeof(grid.i_res));
		ost.write((char*)&grid.j_res, sizeof(grid.j_res));
		ost.write((char*)&grid.k_res, sizeof(grid.k_res));
		ost.write((char*)&grid.x_min, sizeof(grid.x_min));
		ost.write((char*)&grid.y_min, sizeof(grid.y_min));
		ost.write((char*)&grid.z_min, sizeof(grid.z_min));
		ost.write((char*)&grid.x_max, sizeof(grid.x_max));
		ost.write((char*)&grid.y_max, sizeof(grid.y_max));
		ost.write((char*)&grid.z_max, sizeof(grid.z_max));

		// To record width of ghost cell
		ost.write((char*)&ghost_width, sizeof(ghost_width));

		// To record phi values
		GRID_ITERATION_3D(grid)
			ost.write((char*)&phi(i,j,k), sizeof(T));
		
		ost.close();
	}

	void Read(const char* filename, MULTITHREADING* multithreading_input)
	{
		ifstream ist(filename, ios::binary);
		
		if(!ist.is_open()) cout << "LEVELSET_3D::Failed to read" << filename << endl;

		int i_start, j_start, k_start, i_res, j_res, k_res, ghost_width;
		T x_min, y_min, z_min, x_max, y_max, z_max, phi_;

		// To load grid data
		ist.read((char*)&i_start, sizeof(i_start));
		ist.read((char*)&j_start, sizeof(j_start));
		ist.read((char*)&k_start, sizeof(k_start));
		ist.read((char*)&i_res, sizeof(i_res));
		ist.read((char*)&j_res, sizeof(j_res));
		ist.read((char*)&k_res, sizeof(k_res));
		ist.read((char*)&x_min, sizeof(x_min));
		ist.read((char*)&y_min, sizeof(y_min));
		ist.read((char*)&z_min, sizeof(z_min));
		ist.read((char*)&x_max, sizeof(x_max));
		ist.read((char*)&y_max, sizeof(y_max));
		ist.read((char*)&z_max, sizeof(z_max));

		// To load width of ghost cell
		ist.read((char*)&ghost_width, sizeof(ghost_width));

		GRID_STRUCTURE_3D grid_input(i_res, j_res, k_res, i_start, j_start, k_start, x_min, y_min, z_min, x_max, y_max, z_max);

		// To initialize levelset
		Initialize(multithreading_input, grid_input, ghost_width);

		// To load phi values
		GRID_ITERATION_3D(grid_input)
		{
			ist.read((char*)&phi_, sizeof(T));
			phi(i, j, k) = phi_;
		}

		ist.close();
	}

	inline const VT ClosestPoint(const VT& position) const
	{
		const T closest_distance(signed_distance_field.TriLinearInterpolation(position));
		const VT unit_normal(normal.TriLinearInterpolation(position).Normalized());
		return position - unit_normal*closest_distance;
	}

	inline void ClosestPoint(VT& position, const int& itr, const T& phi_target) const
	{
		for(int i = 0 ; i < itr; i++)
		{
			const T phi_old = signed_distance_field.TriLinearInterpolation(position) - phi_target;
			const VT unit_normal = normal.TriLinearInterpolation(position).Normalized();
			const VT new_position = position - unit_normal*phi_old;
			const T phi_new = signed_distance_field.TriLinearInterpolation(position) - phi_target;

			if(ABS(phi_new) < ABS(phi_old)) position = new_position;
			else return;
		}
	}

	const T EstimateVolume(const int& i, const int& j, const int& k, const int& num_samples)
	{
		const VT cell_center = CellCenter(i, j, k);
		const T dx = grid.dx/(T)(num_samples + 1);

		int count = 0;
		for(T sx = cell_center.x - grid.dx*(T)0.5 + (T)0.5*dx; sx <= cell_center.x + (T)0.5*grid.dx; sx += dx)
		{
			for(T sy = cell_center.y - grid.dy*(T)0.5 + (T)0.5*dx; sy <= cell_center.y + (T)0.5*grid.dy; sy += dx)
			{
				for(T sz = cell_center.z - grid.dz*(T)0.5 + (T)0.5*dx; sz <= cell_center.z + (T)0.5*grid.dz; sz += dx)
				{
					if(signed_distance_field.TriLinearInterpolation(VT(sx, sy, sz)) <= (T)0) count++;
				}
			}
		}

		return (T)count/(T)POW3(num_samples);
	}

	void GetMaxValue(const int& thread_id)
	{
		max_abs_value = 0;

		BEGIN_GRID_ITERATION_3D(partial_grids[thread_id])
		{
			if (max_abs_value <= abs(signed_distance_field(i, j, k)))
			{
				max_abs_value = abs(signed_distance_field(i, j, k));
			}
		}
		END_GRID_ITERATION_3D;
	}

};

void Sampling(MULTITHREADING* multithreading, const int& thread_id, const LEVELSET_3D& from_input, LEVELSET_3D& to_input);
void DownSampling(MULTITHREADING* multithreading, const int& thread_id, const LEVELSET_3D& from_input, LEVELSET_3D& to_input);
