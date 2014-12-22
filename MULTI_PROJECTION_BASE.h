#pragma once

#include "LEVELSET_3D.h"
#include "PROJECTION_3D.h"

class MULTI_PROJECTION_BASE
{
public: // Essential Data
	MULTITHREADING*			multithreading;
	bool					place_dirichlet_bc_at_face;
	
	int						projection_type;
	int&					frame;

public: // Constructor and Destructor
	MULTI_PROJECTION_BASE(MULTITHREADING* multithreading_input, int& frame_input)
		: frame(frame_input), place_dirichlet_bc_at_face(false)
	{
		multithreading = multithreading_input;
	}

	~MULTI_PROJECTION_BASE(void)
	{}

public: // Virtual Funtions
	virtual void ExtrapolateFullCells(const int& thread_id, FIELD_STRUCTURE_3D<T>* scalar_field_input, FIELD_STRUCTURE_3D<T>* scalar_field_tmp_input, FIELD_STRUCTURE_3D<int>* bc_field_input, LEVELSET_3D* object_levelset_input)
	{
		// Note : Currently, face Dirichlet version make multithreading artifacts.
		if (place_dirichlet_bc_at_face)
		{
			ExtrapolateFullCellsFaceDirichlet(thread_id, scalar_field_input, scalar_field_tmp_input, bc_field_input, object_levelset_input);
		}
		else
		{
			ExtrapolateFullCellsNodeDirichlet(thread_id, scalar_field_input, scalar_field_tmp_input, bc_field_input, object_levelset_input);
		}
	}

	virtual void ExtrapolateFullCellsFaceDirichlet(const int& thread_id, FIELD_STRUCTURE_3D<T>* scalar_field_input, FIELD_STRUCTURE_3D<T>* scalar_field_tmp_input, FIELD_STRUCTURE_3D<int>* bc_field_input, LEVELSET_3D* object_levelset_input)
	{
		// Note : In this function, DIR boundary condition is applied on the face of the cell
		// Note : this method should be called after down-sampling and before upsampling
		// Before upsampling, pressure in BC_DIR, BC_OBJ should be correct values. (including GHOST CELLS)
		// BC_DIR cells are evaluated according to the Nuemann scheme.
		// BC_OBJ cells are set by given constant value such as zero (in the case of zero Dirichlet boundary condition).
		// Without correct pressure values in BC_OBJ and BC_DIR cells, interpolated pressure values near those cells after upsampling will be incoreect one.
		// Note : Role of this function is similar to FillGhostCellsFrom() but this function considers non-ghost cells as well

		T scalar_ave;
		int num_neighbor;

		// TODO : Check whether below is correct approach or not
		scalar_field_tmp_input->CopyAllValuesFrom(thread_id, *scalar_field_input);

		GRID_STRUCTURE_3D& grid(bc_field_input->grid);

		ARRAY_3D<int>& bc_arr(bc_field_input->array_for_this);

		ARRAY_3D<T>& scalar_arr(scalar_field_input->array_for_this);
		ARRAY_3D<T>& scalar_tmp_arr(scalar_field_tmp_input->array_for_this);

		ARRAY_3D<T>& obj_phi_arr(object_levelset_input->arr);

		if (projection_type == FREE_SURFACE_WATER)
		{
			BEGIN_GRID_ITERATION_3D(bc_field_input->partial_grids_ghost[thread_id])
			{
				if (bc_arr(i, j, k) >= BC_FULL)	// Note : After BuildLinearSystem(), All BC_FULL is changed to its sequence number
				{
					continue;
				}
				else if (i < grid.i_start)		// Wall boundary condition, TODO : Reflecting wall boundary condition in script
				{
					// NOTE : Wall Boundary is also BC_OBJ, so wall boundary handling should be come first
					// For wall boundary, we just copy outermost values. Therefore, it's distinguished from BC_OBJ handling which uses average value of nearby FULL cells.
					// and we need to check boundary range instead of just blindly using "bc == BC_OBJ".
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_start, j, k);
				}
				else if (i > grid.i_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_end, j, k);
				}
				else if (j < grid.j_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, grid.j_start, k);
				}
				else if (j > grid.j_end)
				{
					// Note : Ceiling is also BC_DIR, so wall boundary handling should be come first.
					scalar_arr(i, j, k) = -scalar_tmp_arr(i, grid.j_end, k);	// Face Dirichlet Boundary Condition
				}
				else if (k < grid.k_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_start);
				}
				else if (k > grid.k_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_end);
				}// End of Wall Boundary Condition
				else if (bc_arr(i, j, k) == BC_DIR)		// Air part
				{
					scalar_ave = (T)0;
					num_neighbor = 0;

					if(bc_arr(i-1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k-1); num_neighbor++;}
					if(bc_arr(i-1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k  ); num_neighbor++;}
					if(bc_arr(i-1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k+1); num_neighbor++;}
					if(bc_arr(i-1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k-1); num_neighbor++;}
					if(bc_arr(i-1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k  ); num_neighbor++;}
					if(bc_arr(i-1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k+1); num_neighbor++;}
					if(bc_arr(i-1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k-1); num_neighbor++;}
					if(bc_arr(i-1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k  ); num_neighbor++;}
					if(bc_arr(i-1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k+1); num_neighbor++;}
				
					if(bc_arr(i  ,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k-1); num_neighbor++;}
					if(bc_arr(i  ,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k  ); num_neighbor++;}
					if(bc_arr(i  ,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k+1); num_neighbor++;}
					if(bc_arr(i  ,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k-1); num_neighbor++;}
					if(bc_arr(i  ,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k  ); num_neighbor++;}
					if(bc_arr(i  ,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k+1); num_neighbor++;}
					if(bc_arr(i  ,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k-1); num_neighbor++;}
					if(bc_arr(i  ,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k  ); num_neighbor++;}
					if(bc_arr(i  ,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k+1); num_neighbor++;}

					if(bc_arr(i+1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k-1); num_neighbor++;}
					if(bc_arr(i+1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k  ); num_neighbor++;}
					if(bc_arr(i+1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k+1); num_neighbor++;}
					if(bc_arr(i+1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k-1); num_neighbor++;}
					if(bc_arr(i+1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k  ); num_neighbor++;}
					if(bc_arr(i+1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k+1); num_neighbor++;}
					if(bc_arr(i+1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k-1); num_neighbor++;}
					if(bc_arr(i+1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k  ); num_neighbor++;}
					
					// A BC_DIR cell adjacent to BC_FULL cells
					if (num_neighbor != 0)
					{
						scalar_arr(i, j, k) = -scalar_ave/(T)num_neighbor;
					}
					else
					{
						scalar_arr(i, j, k) = (T)0;
					}
				}
				else if (obj_phi_arr(i, j, k) <= (T)0)		// OBJ (Not Wall)
				{
					scalar_ave = (T)0;
					num_neighbor = 0;

					if(bc_arr(i-1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k-1); num_neighbor++;}
					if(bc_arr(i-1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k  ); num_neighbor++;}
					if(bc_arr(i-1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k+1); num_neighbor++;}
					if(bc_arr(i-1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k-1); num_neighbor++;}
					if(bc_arr(i-1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k  ); num_neighbor++;}
					if(bc_arr(i-1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k+1); num_neighbor++;}
					if(bc_arr(i-1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k-1); num_neighbor++;}
					if(bc_arr(i-1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k  ); num_neighbor++;}
					if(bc_arr(i-1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k+1); num_neighbor++;}
				
					if(bc_arr(i  ,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k-1); num_neighbor++;}
					if(bc_arr(i  ,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k  ); num_neighbor++;}
					if(bc_arr(i  ,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k+1); num_neighbor++;}
					if(bc_arr(i  ,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k-1); num_neighbor++;}
					if(bc_arr(i  ,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k  ); num_neighbor++;}
					if(bc_arr(i  ,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k+1); num_neighbor++;}
					if(bc_arr(i  ,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k-1); num_neighbor++;}
					if(bc_arr(i  ,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k  ); num_neighbor++;}
					if(bc_arr(i  ,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k+1); num_neighbor++;}

					if(bc_arr(i+1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k-1); num_neighbor++;}
					if(bc_arr(i+1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k  ); num_neighbor++;}
					if(bc_arr(i+1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k+1); num_neighbor++;}
					if(bc_arr(i+1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k-1); num_neighbor++;}
					if(bc_arr(i+1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k  ); num_neighbor++;}
					if(bc_arr(i+1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k+1); num_neighbor++;}
					if(bc_arr(i+1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k-1); num_neighbor++;}
					if(bc_arr(i+1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k  ); num_neighbor++;}
					
					// A BC_DIR cell adjacent to BC_FULL cells
					if (num_neighbor != 0)
					{
						scalar_arr(i, j, k) = scalar_ave/(T)num_neighbor;
					}
					else
					{
						scalar_arr(i, j, k) = (T)0;
					}
				}
			}
			END_GRID_ITERATION_3D;
		}
		else if (projection_type == AIR)
		{
			BEGIN_GRID_ITERATION_3D(bc_field_input->partial_grids_ghost[thread_id])
			{
				if (bc_arr(i, j, k) >= BC_FULL)
				{
					continue;
				}
				else if (i < grid.i_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_start, j, k);
				}
				else if (i > grid.i_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_end, j, k);
				}
				else if (j < grid.j_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, grid.j_start, k);
				}
				else if (j > grid.j_end)
				{
					scalar_arr(i, j, k) = -scalar_tmp_arr(i, grid.j_end, k);
				}
				else if (k < grid.k_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_start);
				}
				else if (k > grid.k_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_end);
				}// End of wall boundary condition
				else if (bc_arr(i, j, k) == BC_OBJ)			// OBJ + WALL + WATER
				{
					scalar_ave = (T)0;
					num_neighbor = 0;

					if(bc_arr(i-1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k-1); num_neighbor++;}
 					if(bc_arr(i-1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k  ); num_neighbor++;}
 					if(bc_arr(i-1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k+1); num_neighbor++;}
 					if(bc_arr(i-1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k-1); num_neighbor++;}
 					if(bc_arr(i-1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k  ); num_neighbor++;}
 					if(bc_arr(i-1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k+1); num_neighbor++;}
 					if(bc_arr(i-1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k-1); num_neighbor++;}
 					if(bc_arr(i-1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k  ); num_neighbor++;}
 					if(bc_arr(i-1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k+1); num_neighbor++;}
 							
 					if(bc_arr(i  ,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k-1); num_neighbor++;}
 					if(bc_arr(i  ,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k  ); num_neighbor++;}
 					if(bc_arr(i  ,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k+1); num_neighbor++;}
 					if(bc_arr(i  ,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k-1); num_neighbor++;}
 					if(bc_arr(i  ,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k  ); num_neighbor++;}
 					if(bc_arr(i  ,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k+1); num_neighbor++;}
 					if(bc_arr(i  ,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k-1); num_neighbor++;}
 					if(bc_arr(i  ,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k  ); num_neighbor++;}
 					if(bc_arr(i  ,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k+1); num_neighbor++;}
 			
 					if(bc_arr(i+1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k-1); num_neighbor++;}
 					if(bc_arr(i+1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k  ); num_neighbor++;}
 					if(bc_arr(i+1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k+1); num_neighbor++;}
 					if(bc_arr(i+1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k-1); num_neighbor++;}
 					if(bc_arr(i+1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k  ); num_neighbor++;}
 					if(bc_arr(i+1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k+1); num_neighbor++;}
 					if(bc_arr(i+1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k-1); num_neighbor++;}
 					if(bc_arr(i+1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k  ); num_neighbor++;}
 					if(bc_arr(i+1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k+1); num_neighbor++;}	
				
					// Interfacial cell between BC_OBJ and BC_FULL
					if (num_neighbor != 0)
					{
						scalar_arr(i, j, k) = scalar_ave/(T)num_neighbor;
					}
					else
					{
						scalar_arr(i, j, k) = (T)0;
					}
				}
			}
			END_GRID_ITERATION_3D;
		}
	}

	virtual void ExtrapolateFullCellsNodeDirichlet(const int& thread_id, FIELD_STRUCTURE_3D<T>* scalar_field_input, FIELD_STRUCTURE_3D<T>* scalar_field_tmp_input, FIELD_STRUCTURE_3D<int>* bc_field_input, LEVELSET_3D* object_levelset_input)
	{
		// Note : In this function, DIR boundary condition is applied on the face of the cell
		// Note : this method should be called after down-sampling and before upsampling
		// Before upsampling, pressure in BC_DIR, BC_OBJ should be correct values. (including GHOST CELLS)
		// BC_DIR cells are evaluated according to the Nuemann scheme.
		// BC_OBJ cells are set by given constant value such as zero (in the case of zero Dirichlet boundary condition).
		// Without correct pressure values in BC_OBJ and BC_DIR cells, interpolated pressure values near those cells after upsampling will be incoreect one.
		// Note : Role of this function is similar to FillGhostCellsFrom() but this function considers non-ghost cells as well

		T scalar_ave;
		int num_neighbor;

		// TODO : Check whether below is correct approach or not
		scalar_field_tmp_input->CopyAllValuesFrom(thread_id, *scalar_field_input);

		GRID_STRUCTURE_3D& grid(bc_field_input->grid);

		ARRAY_3D<int>& bc_arr(bc_field_input->array_for_this);

		ARRAY_3D<T>& scalar_arr(scalar_field_input->array_for_this);
		ARRAY_3D<T>& scalar_tmp_arr(scalar_field_tmp_input->array_for_this);

		ARRAY_3D<T>& obj_phi_arr(object_levelset_input->arr);

		if (projection_type == FREE_SURFACE_WATER)
		{
			BEGIN_GRID_ITERATION_3D(bc_field_input->partial_grids_ghost[thread_id])
			{
				if (bc_arr(i, j, k) >= BC_FULL)	// Note : After BuildLinearSystem(), All BC_FULL is changed to its sequence number
				{
					continue;
				}
				else if (bc_arr(i, j, k) == BC_DIR)
				{
					scalar_arr(i, j, k) = 0;	// Node Dirichlet Boundary Condition
				}
				else if (i < grid.i_start)		// Wall boundary condition, TODO : Reflecting wall boundary condition in script
				{
					// NOTE : Wall Boundary is also BC_OBJ, so wall boundary handling should be come first
					// For wall boundary, we just copy outermost values. Therefore, it's distinguished from BC_OBJ handling which uses average value of nearby FULL cells.
					// and we need to check boundary range instead of just blindly using "bc == BC_OBJ".
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_start, j, k);
				}
				else if (i > grid.i_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_end, j, k);
				}
				else if (j < grid.j_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, grid.j_start, k);
				}
				else if (k < grid.k_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_start);
				}
				else if (k > grid.k_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_end);
				}// End of Wall Boundary Condition
				else if (obj_phi_arr(i, j, k) <= (T)0)		// OBJ (Not Wall)
				{
					scalar_ave = (T)0;
					num_neighbor = 0;

					if(bc_arr(i-1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k-1); num_neighbor++;}
					if(bc_arr(i-1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k  ); num_neighbor++;}
					if(bc_arr(i-1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k+1); num_neighbor++;}
					if(bc_arr(i-1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k-1); num_neighbor++;}
					if(bc_arr(i-1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k  ); num_neighbor++;}
					if(bc_arr(i-1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k+1); num_neighbor++;}
					if(bc_arr(i-1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k-1); num_neighbor++;}
					if(bc_arr(i-1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k  ); num_neighbor++;}
					if(bc_arr(i-1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k+1); num_neighbor++;}
				
					if(bc_arr(i  ,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k-1); num_neighbor++;}
					if(bc_arr(i  ,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k  ); num_neighbor++;}
					if(bc_arr(i  ,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k+1); num_neighbor++;}
					if(bc_arr(i  ,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k-1); num_neighbor++;}
					if(bc_arr(i  ,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k  ); num_neighbor++;}
					if(bc_arr(i  ,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k+1); num_neighbor++;}
					if(bc_arr(i  ,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k-1); num_neighbor++;}
					if(bc_arr(i  ,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k  ); num_neighbor++;}
					if(bc_arr(i  ,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k+1); num_neighbor++;}

					if(bc_arr(i+1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k-1); num_neighbor++;}
					if(bc_arr(i+1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k  ); num_neighbor++;}
					if(bc_arr(i+1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k+1); num_neighbor++;}
					if(bc_arr(i+1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k-1); num_neighbor++;}
					if(bc_arr(i+1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k  ); num_neighbor++;}
					if(bc_arr(i+1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k+1); num_neighbor++;}
					if(bc_arr(i+1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k-1); num_neighbor++;}
					if(bc_arr(i+1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k  ); num_neighbor++;}
					
					// A BC_DIR cell adjacent to BC_FULL cells
					if (num_neighbor != 0)
					{
						scalar_arr(i, j, k) = scalar_ave/(T)num_neighbor;
					}
					else
					{
						scalar_arr(i, j, k) = (T)0;
					}
				}
			}
			END_GRID_ITERATION_3D;
		}
		else if (projection_type == AIR)
		{
			BEGIN_GRID_ITERATION_3D(bc_field_input->partial_grids_ghost[thread_id])
			{
				if (bc_arr(i, j, k) >= BC_FULL)
				{
					continue;
				}
				else if (i < grid.i_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_start, j, k);
				}
				else if (i > grid.i_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(grid.i_end, j, k);
				}
				else if (j < grid.j_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, grid.j_start, k);
				}
				else if (j > grid.j_end)
				{
					scalar_arr(i, j, k) = (T)0;
				}
				else if (k < grid.k_start)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_start);
				}
				else if (k > grid.k_end)
				{
					scalar_arr(i, j, k) = scalar_tmp_arr(i, j, grid.k_end);
				}// End of wall boundary condition
				else if (bc_arr(i, j, k) == BC_OBJ)			// OBJ + WALL + WATER
				{
					scalar_ave = (T)0;
					num_neighbor = 0;

					if(bc_arr(i-1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k-1); num_neighbor++;}
 					if(bc_arr(i-1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k  ); num_neighbor++;}
 					if(bc_arr(i-1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j-1,k+1); num_neighbor++;}
 					if(bc_arr(i-1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k-1); num_neighbor++;}
 					if(bc_arr(i-1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k  ); num_neighbor++;}
 					if(bc_arr(i-1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j  ,k+1); num_neighbor++;}
 					if(bc_arr(i-1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k-1); num_neighbor++;}
 					if(bc_arr(i-1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k  ); num_neighbor++;}
 					if(bc_arr(i-1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i-1,j+1,k+1); num_neighbor++;}
 							
 					if(bc_arr(i  ,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k-1); num_neighbor++;}
 					if(bc_arr(i  ,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k  ); num_neighbor++;}
 					if(bc_arr(i  ,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j-1,k+1); num_neighbor++;}
 					if(bc_arr(i  ,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k-1); num_neighbor++;}
 					if(bc_arr(i  ,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k  ); num_neighbor++;}
 					if(bc_arr(i  ,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j  ,k+1); num_neighbor++;}
 					if(bc_arr(i  ,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k-1); num_neighbor++;}
 					if(bc_arr(i  ,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k  ); num_neighbor++;}
 					if(bc_arr(i  ,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i  ,j+1,k+1); num_neighbor++;}
 			
 					if(bc_arr(i+1,j-1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k-1); num_neighbor++;}
 					if(bc_arr(i+1,j-1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k  ); num_neighbor++;}
 					if(bc_arr(i+1,j-1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j-1,k+1); num_neighbor++;}
 					if(bc_arr(i+1,j,  k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k-1); num_neighbor++;}
 					if(bc_arr(i+1,j,  k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k  ); num_neighbor++;}
 					if(bc_arr(i+1,j,  k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j  ,k+1); num_neighbor++;}
 					if(bc_arr(i+1,j+1,k-1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k-1); num_neighbor++;}
 					if(bc_arr(i+1,j+1,k  ) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k  ); num_neighbor++;}
 					if(bc_arr(i+1,j+1,k+1) == BC_FULL){scalar_ave += scalar_tmp_arr(i+1,j+1,k+1); num_neighbor++;}	
				
					// Interfacial cell between BC_OBJ and BC_FULL
					if (num_neighbor != 0)
					{
						scalar_arr(i, j, k) = scalar_ave/(T)num_neighbor;
					}
					else
					{
						scalar_arr(i, j, k) = (T)0;
					}
				}
			}
			END_GRID_ITERATION_3D;
		}
	}

	virtual void FillNonFullCellsWithZero(const int& thread_id, FIELD_STRUCTURE_3D<T>* scalar_field_input, FIELD_STRUCTURE_3D<int>* bc_field_input)
	{
		GRID_STRUCTURE_3D& grid(bc_field_input->grid);
		ARRAY_3D<int>& bc_arr(bc_field_input->array_for_this);
		ARRAY_3D<T>& scalar_arr(scalar_field_input->array_for_this);

		BEGIN_GRID_ITERATION_3D(bc_field_input->partial_grids_ghost[thread_id])
		{
			if (bc_arr(i, j, k) >= BC_FULL)
			{
				continue;
			}
			else
			{
				scalar_arr(i, j, k) = 0;
			}
		}
		END_GRID_ITERATION_3D;
	}
};