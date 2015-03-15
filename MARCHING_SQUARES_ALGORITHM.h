#pragma once

#include "COMMON_DEFINITIONS.h"
#include "MULTITHREADING.h"
#include "GRID_STRUCTURE_2D.h"
#include "FIELD_STRUCTURE_2D.h"
#include "ARRAY_2D.h"
#include "DYNAMIC_ARRAY.h"
#include "MARCHING_SQUARES_ALGORITHM_TABLE.h"

class MARCHING_SQUARES_ALGORITHM
{
public: // Essential Data
    MULTITHREADING*             multithreading;
    GRID_STRUCTURE_2D           grid;
    ARRAY<GRID_STRUCTURE_2D>    partial_grids;

    ARRAY<DYNAMIC_ARRAY<VT2>>   partial_vertex_array;
    ARRAY<DYNAMIC_ARRAY<VT2>>   contour_vertex_array;

public: // Constructor and Destructor
    MARCHING_SQUARES_ALGORITHM(void)
		: multithreading(0)
	{}

    ~MARCHING_SQUARES_ALGORITHM(void)
	{}

public: // Initialization Function
    void Initialize(MULTITHREADING* multithreading_input, const GRID_STRUCTURE_2D& grid_input, const int num_contour = 1)
	{
	    multithreading = multithreading_input;
        grid = grid_input;
        
        int num_threads = multithreading->num_threads;

        grid.SplitInHeight(num_threads, partial_grids);

        partial_vertex_array.Initialize(num_threads);

		for (int i = 0; i < num_threads; i++)
		{
            partial_vertex_array[i].Initialize(grid.i_res, grid.i_res);
		}

        contour_vertex_array.Initialize(num_contour);

		for (int i = 0; i < num_contour; i++)
		{
            contour_vertex_array[i].Initialize(grid.i_res, grid.i_res);
		}
	}

public: // Member Functions
    void EdgePoint(VT2 points[4], const T phi_values[4], FIELD_STRUCTURE_2D<T>& arr, const int& i, const int& j, const T& iso_value)
	{
        const VT2 cell_center = grid.CellCenter(i, j);

        const T* v = phi_values;

        ARRAY<VT2> p(4);

        p[0] = grid.CellCenter(i    , j + 1);
        p[1] = grid.CellCenter(i + 1, j + 1);
        p[2] = grid.CellCenter(i + 1, j    );
        p[3] = grid.CellCenter(i    , j    );

		for (int k = 0; k < 4; k++)
		{
			T theta = (iso_value - v[k])/(v[(k+1)%4] - v[k]);
            points[k] = p[k]*((T)1 - theta) + p[(k+1)%4]*theta;
		}
	}

    void Polygonize(const int& thread_id, FIELD_STRUCTURE_2D<T>& arr, const int contour_index, const T iso_value)
	{
	    int num_threads;

        BEGIN_HEAD_THREAD_WORK
		{
		    num_threads = multithreading->num_threads;
			for (int i = 0; i < num_threads; i++)
			{
                partial_vertex_array[i].Emptify();
			}
            contour_vertex_array[contour_index].Emptify();
		}
        END_HEAD_THREAD_WORK;

        VT2 vertices[4];
        T   phi_values[4];

        const GRID_STRUCTURE_2D& grid_p = partial_grids[thread_id];
        DYNAMIC_ARRAY<VT2>& vertex_array = partial_vertex_array[thread_id];

        int i_start = grid_p.i_start, i_end = grid_p.i_end;
        int j_start = grid_p.j_start, j_end = grid_p.j_end;

        i_end--;
		if (thread_id == num_threads - 1)
		{
            j_end--;
		}

		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_start; i <= i_end; i++)
			{
                phi_values[0] = arr(i    , j + 1);
                phi_values[1] = arr(i + 1, j + 1);
                phi_values[2] = arr(i + 1, j    );
                phi_values[3] = arr(i    , j    );

                int square_index(0);

				if (phi_values[0] < iso_value)
				{
                    square_index |= 8;
				}
                if (phi_values[1] < iso_value)
				{
                    square_index |= 4;
				}
                if (phi_values[2] < iso_value)
				{
                    square_index |= 2;
				}
                if (phi_values[3] < iso_value)
				{
                    square_index |= 1;
				}

				if (square_index > 0)
				{
                    EdgePoint(vertices, phi_values, arr, i, j, iso_value);
				}

                int edge_index = MS_TABLES::MS_EDGE_TABLE[square_index];

				if (edge_index != -1)
				{
                    int* contour_index = MS_TABLES::MS_CONTUOR_LINE_TABLE[edge_index];

					for (int k = 0; k < 4; k++)
					{
                        int ix = contour_index[k];
						if (ix != -1)
						{
                            vertex_array.Push(vertices[ix]);
						}
					}
				}
				else
				{
                    T pc = arr.BiLinearInterpolation((vertices[0] + vertices[1] + vertices[2] + vertices[3])*(T)0.25);

                    if(square_index == 5 && pc < (T)iso_value)
				    {
				    	vertex_array.Push(vertices[0]); 
						vertex_array.Push(vertices[3]);
				    	vertex_array.Push(vertices[1]); 
						vertex_array.Push(vertices[2]);
    				}
                    if(square_index == 5 && pc > (T)iso_value)
			    	{
			    		vertex_array.Push(vertices[0]); 
						vertex_array.Push(vertices[1]);
			    		vertex_array.Push(vertices[2]); 
						vertex_array.Push(vertices[3]);
			    	}
				    if(square_index == 10 && pc < (T)iso_value)
				    {
				    	vertex_array.Push(vertices[0]); 
						vertex_array.Push(vertices[1]);
					    vertex_array.Push(vertices[2]); 
						vertex_array.Push(vertices[3]);
				    }
				    if(square_index == 10 && pc > (T)iso_value)
				    {
				    	vertex_array.Push(vertices[0]); 
						vertex_array.Push(vertices[3]);
					    vertex_array.Push(vertices[1]);	
						vertex_array.Push(vertices[2]);
			    	}
				}
			}
		}
        multithreading->Sync(thread_id);
	    
        BEGIN_HEAD_THREAD_WORK
		{
			for (int k = 0; k < partial_vertex_array.num_elements; k++)
			{
                for (int l = 0; l < vertex_array.num_of_elements; l++)
				{
                    contour_vertex_array[contour_index].Push((partial_vertex_array[k])[l]);
				}
			}
		}
        END_HEAD_THREAD_WORK;
	}

    DYNAMIC_ARRAY<VT2>* ContourVertexArray(const int contour_index = 0)
	{
	    return &(contour_vertex_array[contour_index]);
	}
};