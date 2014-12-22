#include "stdafx.h"
#include "LEVELSET_3D.h"

void Sampling(MULTITHREADING* multithreading, const int& thread_id, const LEVELSET_3D& from_input, LEVELSET_3D& to_input)
{
	// Note: Only levelset upsampling handles ghost cells, whereas FIELD_STRUCTURE_3D doesn't.
	ARRAY_3D<T>& to_arr(to_input.signed_distance_field.array_for_this);

	BEGIN_GRID_ITERATION_3D(to_input.signed_distance_field.partial_grids_ghost[thread_id])
	{
			to_arr(i, j, k) = from_input.TriLinearInterpolation(to_input.CellCenter(i, j, k));
	}
	END_GRID_ITERATION_3D;
}

void DownSampling(MULTITHREADING* multithreading, const int& thread_id, const LEVELSET_3D& from_input, LEVELSET_3D& to_input)
{
	// Note: This is a fast, optimized implementation of downsampling, called restriction
	//		 This is for 8 to 1 cell downsampling
	const ARRAY_3D<T>& from_arr(from_input.signed_distance_field.array_for_this);
	ARRAY_3D<T>& to_arr(to_input.signed_distance_field.array_for_this);

	GRID_ITERATION_3D(to_input.partial_grids[thread_id])
	{
		to_arr(i, j, k) = ( from_arr(i*2,     j*2,     k*2    )
						  +	from_arr(i*2,     j*2,     k*2 + 1)
						  + from_arr(i*2,     j*2 + 1, k*2    )
						  + from_arr(i*2,     j*2 + 1, k*2 + 1)
						  + from_arr(i*2 + 1, j*2    , k*2    )
						  + from_arr(i*2 + 1, j*2 ,    k*2 + 1)
						  + from_arr(i*2 + 1, j*2 + 1, k*2    )
						  + from_arr(i*2 + 1, j*2 + 1, k*2 + 1)) / (T)8;
	}
	multithreading->Sync(thread_id);
}