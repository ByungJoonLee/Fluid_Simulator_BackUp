#include "stdafx.h"
#include "GRID_STRUCTURE_3D.h"

void GRID_STRUCTURE_3D::SplitInXDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids)
{
	partial_grids.Initialize(num_threads);

	const int quotient = i_res / num_threads;
	const int remainder = i_res % num_threads;

	int i_start_p = i_start;
	T x_min_p = x_min;

	for(int i = 0; i < num_threads; i++)
	{
		const int i_res_p = (i < remainder) ? (quotient + 1) : (quotient);
		const T x_max_p = dx*(T)i_res_p + x_min_p;
		partial_grids[i].Initialize(i_res_p, j_res, k_res, i_start_p, j_start, k_start, x_min_p, y_min, z_min, x_max_p, y_max, z_max);
		i_start_p += i_res_p;
		x_min_p = x_max_p;
	}
}

void GRID_STRUCTURE_3D::SplitInYDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids)
{
	partial_grids.Initialize(num_threads);

	const int quotient = j_res / num_threads;
	const int remainder = j_res % num_threads;

	int j_start_p = j_start;
	T y_min_p = y_min;

	for(int j = 0; j < num_threads; j++)
	{
		const int j_res_p = (j < remainder) ? (quotient + 1) : (quotient);
		const T y_max_p = dy*(T)j_res_p + y_min_p;
		partial_grids[j].Initialize(i_res, j_res_p, k_res, i_start, j_start_p, k_start, x_min, y_min_p, z_min, x_max, y_max_p, z_max);
		j_start_p += j_res_p;
		y_min_p = y_max_p;
	}
}

void GRID_STRUCTURE_3D::SplitInZDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids)
{
	partial_grids.Initialize(num_threads);

	const int quotient = k_res / num_threads;
	const int remainder = k_res % num_threads;
	
	int k_start_p = k_start;
	T z_min_p = z_min;

	for(int k = 0; k < num_threads; k++)
	{
		const int k_res_p = (k < remainder) ? (quotient + 1) : (quotient);
		const T z_max_p = dz*(T)k_res_p + z_min_p;
		partial_grids[k].Initialize(i_res, j_res, k_res_p, i_start, j_start, k_start_p, x_min, y_min, z_min_p, x_max, y_max, z_max_p);
		k_start_p += k_res_p;
		z_min_p = z_max_p;
	}
}

void GRID_STRUCTURE_3D::SplitInMaxDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids)
{
	if(num_threads <= k_res) SplitInZDirection(num_threads, partial_grids);
	else if(k_res >= j_res && k_res >= i_res) SplitInZDirection(num_threads, partial_grids);
	else if(j_res >= i_res) SplitInYDirection(num_threads, partial_grids);
	else SplitInXDirection(num_threads, partial_grids);
}

void GRID_STRUCTURE_3D::SplitInZDierctionPartially(const int& num_threads, ARRAY<GRID_STRUCTURE_3D>& partial_grids)
{
	partial_grids.Initialize(num_threads);

	int i_idx = i_end - i_start + 1;
	int j_idx = j_end - j_start + 1;
	int k_idx = k_end - k_start + 1;

	const int quotient = k_idx / num_threads;
	const int remainder = k_idx % num_threads;

	int k_start_p = k_start;
	T z_min_p = z_min;

	for(int i = 0; i < num_threads; i++)
	{
		const int k_res_p = (i < remainder) ? (quotient + 1) : (quotient);
		const T z_max_p = dz*(T)k_res_p + z_min_p;
		partial_grids[i].Initialize(i_idx, j_idx, k_res_p, i_start, j_start, k_start, x_min, y_min, z_min_p, x_max, y_max, z_max_p);
		k_start_p += k_res_p;
		z_min_p = z_max_p;
	}
}

int GRID_STRUCTURE_3D::RecommendMaxMultigridLevel(const int& lowest_level_res)
{
	int min_res = MIN3(i_res, j_res, k_res);
	int max_level = 0;
	while(true)
	{
		max_level++;
		min_res /= 2;
		if(min_res < lowest_level_res)
			break;
	}

	return max_level;
}




