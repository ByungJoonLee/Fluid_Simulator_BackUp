#pragma once

#include "LEVELSET_3D.h"
// #include "NETWORKING_NODE.h"

template<class TT>
class ADVECTION_METHOD_3D
{
public: 
	static void SL1stOrder(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<VT>& velocity, const T& dt, MULTITHREADING& multithreading)
	{
		const ARRAY_3D<VT>& velocity_array(velocity.array_for_this);
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(thread_id, rho_array, true);

		multithreading.Sync(thread_id);

		GRID_ITERATION_3D(rho.partial_grids[thread_id])
			rho_array(i, j, k) = rho_ghost.TriLinearInterpolation(rho.CellCenter(i, j, k) - velocity_array(i, j, k)*dt);
		
		multithreading.Sync(thread_id);
	}

	static void SL1stOrder(const int& thread_id, LEVELSET_3D& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<VT>& velocity, const T& dt, MULTITHREADING& multithreading)
	{
		const ARRAY_3D<VT>& velocity_array(velocity.array_for_this);
		ARRAY_3D<TT>& rho_array(rho.arr);

		rho_ghost.FillGhostCellsFrom(thread_id, rho_array, true);

		multithreading.Sync(thread_id);

		int i, j, k;
		INIT_GRIDRANGE_3D(rho.partial_grids[thread_id], i_start, j_start, k_start, i_end, j_end, k_end);
		LOOPS_3D(i, j, k, i_start - 1, j_start - 1, k_start - 1, i_end + 1, j_end + 1, k_end + 1)
			rho_array(i, j, k) = rho_ghost.TriLinearInterpolation(rho.CellCenter(i, j, k) - velocity_array(i, j, k)*dt);

		multithreading.Sync(thread_id);
	}

	static void WENO5th(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<TT>& velocity_x, const FIELD_STRUCTURE_3D<TT>& velocity_y, const FIELD_STRUCTURE_3D<TT>& velocity_z, const T& dt, MULTITHREADING& multithreading, const T& epsilon)
	{
		const ARRAY_3D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_3D<TT>& velocity_array_y(velocity_y.array_for_this);
		const ARRAY_3D<TT>& velocity_array_z(velocity_z.array_for_this);
		
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			// x-components
			TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
			TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
			TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));

			TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));
			TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
			TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_p_x_4 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));

			// y-components
			TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
			TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
			TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));

			TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));
			TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
			TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));

			// z-components
			TT rho_m_z_1 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
			TT rho_m_z_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
			TT rho_m_z_3 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_m_z_4 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_m_z_5 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));

			TT rho_p_z_1 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));
			TT rho_p_z_2 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
			TT rho_p_z_3 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_p_z_4 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_p_z_5 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));

			// Smoothness 
			// x-component
			TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
			TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
			TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

			TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
			TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
			TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

			// y-component
			TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
			TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
			TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

			TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
			TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
			TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

			// z-component
			TT s_m_z_1 = (T)13/12*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3)*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3) + (T)1/4*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3)*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3);
			TT s_m_z_2 = (T)13/12*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4)*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4) + (T)1/4*(rho_m_z_2 - rho_m_z_4)*(rho_m_z_2 - rho_m_z_4);
			TT s_m_z_3 = (T)13/12*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5)*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5) + (T)1/4*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5)*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5);

			TT s_p_z_1 = (T)13/12*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3)*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3) + (T)1/4*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3)*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3);
			TT s_p_z_2 = (T)13/12*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4)*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4) + (T)1/4*(rho_p_z_2 - rho_p_z_4)*(rho_p_z_2 - rho_p_z_4);
			TT s_p_z_3 = (T)13/12*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5)*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5) + (T)1/4*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5)*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5);

			// Weights
			// x-component
			TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
			TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
			TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

			TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

			TT w_m_x_1 = a_m_x_1/sum_m_x;
			TT w_m_x_2 = a_m_x_2/sum_m_x;
			TT w_m_x_3 = a_m_x_3/sum_m_x;

			TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
			TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
			TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

			TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

			TT w_p_x_1 = a_p_x_1/sum_p_x;
			TT w_p_x_2 = a_p_x_2/sum_p_x;
			TT w_p_x_3 = a_p_x_3/sum_p_x;

			// y-component
			TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
			TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
			TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

			TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

			TT w_m_y_1 = a_m_y_1/sum_m_y;
			TT w_m_y_2 = a_m_y_2/sum_m_y;
			TT w_m_y_3 = a_m_y_3/sum_m_y;

			TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
			TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
			TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

			TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

			TT w_p_y_1 = a_p_y_1/sum_p_y;
			TT w_p_y_2 = a_p_y_2/sum_p_y;
			TT w_p_y_3 = a_p_y_3/sum_p_y;

			// z-component
			TT a_m_z_1 = (T)1/10*(T)1/((epsilon + s_m_z_1)*(epsilon + s_m_z_1));
			TT a_m_z_2 = (T)6/10*(T)1/((epsilon + s_m_z_2)*(epsilon + s_m_z_2));
			TT a_m_z_3 = (T)3/10*(T)1/((epsilon + s_m_z_3)*(epsilon + s_m_z_3));

			TT sum_m_z = a_m_z_1 + a_m_z_2 + a_m_z_3;

			TT w_m_z_1 = a_m_z_1/sum_m_z;
			TT w_m_z_2 = a_m_z_2/sum_m_z;
			TT w_m_z_3 = a_m_z_3/sum_m_z;

			TT a_p_z_1 = (T)1/10*(T)1/((epsilon + s_p_z_1)*(epsilon + s_p_z_1));
			TT a_p_z_2 = (T)6/10*(T)1/((epsilon + s_p_z_2)*(epsilon + s_p_z_2));
			TT a_p_z_3 = (T)3/10*(T)1/((epsilon + s_p_z_3)*(epsilon + s_p_z_3));

			TT sum_p_z = a_p_z_1 + a_p_z_2 + a_p_z_3;

			TT w_p_z_1 = a_p_z_1/sum_p_z;
			TT w_p_z_2 = a_p_z_2/sum_p_z;
			TT w_p_z_3 = a_p_z_3/sum_p_z;

			// Approximation of derivatives
			// rho_x
			TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6); 
			TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6); 
			
			// rho_y
			TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6); 
			TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6); 
			
			// rho_z
			TT rho_m_z = w_m_z_1*(rho_m_z_1*(T)1/3 - rho_m_z_2*(T)7/6 + rho_m_z_3*(T)11/6) + w_m_z_2*(rho_m_z_2*(-(T)1/6) + rho_m_z_3*(T)5/6 + rho_m_z_4*(T)1/3) + w_m_z_3*(rho_m_z_3*(T)1/3 + rho_m_z_4*(T)5/6 - rho_m_z_5*(T)1/6); 
			TT rho_p_z = w_p_z_1*(rho_p_z_1*(T)1/3 - rho_p_z_2*(T)7/6 + rho_p_z_3*(T)11/6) + w_p_z_2*(rho_p_z_2*(-(T)1/6) + rho_p_z_3*(T)5/6 + rho_p_z_4*(T)1/3) + w_p_z_3*(rho_p_z_3*(T)1/3 + rho_p_z_4*(T)5/6 - rho_p_z_5*(T)1/6); 

			T u_vel, v_vel, w_vel;
			
			// Velocity interpolation - Using MAC grid
			if (rho.is_scalar)
			{
				u_vel = (velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k))*(T)0.5;
				v_vel = (velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k))*(T)0.5;
				w_vel = (velocity_array_z(i, j, k + 1) + velocity_array_z(i, j, k))*(T)0.5;
			}
			else if (rho.is_x_component)
			{
				u_vel = velocity_array_x(i, j, k);
				v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i - 1, j, k) + velocity_array_y(i - 1, j + 1, k))*(T)0.25;
				w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i - 1, j, k) + velocity_array_z(i - 1, j, k + 1))*(T)0.25;
			}
			else if (rho.is_y_component)
			{
				u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j - 1, k) + velocity_array_x(i + 1, j - 1, k))*(T)0.25;
				v_vel = velocity_array_y(i, j, k);
				w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i, j - 1, k) + velocity_array_z(i, j - 1, k + 1))*(T)0.25;
			}
			else if (rho.is_z_component)
			{
				u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k - 1) + velocity_array_x(i + 1, j, k - 1))*(T)0.25;
				v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k - 1) + velocity_array_y(i, j + 1, k - 1))*(T)0.25;
				w_vel = velocity_array_z(i, j, k);
			}
			
			if (u_vel > 0)
			{
				if (v_vel > 0)
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
					}
				}
				else
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
					}
				}
			}
			else
			{
				if (v_vel > 0)
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
					}
				}
				else
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
					}
				}
			}
		}
	
		multithreading.Sync(thread_id);
	}

	static void WENO5th(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<TT>& boundary_levelset, const FIELD_STRUCTURE_3D<TT>& velocity_x, const FIELD_STRUCTURE_3D<TT>& velocity_y, const FIELD_STRUCTURE_3D<TT>& velocity_z, const T& dt, MULTITHREADING& multithreading, const T& epsilon)
	{
		const ARRAY_3D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_3D<TT>& velocity_array_y(velocity_y.array_for_this);
		const ARRAY_3D<TT>& velocity_array_z(velocity_z.array_for_this);
		
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			if (boundary_levelset(i, j, k) > 0)
			{
				continue;
			}

			// x-components
			TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
			TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
			TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));

			TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));
			TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
			TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_p_x_4 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));

			// y-components
			TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
			TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
			TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));

			TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));
			TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
			TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));

			// z-components
			TT rho_m_z_1 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
			TT rho_m_z_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
			TT rho_m_z_3 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_m_z_4 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_m_z_5 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));

			TT rho_p_z_1 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));
			TT rho_p_z_2 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
			TT rho_p_z_3 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_p_z_4 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_p_z_5 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));

			// Smoothness 
			// x-component
			TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
			TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
			TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

			TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
			TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
			TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

			// y-component
			TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
			TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
			TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

			TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
			TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
			TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

			// z-component
			TT s_m_z_1 = (T)13/12*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3)*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3) + (T)1/4*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3)*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3);
			TT s_m_z_2 = (T)13/12*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4)*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4) + (T)1/4*(rho_m_z_2 - rho_m_z_4)*(rho_m_z_2 - rho_m_z_4);
			TT s_m_z_3 = (T)13/12*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5)*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5) + (T)1/4*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5)*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5);

			TT s_p_z_1 = (T)13/12*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3)*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3) + (T)1/4*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3)*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3);
			TT s_p_z_2 = (T)13/12*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4)*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4) + (T)1/4*(rho_p_z_2 - rho_p_z_4)*(rho_p_z_2 - rho_p_z_4);
			TT s_p_z_3 = (T)13/12*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5)*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5) + (T)1/4*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5)*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5);

			// Weights
			// x-component
			TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
			TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
			TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

			TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

			TT w_m_x_1 = a_m_x_1/sum_m_x;
			TT w_m_x_2 = a_m_x_2/sum_m_x;
			TT w_m_x_3 = a_m_x_3/sum_m_x;

			TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
			TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
			TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

			TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

			TT w_p_x_1 = a_p_x_1/sum_p_x;
			TT w_p_x_2 = a_p_x_2/sum_p_x;
			TT w_p_x_3 = a_p_x_3/sum_p_x;

			// y-component
			TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
			TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
			TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

			TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

			TT w_m_y_1 = a_m_y_1/sum_m_y;
			TT w_m_y_2 = a_m_y_2/sum_m_y;
			TT w_m_y_3 = a_m_y_3/sum_m_y;

			TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
			TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
			TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

			TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

			TT w_p_y_1 = a_p_y_1/sum_p_y;
			TT w_p_y_2 = a_p_y_2/sum_p_y;
			TT w_p_y_3 = a_p_y_3/sum_p_y;

			// z-component
			TT a_m_z_1 = (T)1/10*(T)1/((epsilon + s_m_z_1)*(epsilon + s_m_z_1));
			TT a_m_z_2 = (T)6/10*(T)1/((epsilon + s_m_z_2)*(epsilon + s_m_z_2));
			TT a_m_z_3 = (T)3/10*(T)1/((epsilon + s_m_z_3)*(epsilon + s_m_z_3));

			TT sum_m_z = a_m_z_1 + a_m_z_2 + a_m_z_3;

			TT w_m_z_1 = a_m_z_1/sum_m_z;
			TT w_m_z_2 = a_m_z_2/sum_m_z;
			TT w_m_z_3 = a_m_z_3/sum_m_z;

			TT a_p_z_1 = (T)1/10*(T)1/((epsilon + s_p_z_1)*(epsilon + s_p_z_1));
			TT a_p_z_2 = (T)6/10*(T)1/((epsilon + s_p_z_2)*(epsilon + s_p_z_2));
			TT a_p_z_3 = (T)3/10*(T)1/((epsilon + s_p_z_3)*(epsilon + s_p_z_3));

			TT sum_p_z = a_p_z_1 + a_p_z_2 + a_p_z_3;

			TT w_p_z_1 = a_p_z_1/sum_p_z;
			TT w_p_z_2 = a_p_z_2/sum_p_z;
			TT w_p_z_3 = a_p_z_3/sum_p_z;

			// Approximation of derivatives
			// rho_x
			TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6); 
			TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6); 
			
			// rho_y
			TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6); 
			TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6); 
			
			// rho_z
			TT rho_m_z = w_m_z_1*(rho_m_z_1*(T)1/3 - rho_m_z_2*(T)7/6 + rho_m_z_3*(T)11/6) + w_m_z_2*(rho_m_z_2*(-(T)1/6) + rho_m_z_3*(T)5/6 + rho_m_z_4*(T)1/3) + w_m_z_3*(rho_m_z_3*(T)1/3 + rho_m_z_4*(T)5/6 - rho_m_z_5*(T)1/6); 
			TT rho_p_z = w_p_z_1*(rho_p_z_1*(T)1/3 - rho_p_z_2*(T)7/6 + rho_p_z_3*(T)11/6) + w_p_z_2*(rho_p_z_2*(-(T)1/6) + rho_p_z_3*(T)5/6 + rho_p_z_4*(T)1/3) + w_p_z_3*(rho_p_z_3*(T)1/3 + rho_p_z_4*(T)5/6 - rho_p_z_5*(T)1/6); 

			T u_vel, v_vel, w_vel;
			
			// Velocity interpolation - Using MAC grid
			if (rho.is_scalar)
			{
				u_vel = (velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k))*(T)0.5;
				v_vel = (velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k))*(T)0.5;
				w_vel = (velocity_array_z(i, j, k + 1) + velocity_array_z(i, j, k))*(T)0.5;
			}
			else if (rho.is_x_component)
			{
				u_vel = velocity_array_x(i, j, k);
				v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i - 1, j, k) + velocity_array_y(i - 1, j + 1, k))*(T)0.25;
				w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i - 1, j, k) + velocity_array_z(i - 1, j, k + 1))*(T)0.25;
			}
			else if (rho.is_y_component)
			{
				u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j - 1, k) + velocity_array_x(i + 1, j - 1, k))*(T)0.25;
				v_vel = velocity_array_y(i, j, k);
				w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i, j - 1, k) + velocity_array_z(i, j - 1, k + 1))*(T)0.25;
			}
			else if (rho.is_z_component)
			{
				u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k - 1) + velocity_array_x(i + 1, j, k - 1))*(T)0.25;
				v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k - 1) + velocity_array_y(i, j + 1, k - 1))*(T)0.25;
				w_vel = velocity_array_z(i, j, k);
			}
			
			if (u_vel > 0)
			{
				if (v_vel > 0)
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
					}
				}
				else
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
					}
				}
			}
			else
			{
				if (v_vel > 0)
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
					}
				}
				else
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
					}
				}
			}
		}
	
		multithreading.Sync(thread_id);
	}
	// Use this only for updating the levelset method
	static void WENO5th(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<VT>& velocity, const T& dt, MULTITHREADING& multithreading, const T& epsilon, const bool& use_rk2 = false)
	{
		const ARRAY_3D<VT>& velocity_array(velocity.array_for_this);
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(thread_id, rho_array, true);

		multithreading.Sync(thread_id);
		
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			// x-components
			TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
			TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
			TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));

			TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));
			TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
			TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_p_x_4 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));

			// y-components
			TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
			TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
			TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));

			TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));
			TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
			TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));

			// z-components
			TT rho_m_z_1 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
			TT rho_m_z_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
			TT rho_m_z_3 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_m_z_4 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_m_z_5 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));

			TT rho_p_z_1 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));
			TT rho_p_z_2 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
			TT rho_p_z_3 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_p_z_4 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_p_z_5 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));

			// Smoothness 
			// x-component
			TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
			TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
			TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

			TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
			TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
			TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

			// y-component
			TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
			TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
			TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

			TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
			TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
			TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

			// z-component
			TT s_m_z_1 = (T)13/12*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3)*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3) + (T)1/4*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3)*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3);
			TT s_m_z_2 = (T)13/12*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4)*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4) + (T)1/4*(rho_m_z_2 - rho_m_z_4)*(rho_m_z_2 - rho_m_z_4);
			TT s_m_z_3 = (T)13/12*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5)*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5) + (T)1/4*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5)*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5);

			TT s_p_z_1 = (T)13/12*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3)*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3) + (T)1/4*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3)*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3);
			TT s_p_z_2 = (T)13/12*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4)*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4) + (T)1/4*(rho_p_z_2 - rho_p_z_4)*(rho_p_z_2 - rho_p_z_4);
			TT s_p_z_3 = (T)13/12*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5)*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5) + (T)1/4*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5)*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5);

			// Weights
			// x-component
			TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
			TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
			TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

			TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

			TT w_m_x_1 = a_m_x_1/sum_m_x;
			TT w_m_x_2 = a_m_x_2/sum_m_x;
			TT w_m_x_3 = a_m_x_3/sum_m_x;

			TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
			TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
			TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

			TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

			TT w_p_x_1 = a_p_x_1/sum_p_x;
			TT w_p_x_2 = a_p_x_2/sum_p_x;
			TT w_p_x_3 = a_p_x_3/sum_p_x;

			// y-component
			TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
			TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
			TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

			TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

			TT w_m_y_1 = a_m_y_1/sum_m_y;
			TT w_m_y_2 = a_m_y_2/sum_m_y;
			TT w_m_y_3 = a_m_y_3/sum_m_y;

			TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
			TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
			TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

			TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

			TT w_p_y_1 = a_p_y_1/sum_p_y;
			TT w_p_y_2 = a_p_y_2/sum_p_y;
			TT w_p_y_3 = a_p_y_3/sum_p_y;

			// z-component
			TT a_m_z_1 = (T)1/10*(T)1/((epsilon + s_m_z_1)*(epsilon + s_m_z_1));
			TT a_m_z_2 = (T)6/10*(T)1/((epsilon + s_m_z_2)*(epsilon + s_m_z_2));
			TT a_m_z_3 = (T)3/10*(T)1/((epsilon + s_m_z_3)*(epsilon + s_m_z_3));

			TT sum_m_z = a_m_z_1 + a_m_z_2 + a_m_z_3;

			TT w_m_z_1 = a_m_z_1/sum_m_z;
			TT w_m_z_2 = a_m_z_2/sum_m_z;
			TT w_m_z_3 = a_m_z_3/sum_m_z;

			TT a_p_z_1 = (T)1/10*(T)1/((epsilon + s_p_z_1)*(epsilon + s_p_z_1));
			TT a_p_z_2 = (T)6/10*(T)1/((epsilon + s_p_z_2)*(epsilon + s_p_z_2));
			TT a_p_z_3 = (T)3/10*(T)1/((epsilon + s_p_z_3)*(epsilon + s_p_z_3));

			TT sum_p_z = a_p_z_1 + a_p_z_2 + a_p_z_3;

			TT w_p_z_1 = a_p_z_1/sum_p_z;
			TT w_p_z_2 = a_p_z_2/sum_p_z;
			TT w_p_z_3 = a_p_z_3/sum_p_z;

			// Approximation of derivatives
			// rho_x
			TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6); 
			TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6); 
			
			// rho_y
			TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6); 
			TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6); 
			
			// rho_z
			TT rho_m_z = w_m_z_1*(rho_m_z_1*(T)1/3 - rho_m_z_2*(T)7/6 + rho_m_z_3*(T)11/6) + w_m_z_2*(rho_m_z_2*(-(T)1/6) + rho_m_z_3*(T)5/6 + rho_m_z_4*(T)1/3) + w_m_z_3*(rho_m_z_3*(T)1/3 + rho_m_z_4*(T)5/6 - rho_m_z_5*(T)1/6); 
			TT rho_p_z = w_p_z_1*(rho_p_z_1*(T)1/3 - rho_p_z_2*(T)7/6 + rho_p_z_3*(T)11/6) + w_p_z_2*(rho_p_z_2*(-(T)1/6) + rho_p_z_3*(T)5/6 + rho_p_z_4*(T)1/3) + w_p_z_3*(rho_p_z_3*(T)1/3 + rho_p_z_4*(T)5/6 - rho_p_z_5*(T)1/6); 

			// Velocity interpolation - Using MAC grid
			T u_vel = (velocity_array(i - 1, j, k).x + velocity_array(i + 1, j, k).x)*(T)0.5;
			T v_vel = (velocity_array(i, j - 1, k).y + velocity_array(i, j + 1, k).y)*(T)0.5;
			T w_vel = (velocity_array(i, j, k - 1).z + velocity_array(i, j, k + 1).z)*(T)0.5;

			if (!use_rk2)
			{
				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				}
				else
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				}
			}
		}

		multithreading.Sync(thread_id);
	}

	static void WENO5thReinitialization(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const T& dt, MULTITHREADING& multithreading, const T& epsilon, const FIELD_STRUCTURE_3D<T>& sign_function)
	{
		ARRAY_3D<TT>& rho_array(rho.array_for_this);
		
		//rho_ghost.FillGhostCellsContinuousDerivativesFrom(thread_id, rho.array_for_this, true);
		//rho_ghost.FillGhostCellsPeriodicInYDirection(thread_id, rho.array_for_this, true);

		multithreading.Sync(thread_id);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			// x-components
			TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
			TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
			TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));

			TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));
			TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
			TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i, j, k));
			TT rho_p_x_4 = one_over_dx*(rho_ghost(i, j, k) - rho_ghost(i - 1, j, k));
			TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));

			// y-components
			TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
			TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
			TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));

			TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));
			TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
			TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i, j, k));
			TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j, k) - rho_ghost(i, j - 1, k));
			TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));

			// z-components
			TT rho_m_z_1 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
			TT rho_m_z_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
			TT rho_m_z_3 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_m_z_4 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_m_z_5 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));

			TT rho_p_z_1 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));
			TT rho_p_z_2 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
			TT rho_p_z_3 = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i, j, k));
			TT rho_p_z_4 = one_over_dz*(rho_ghost(i, j, k) - rho_ghost(i, j, k - 1));
			TT rho_p_z_5 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));

			// Smoothness 
			// x-component
			TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
			TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
			TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

			TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
			TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
			TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

			// y-component
			TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
			TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
			TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

			TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
			TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
			TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

			// z-component
			TT s_m_z_1 = (T)13/12*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3)*(rho_m_z_1 - 2*rho_m_z_2 + rho_m_z_3) + (T)1/4*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3)*(rho_m_z_1 - 4*rho_m_z_2 + 3*rho_m_z_3);
			TT s_m_z_2 = (T)13/12*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4)*(rho_m_z_2 - 2*rho_m_z_3 + rho_m_z_4) + (T)1/4*(rho_m_z_2 - rho_m_z_4)*(rho_m_z_2 - rho_m_z_4);
			TT s_m_z_3 = (T)13/12*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5)*(rho_m_z_3 - 2*rho_m_z_4 + rho_m_z_5) + (T)1/4*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5)*(3*rho_m_z_3 - 4*rho_m_z_4 + rho_m_z_5);

			TT s_p_z_1 = (T)13/12*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3)*(rho_p_z_1 - 2*rho_p_z_2 + rho_p_z_3) + (T)1/4*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3)*(rho_p_z_1 - 4*rho_p_z_2 + 3*rho_p_z_3);
			TT s_p_z_2 = (T)13/12*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4)*(rho_p_z_2 - 2*rho_p_z_3 + rho_p_z_4) + (T)1/4*(rho_p_z_2 - rho_p_z_4)*(rho_p_z_2 - rho_p_z_4);
			TT s_p_z_3 = (T)13/12*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5)*(rho_p_z_3 - 2*rho_p_z_4 + rho_p_z_5) + (T)1/4*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5)*(3*rho_p_z_3 - 4*rho_p_z_4 + rho_p_z_5);

			// Weights
			// x-component
			TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
			TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
			TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

			TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

			TT w_m_x_1 = a_m_x_1/sum_m_x;
			TT w_m_x_2 = a_m_x_2/sum_m_x;
			TT w_m_x_3 = a_m_x_3/sum_m_x;

			TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
			TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
			TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

			TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

			TT w_p_x_1 = a_p_x_1/sum_p_x;
			TT w_p_x_2 = a_p_x_2/sum_p_x;
			TT w_p_x_3 = a_p_x_3/sum_p_x;

			// y-component
			TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
			TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
			TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

			TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

			TT w_m_y_1 = a_m_y_1/sum_m_y;
			TT w_m_y_2 = a_m_y_2/sum_m_y;
			TT w_m_y_3 = a_m_y_3/sum_m_y;

			TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
			TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
			TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

			TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

			TT w_p_y_1 = a_p_y_1/sum_p_y;
			TT w_p_y_2 = a_p_y_2/sum_p_y;
			TT w_p_y_3 = a_p_y_3/sum_p_y;

			// z-component
			TT a_m_z_1 = (T)1/10*(T)1/((epsilon + s_m_z_1)*(epsilon + s_m_z_1));
			TT a_m_z_2 = (T)6/10*(T)1/((epsilon + s_m_z_2)*(epsilon + s_m_z_2));
			TT a_m_z_3 = (T)3/10*(T)1/((epsilon + s_m_z_3)*(epsilon + s_m_z_3));

			TT sum_m_z = a_m_z_1 + a_m_z_2 + a_m_z_3;

			TT w_m_z_1 = a_m_z_1/sum_m_z;
			TT w_m_z_2 = a_m_z_2/sum_m_z;
			TT w_m_z_3 = a_m_z_3/sum_m_z;

			TT a_p_z_1 = (T)1/10*(T)1/((epsilon + s_p_z_1)*(epsilon + s_p_z_1));
			TT a_p_z_2 = (T)6/10*(T)1/((epsilon + s_p_z_2)*(epsilon + s_p_z_2));
			TT a_p_z_3 = (T)3/10*(T)1/((epsilon + s_p_z_3)*(epsilon + s_p_z_3));

			TT sum_p_z = a_p_z_1 + a_p_z_2 + a_p_z_3;

			TT w_p_z_1 = a_p_z_1/sum_p_z;
			TT w_p_z_2 = a_p_z_2/sum_p_z;
			TT w_p_z_3 = a_p_z_3/sum_p_z;

			// Approximation of derivatives
			// rho_x
			TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6); 
			TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6); 
			
			// rho_y
			TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6); 
			TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6); 
			
			// rho_z
			TT rho_m_z = w_m_z_1*(rho_m_z_1*(T)1/3 - rho_m_z_2*(T)7/6 + rho_m_z_3*(T)11/6) + w_m_z_2*(rho_m_z_2*(-(T)1/6) + rho_m_z_3*(T)5/6 + rho_m_z_4*(T)1/3) + w_m_z_3*(rho_m_z_3*(T)1/3 + rho_m_z_4*(T)5/6 - rho_m_z_5*(T)1/6); 
			TT rho_p_z = w_p_z_1*(rho_p_z_1*(T)1/3 - rho_p_z_2*(T)7/6 + rho_p_z_3*(T)11/6) + w_p_z_2*(rho_p_z_2*(-(T)1/6) + rho_p_z_3*(T)5/6 + rho_p_z_4*(T)1/3) + w_p_z_3*(rho_p_z_3*(T)1/3 + rho_p_z_4*(T)5/6 - rho_p_z_5*(T)1/6); 

			// Define the sign functions
			TT s = sign_function(i, j, k);

			// Define Sub functions
			TT smx = s*rho_m_x, spx = s*rho_p_x, smy = s*rho_m_y, spy = s*rho_p_y, smz = s*rho_m_z, spz = s*rho_p_z;

			// 1st case
			if ((spx >= 0) && (smx >= 0))
			{
				if ((spy >= 0) && (smy >= 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy <= 0) && (smy <= 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy > 0) && (smy < 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy < 0) && (smy > 0))
				{
					T sss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
					if (sss > 0)
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else 
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
				}
			}
			else if ((spx <= 0) && (smx <= 0))
			{
				if ((spy >= 0) && (smy >= 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy <= 0) && (smy <= 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy > 0) && (smy < 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy < 0) && (smy > 0))
				{
					T sss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
					if (sss > 0)
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else 
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
				}
			}
			else if ((spx > 0) && (smx < 0))
			{
				if ((spy >= 0) && (smy >= 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy <= 0) && (smy <= 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y)) - 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy > 0) && (smy < 0))
				{
					if ((spz >= 0) && (smz >= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_z)) - 1);
					}
					else if ((spz <= 0) && (smz <= 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_z)) - 1);
					}
					else if ((spz > 0) && (smz < 0))
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(- 1);
					}
					else if ((spz < 0) && (smz > 0))
					{
						T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
						if (ss > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_z)) - 1);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_z)) - 1);
						}
					}
				}
				else if ((spy < 0) && (smy > 0))
				{
					T sss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
					if (sss > 0)
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else 
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
				}
			}
			else if ((spx < 0) && (smx > 0))
			{
				T ssss = s*(abs(rho_p_x) - abs(rho_m_x))/(rho_p_x - rho_m_x);
				if (ssss > 0)
				{
					if ((spy >= 0) && (smy >= 0))
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else if ((spy > 0) && (smy < 0))
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else if ((spy < 0) && (smy > 0))
					{
						T sss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (sss > 0)
						{
							if ((spz >= 0) && (smz >= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else if ((spz <= 0) && (smz <= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
							else if ((spz > 0) && (smz < 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
							}
							else if ((spz < 0) && (smz > 0))
							{
								T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
								if (ss > 0)
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
								}
								else
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
								}
							}
						}
						else 
						{
							if ((spz >= 0) && (smz >= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else if ((spz <= 0) && (smz <= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
							else if ((spz > 0) && (smz < 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x)) + POW2(rho_p_y) - 1);
							}
							else if ((spz < 0) && (smz > 0))
							{
								T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
								if (ss > 0)
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
								}
								else
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
								}
							}
						}
					}
				}
				else
				{
					if ((spy >= 0) && (smy >= 0))
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else if ((spy > 0) && (smy < 0))
					{
						if ((spz >= 0) && (smz >= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_z)) - 1);
						}
						else if ((spz <= 0) && (smz <= 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_z)) - 1);
						}
						else if ((spz > 0) && (smz < 0))
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x)) - 1);
						}
						else if ((spz < 0) && (smz > 0))
						{
							T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
							if (ss > 0)
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_z)) - 1);
							}
							else
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_z)) - 1);
							}
						}
					}
					else if ((spy < 0) && (smy > 0))
					{
						T sss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (sss > 0)
						{
							if ((spz >= 0) && (smz >= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
							}
							else if ((spz <= 0) && (smz <= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
							}
							else if ((spz > 0) && (smz < 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
							}
							else if ((spz < 0) && (smz > 0))
							{
								T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
								if (ss > 0)
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_m_z)) - 1);
								}
								else
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y) + POW2(rho_p_z)) - 1);
								}
							}
						}
						else 
						{
							if ((spz >= 0) && (smz >= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
							}
							else if ((spz <= 0) && (smz <= 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
							}
							else if ((spz > 0) && (smz < 0))
							{
								rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x)) + POW2(rho_p_y) - 1);
							}
							else if ((spz < 0) && (smz > 0))
							{
								T ss = s*(abs(rho_p_z) - abs(rho_m_z))/(rho_p_z - rho_m_z);
								if (ss > 0)
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_m_z)) - 1);
								}
								else
								{
									rho_array(i, j, k) = rho_array(i, j, k) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y) + POW2(rho_p_z)) - 1);
								}
							}
						}
					}
				}
			}

			if (s*rho_array(i, j, k) < (T)0)
			{
				if (s < 0)
				{
					rho_array(i, j, k) = -epsilon;
				}
				if (s > 0)
				{
					rho_array(i, j, k) = epsilon;
				}
			}
		}
		multithreading.Sync(thread_id);
	}

	static void ENO3rd(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<TT>& velocity_x, const FIELD_STRUCTURE_3D<TT>& velocity_y, const FIELD_STRUCTURE_3D<TT>& velocity_z, const T& dt, MULTITHREADING& multithreading)
	{
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		const ARRAY_3D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_3D<TT>& velocity_array_y(velocity_y.array_for_this);
		const ARRAY_3D<TT>& velocity_array_z(velocity_z.array_for_this);
		
		T dx(rho.dx), dy(rho.dy), dz(rho.dz);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy), one_over_2dz(rho.one_over_2dz);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3), one_over_3dz(one_over_dx*(T)1/3);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			// x-component
			TT diff_1_x_n_3 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
			TT diff_1_x_n_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
			TT diff_1_x_n_1 = one_over_dx*(rho_ghost(i,     j, k) - rho_ghost(i - 1, j, k));
			TT diff_1_x_0   = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i,     j, k));
			TT diff_1_x_p_1 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
			TT diff_1_x_p_2 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));

			TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
			TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
			TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
			TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
			TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

			TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
			TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
			TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
			TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);
			
			TT rho_m_x;
			if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
			{
				if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
				}
				else
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
				}
				else
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}

			TT rho_p_x;
			if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;		
				}
				else
				{
					rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
				{
					rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
				}
				else
				{
					rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
				}
			}

			// y-component
			TT diff_1_y_n_3 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
			TT diff_1_y_n_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
			TT diff_1_y_n_1 = one_over_dy*(rho_ghost(i,     j, k) - rho_ghost(i, j - 1, k));
			TT diff_1_y_0   = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i,     j, k));
			TT diff_1_y_p_1 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
			TT diff_1_y_p_2 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));

			TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
			TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
			TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
			TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
			TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

			TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
			TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
			TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
			TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);
			
			TT rho_m_y;
			if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
			{
				if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
				}
				else
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
				}
				else
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}

			TT rho_p_y;
			if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;		
				}
				else
				{
					rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
				{
					rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
				}
				else
				{
					rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
				}
			}

			// z-component
			TT diff_1_z_n_3 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
			TT diff_1_z_n_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
			TT diff_1_z_n_1 = one_over_dz*(rho_ghost(i,     j, k) - rho_ghost(i, j, k - 1));
			TT diff_1_z_0   = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i,     j, k));
			TT diff_1_z_p_1 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
			TT diff_1_z_p_2 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));

			TT diff_2_z_n_2 = one_over_2dz*(diff_1_z_n_2 - diff_1_z_n_3);
			TT diff_2_z_n_1 = one_over_2dz*(diff_1_z_n_1 - diff_1_z_n_2);
			TT diff_2_z_0   = one_over_2dz*(diff_1_z_0   - diff_1_z_n_1);
			TT diff_2_z_p_1 = one_over_2dz*(diff_1_z_p_1 - diff_1_z_0  );
			TT diff_2_z_p_2 = one_over_2dz*(diff_1_z_p_2 - diff_1_z_p_1);

			TT diff_3_z_n_2 = one_over_3dz*(diff_2_z_n_1 - diff_2_z_n_2);
			TT diff_3_z_n_1 = one_over_3dz*(diff_2_z_0   - diff_2_z_n_1);
			TT diff_3_z_0   = one_over_3dz*(diff_2_z_p_1 - diff_2_z_0  );
			TT diff_3_z_p_1 = one_over_3dz*(diff_2_z_p_2 - diff_2_z_p_1);
			
			TT rho_m_z;
			if (abs(diff_2_z_n_1) <= abs(diff_2_z_0))
			{
				if (abs(diff_3_z_n_2) <= abs(diff_3_z_n_1))
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_2*dz*dz;
				}
				else
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_1*dz*dz;
				}
			}
			else
			{
				if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_n_1*dz*dz;
				}
				else
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_0*dz*dz;
				}
			}

			TT rho_p_z;
			if (abs(diff_2_z_0) <= abs(diff_2_z_p_1))
			{
				if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
				{
					rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_n_1*dz*dz;		
				}
				else
				{
					rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_0*dz*dz;
				}
			}
			else
			{
				if (abs(diff_3_z_0) <= abs(diff_3_z_p_1))
				{
					rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_0*dz*dz;
				}
				else
				{
					rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_p_1*dz*dz;
				}
			}

			T u_vel, v_vel, w_vel;
			// Velocity interpolation - Using MAC grid
			if (rho.is_scalar)
			{
				u_vel = (velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k))*(T)0.5;
				v_vel = (velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k))*(T)0.5;
				w_vel = (velocity_array_z(i, j, k + 1) + velocity_array_z(i, j, k))*(T)0.5;
			}
			else if (rho.is_x_component)
			{
				u_vel = velocity_array_x(i, j, k);
				v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i - 1, j, k) + velocity_array_y(i - 1, j + 1, k))*(T)0.25;
				w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i - 1, j, k) + velocity_array_z(i - 1, j, k + 1))*(T)0.25;
			}
			else if (rho.is_y_component)
			{
				u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j - 1, k) + velocity_array_x(i + 1, j - 1, k))*(T)0.25;
				v_vel = velocity_array_y(i, j, k);
				w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i, j - 1, k) + velocity_array_z(i, j - 1, k + 1))*(T)0.25;
			}
			else if (rho.is_z_component)
			{
				u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k - 1) + velocity_array_x(i + 1, j, k - 1))*(T)0.25;
				v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k - 1) + velocity_array_y(i, j + 1, k - 1))*(T)0.25;
				w_vel = velocity_array_z(i, j, k);
			}

			if (u_vel > 0)
			{
				if (v_vel > 0)
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
					}
				}
				else
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
					}
				}
			}
			else
			{
				if (v_vel > 0)
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
					}
				}
				else
				{
					if (w_vel > 0)
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
					}
					else
					{
						rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
					}
				}
			}
		}
		
		multithreading.Sync(thread_id);
	}

	static void ENO3rd(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, FIELD_STRUCTURE_3D<TT>& rho_ghost_x, FIELD_STRUCTURE_3D<TT>& rho_ghost_y, FIELD_STRUCTURE_3D<TT>& rho_ghost_z, const FIELD_STRUCTURE_3D<TT>& boundary_levelset, const FIELD_STRUCTURE_3D<TT>& velocity_x, const FIELD_STRUCTURE_3D<TT>& velocity_y, const FIELD_STRUCTURE_3D<TT>& velocity_z, const T& dt, MULTITHREADING& multithreading)
	{
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		const ARRAY_3D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_3D<TT>& velocity_array_y(velocity_y.array_for_this);
		const ARRAY_3D<TT>& velocity_array_z(velocity_z.array_for_this);
		
		T dx(rho.dx), dy(rho.dy), dz(rho.dz);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy), one_over_2dz(rho.one_over_2dz);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3), one_over_3dz(one_over_dz*(T)1/3);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			if (rho.fixed(i, j, k) == true)
			{
				continue;
			}

			if (boundary_levelset(i, j, k) < (T)0)
			{
				// x-component
				TT diff_1_x_n_3 = one_over_dx*(rho_ghost_x(i - 2, j, k) - rho_ghost_x(i - 3, j, k));
				TT diff_1_x_n_2 = one_over_dx*(rho_ghost_x(i - 1, j, k) - rho_ghost_x(i - 2, j, k));
				TT diff_1_x_n_1 = one_over_dx*(rho_ghost_x(i,     j, k) - rho_ghost_x(i - 1, j, k));
				TT diff_1_x_0   = one_over_dx*(rho_ghost_x(i + 1, j, k) - rho_ghost_x(i,     j, k));
				TT diff_1_x_p_1 = one_over_dx*(rho_ghost_x(i + 2, j, k) - rho_ghost_x(i + 1, j, k));
				TT diff_1_x_p_2 = one_over_dx*(rho_ghost_x(i + 3, j, k) - rho_ghost_x(i + 2, j, k));

				TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
				TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
				TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
				TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
				TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

				TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
				TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
				TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
				TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);

				TT rho_m_x;
				if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
				{
					if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}

				TT rho_p_x;
				if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;		
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
					}
				}

				// y-component
				TT diff_1_y_n_3 = one_over_dy*(rho_ghost_y(i, j - 2, k) - rho_ghost_y(i, j - 3, k));
				TT diff_1_y_n_2 = one_over_dy*(rho_ghost_y(i, j - 1, k) - rho_ghost_y(i, j - 2, k));
				TT diff_1_y_n_1 = one_over_dy*(rho_ghost_y(i,     j, k) - rho_ghost_y(i, j - 1, k));
				TT diff_1_y_0   = one_over_dy*(rho_ghost_y(i, j + 1, k) - rho_ghost_y(i,     j, k));
				TT diff_1_y_p_1 = one_over_dy*(rho_ghost_y(i, j + 2, k) - rho_ghost_y(i, j + 1, k));
				TT diff_1_y_p_2 = one_over_dy*(rho_ghost_y(i, j + 3, k) - rho_ghost_y(i, j + 2, k));

				TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
				TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
				TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
				TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
				TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

				TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
				TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
				TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
				TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);

				TT rho_m_y;
				if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
				{
					if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}

				TT rho_p_y;
				if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;		
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
					}
				}

				// z-component
				TT diff_1_z_n_3 = one_over_dz*(rho_ghost_z(i, j, k - 2) - rho_ghost_z(i, j, k - 3));
				TT diff_1_z_n_2 = one_over_dz*(rho_ghost_z(i, j, k - 1) - rho_ghost_z(i, j, k - 2));
				TT diff_1_z_n_1 = one_over_dz*(rho_ghost_z(i,     j, k) - rho_ghost_z(i, j, k - 1));
				TT diff_1_z_0   = one_over_dz*(rho_ghost_z(i, j, k + 1) - rho_ghost_z(i,     j, k));
				TT diff_1_z_p_1 = one_over_dz*(rho_ghost_z(i, j, k + 2) - rho_ghost_z(i, j, k + 1));
				TT diff_1_z_p_2 = one_over_dz*(rho_ghost_z(i, j, k + 3) - rho_ghost_z(i, j, k + 2));

				TT diff_2_z_n_2 = one_over_2dz*(diff_1_z_n_2 - diff_1_z_n_3);
				TT diff_2_z_n_1 = one_over_2dz*(diff_1_z_n_1 - diff_1_z_n_2);
				TT diff_2_z_0   = one_over_2dz*(diff_1_z_0   - diff_1_z_n_1);
				TT diff_2_z_p_1 = one_over_2dz*(diff_1_z_p_1 - diff_1_z_0  );
				TT diff_2_z_p_2 = one_over_2dz*(diff_1_z_p_2 - diff_1_z_p_1);

				TT diff_3_z_n_2 = one_over_3dz*(diff_2_z_n_1 - diff_2_z_n_2);
				TT diff_3_z_n_1 = one_over_3dz*(diff_2_z_0   - diff_2_z_n_1);
				TT diff_3_z_0   = one_over_3dz*(diff_2_z_p_1 - diff_2_z_0  );
				TT diff_3_z_p_1 = one_over_3dz*(diff_2_z_p_2 - diff_2_z_p_1);

				TT rho_m_z;
				if (abs(diff_2_z_n_1) <= abs(diff_2_z_0))
				{
					if (abs(diff_3_z_n_2) <= abs(diff_3_z_n_1))
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_2*dz*dz;
					}
					else
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_1*dz*dz;
					}
				}
				else
				{
					if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_n_1*dz*dz;
					}
					else
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_0*dz*dz;
					}
				}

				TT rho_p_z;
				if (abs(diff_2_z_0) <= abs(diff_2_z_p_1))
				{
					if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
					{
						rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_n_1*dz*dz;		
					}
					else
					{
						rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_0*dz*dz;
					}
				}
				else
				{
					if (abs(diff_3_z_0) <= abs(diff_3_z_p_1))
					{
						rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_0*dz*dz;
					}
					else
					{
						rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_p_1*dz*dz;
					}
				}

				T u_vel, v_vel, w_vel;
				// Velocity interpolation - Using MAC grid
				if (rho.is_scalar)
				{
					u_vel = (velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k))*(T)0.5;
					v_vel = (velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k))*(T)0.5;
					w_vel = (velocity_array_z(i, j, k + 1) + velocity_array_z(i, j, k))*(T)0.5;
				}
				else if (rho.is_x_component)
				{
					u_vel = velocity_array_x(i, j, k);
					v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i - 1, j, k) + velocity_array_y(i - 1, j + 1, k))*(T)0.25;
					w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i - 1, j, k) + velocity_array_z(i - 1, j, k + 1))*(T)0.25;
				}
				else if (rho.is_y_component)
				{
					u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j - 1, k) + velocity_array_x(i + 1, j - 1, k))*(T)0.25;
					v_vel = velocity_array_y(i, j, k);
					w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i, j - 1, k) + velocity_array_z(i, j - 1, k + 1))*(T)0.25;
				}
				else if (rho.is_z_component)
				{
					u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k - 1) + velocity_array_x(i + 1, j, k - 1))*(T)0.25;
					v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k - 1) + velocity_array_y(i, j + 1, k - 1))*(T)0.25;
					w_vel = velocity_array_z(i, j, k);
				}

				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				}
				else
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				} 
			}
			else
			{
				continue;
			}
		}
		
		multithreading.Sync(thread_id);
	}

	static void ENO3rd(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<TT>& boundary_levelset, const FIELD_STRUCTURE_3D<TT>& velocity_x, const FIELD_STRUCTURE_3D<TT>& velocity_y, const FIELD_STRUCTURE_3D<TT>& velocity_z, const T& dt, MULTITHREADING& multithreading)
	{
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		const ARRAY_3D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_3D<TT>& velocity_array_y(velocity_y.array_for_this);
		const ARRAY_3D<TT>& velocity_array_z(velocity_z.array_for_this);
		
		T dx(rho.dx), dy(rho.dy), dz(rho.dz);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy), one_over_2dz(rho.one_over_2dz);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3), one_over_3dz(one_over_dx*(T)1/3);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			if (rho.fixed(i, j, k) == true)
			{
				continue;
			}
			if (boundary_levelset(i, j, k) <= (T)0)
			{
				// x-component
				TT diff_1_x_n_3 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
				TT diff_1_x_n_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
				TT diff_1_x_n_1 = one_over_dx*(rho_ghost(i,     j, k) - rho_ghost(i - 1, j, k));
				TT diff_1_x_0   = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i,     j, k));
				TT diff_1_x_p_1 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
				TT diff_1_x_p_2 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));

				TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
				TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
				TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
				TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
				TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

				TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
				TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
				TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
				TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);

				TT rho_m_x;
				if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
				{
					if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}

				TT rho_p_x;
				if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;		
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
					}
				}

				// y-component
				TT diff_1_y_n_3 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
				TT diff_1_y_n_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
				TT diff_1_y_n_1 = one_over_dy*(rho_ghost(i,     j, k) - rho_ghost(i, j - 1, k));
				TT diff_1_y_0   = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i,     j, k));
				TT diff_1_y_p_1 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
				TT diff_1_y_p_2 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));

				TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
				TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
				TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
				TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
				TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

				TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
				TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
				TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
				TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);

				TT rho_m_y;
				if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
				{
					if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}

				TT rho_p_y;
				if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;		
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
					}
				}

				// z-component
				TT diff_1_z_n_3 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
				TT diff_1_z_n_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
				TT diff_1_z_n_1 = one_over_dz*(rho_ghost(i,     j, k) - rho_ghost(i, j, k - 1));
				TT diff_1_z_0   = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i,     j, k));
				TT diff_1_z_p_1 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
				TT diff_1_z_p_2 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));

				TT diff_2_z_n_2 = one_over_2dz*(diff_1_z_n_2 - diff_1_z_n_3);
				TT diff_2_z_n_1 = one_over_2dz*(diff_1_z_n_1 - diff_1_z_n_2);
				TT diff_2_z_0   = one_over_2dz*(diff_1_z_0   - diff_1_z_n_1);
				TT diff_2_z_p_1 = one_over_2dz*(diff_1_z_p_1 - diff_1_z_0  );
				TT diff_2_z_p_2 = one_over_2dz*(diff_1_z_p_2 - diff_1_z_p_1);

				TT diff_3_z_n_2 = one_over_3dz*(diff_2_z_n_1 - diff_2_z_n_2);
				TT diff_3_z_n_1 = one_over_3dz*(diff_2_z_0   - diff_2_z_n_1);
				TT diff_3_z_0   = one_over_3dz*(diff_2_z_p_1 - diff_2_z_0  );
				TT diff_3_z_p_1 = one_over_3dz*(diff_2_z_p_2 - diff_2_z_p_1);

				TT rho_m_z;
				if (abs(diff_2_z_n_1) <= abs(diff_2_z_0))
				{
					if (abs(diff_3_z_n_2) <= abs(diff_3_z_n_1))
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_2*dz*dz;
					}
					else
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_1*dz*dz;
					}
				}
				else
				{
					if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_n_1*dz*dz;
					}
					else
					{
						rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_0*dz*dz;
					}
				}

				TT rho_p_z;
				if (abs(diff_2_z_0) <= abs(diff_2_z_p_1))
				{
					if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
					{
						rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_n_1*dz*dz;		
					}
					else
					{
						rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_0*dz*dz;
					}
				}
				else
				{
					if (abs(diff_3_z_0) <= abs(diff_3_z_p_1))
					{
						rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_0*dz*dz;
					}
					else
					{
						rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_p_1*dz*dz;
					}
				}

				T u_vel, v_vel, w_vel;
				// Velocity interpolation - Using MAC grid
				if (rho.is_scalar)
				{
					u_vel = (velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k))*(T)0.5;
					v_vel = (velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k))*(T)0.5;
					w_vel = (velocity_array_z(i, j, k + 1) + velocity_array_z(i, j, k))*(T)0.5;
				}
				else if (rho.is_x_component)
				{
					u_vel = velocity_array_x(i, j, k);
					v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i - 1, j, k) + velocity_array_y(i - 1, j + 1, k))*(T)0.25;
					w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i - 1, j, k) + velocity_array_z(i - 1, j, k + 1))*(T)0.25;
				}
				else if (rho.is_y_component)
				{
					u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j - 1, k) + velocity_array_x(i + 1, j - 1, k))*(T)0.25;
					v_vel = velocity_array_y(i, j, k);
					w_vel = (velocity_array_z(i, j, k) + velocity_array_z(i, j, k + 1) + velocity_array_z(i, j - 1, k) + velocity_array_z(i, j - 1, k + 1))*(T)0.25;
				}
				else if (rho.is_z_component)
				{
					u_vel = (velocity_array_x(i, j, k) + velocity_array_x(i + 1, j, k) + velocity_array_x(i, j, k - 1) + velocity_array_x(i + 1, j, k - 1))*(T)0.25;
					v_vel = (velocity_array_y(i, j, k) + velocity_array_y(i, j + 1, k) + velocity_array_y(i, j, k - 1) + velocity_array_y(i, j + 1, k - 1))*(T)0.25;
					w_vel = velocity_array_z(i, j, k);
				}

				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				}
				else
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				} 
			}
			else
			{
				continue;
			}
		}
		
		multithreading.Sync(thread_id);
	}

	static void ENO3rd(const int& thread_id, FIELD_STRUCTURE_3D<TT>& rho, FIELD_STRUCTURE_3D<TT>& rho_ghost, const FIELD_STRUCTURE_3D<VT>& velocity, const T& dt, MULTITHREADING& multithreading, const T& epsilon, const bool& use_rk2 = false)
	{
		const ARRAY_3D<VT>& velocity_array(velocity.array_for_this);
		ARRAY_3D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(thread_id, rho_array, true);

		multithreading.Sync(thread_id);
		
		T dx(rho.dx), dy(rho.dy), dz(rho.dz);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_dz(rho.one_over_dz), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy), one_over_2dz(rho.one_over_2dz);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3), one_over_3dz(one_over_dx*(T)1/3);
		
		GRID_ITERATION_3D(rho.partial_grids[thread_id])
		{
			// x-component
			TT diff_1_x_n_3 = one_over_dx*(rho_ghost(i - 2, j, k) - rho_ghost(i - 3, j, k));
			TT diff_1_x_n_2 = one_over_dx*(rho_ghost(i - 1, j, k) - rho_ghost(i - 2, j, k));
			TT diff_1_x_n_1 = one_over_dx*(rho_ghost(i,     j, k) - rho_ghost(i - 1, j, k));
			TT diff_1_x_0   = one_over_dx*(rho_ghost(i + 1, j, k) - rho_ghost(i,     j, k));
			TT diff_1_x_p_1 = one_over_dx*(rho_ghost(i + 2, j, k) - rho_ghost(i + 1, j, k));
			TT diff_1_x_p_2 = one_over_dx*(rho_ghost(i + 3, j, k) - rho_ghost(i + 2, j, k));

			TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
			TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
			TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
			TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
			TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

			TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
			TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
			TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
			TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);
			
			TT rho_m_x;
			if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
			{
				if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
				}
				else
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
				}
				else
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}

			TT rho_p_x;
			if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;		
				}
				else
				{
					rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
				{
					rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
				}
				else
				{
					rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
				}
			}

			// y-component
			TT diff_1_y_n_3 = one_over_dy*(rho_ghost(i, j - 2, k) - rho_ghost(i, j - 3, k));
			TT diff_1_y_n_2 = one_over_dy*(rho_ghost(i, j - 1, k) - rho_ghost(i, j - 2, k));
			TT diff_1_y_n_1 = one_over_dy*(rho_ghost(i,     j, k) - rho_ghost(i, j - 1, k));
			TT diff_1_y_0   = one_over_dy*(rho_ghost(i, j + 1, k) - rho_ghost(i,     j, k));
			TT diff_1_y_p_1 = one_over_dy*(rho_ghost(i, j + 2, k) - rho_ghost(i, j + 1, k));
			TT diff_1_y_p_2 = one_over_dy*(rho_ghost(i, j + 3, k) - rho_ghost(i, j + 2, k));

			TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
			TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
			TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
			TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
			TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

			TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
			TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
			TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
			TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);
			
			TT rho_m_y;
			if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
			{
				if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
				}
				else
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
				}
				else
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}

			TT rho_p_y;
			if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;		
				}
				else
				{
					rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
				{
					rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
				}
				else
				{
					rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
				}
			}

			// z-component
			TT diff_1_z_n_3 = one_over_dz*(rho_ghost(i, j, k - 2) - rho_ghost(i, j, k - 3));
			TT diff_1_z_n_2 = one_over_dz*(rho_ghost(i, j, k - 1) - rho_ghost(i, j, k - 2));
			TT diff_1_z_n_1 = one_over_dz*(rho_ghost(i,     j, k) - rho_ghost(i, j, k - 1));
			TT diff_1_z_0   = one_over_dz*(rho_ghost(i, j, k + 1) - rho_ghost(i,     j, k));
			TT diff_1_z_p_1 = one_over_dz*(rho_ghost(i, j, k + 2) - rho_ghost(i, j, k + 1));
			TT diff_1_z_p_2 = one_over_dz*(rho_ghost(i, j, k + 3) - rho_ghost(i, j, k + 2));

			TT diff_2_z_n_2 = one_over_2dz*(diff_1_z_n_2 - diff_1_z_n_3);
			TT diff_2_z_n_1 = one_over_2dz*(diff_1_z_n_1 - diff_1_z_n_2);
			TT diff_2_z_0   = one_over_2dz*(diff_1_z_0   - diff_1_z_n_1);
			TT diff_2_z_p_1 = one_over_2dz*(diff_1_z_p_1 - diff_1_z_0  );
			TT diff_2_z_p_2 = one_over_2dz*(diff_1_z_p_2 - diff_1_z_p_1);

			TT diff_3_z_n_2 = one_over_3dz*(diff_2_z_n_1 - diff_2_z_n_2);
			TT diff_3_z_n_1 = one_over_3dz*(diff_2_z_0   - diff_2_z_n_1);
			TT diff_3_z_0   = one_over_3dz*(diff_2_z_p_1 - diff_2_z_0  );
			TT diff_3_z_p_1 = one_over_3dz*(diff_2_z_p_2 - diff_2_z_p_1);
			
			TT rho_m_z;
			if (abs(diff_2_z_n_1) <= abs(diff_2_z_0))
			{
				if (abs(diff_3_z_n_2) <= abs(diff_3_z_n_1))
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_2*dz*dz;
				}
				else
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_n_1*dz + 2*diff_3_z_n_1*dz*dz;
				}
			}
			else
			{
				if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_n_1*dz*dz;
				}
				else
				{
					rho_m_z = diff_1_z_n_1 + diff_2_z_0*dz - diff_3_z_0*dz*dz;
				}
			}

			TT rho_p_z;
			if (abs(diff_2_z_0) <= abs(diff_2_z_p_1))
			{
				if (abs(diff_3_z_n_1) <= abs(diff_3_z_0))
				{
					rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_n_1*dz*dz;		
				}
				else
				{
					rho_p_z = diff_1_z_0 - diff_2_z_0*dz - diff_3_z_0*dz*dz;
				}
			}
			else
			{
				if (abs(diff_3_z_0) <= abs(diff_3_z_p_1))
				{
					rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_0*dz*dz;
				}
				else
				{
					rho_p_z = diff_1_z_0 - diff_2_z_p_1*dz + 2*diff_3_z_p_1*dz*dz;
				}
			}

			// Velocity interpolation - Using MAC grid
			T u_vel = (velocity_array(i - 1, j, k).x + velocity_array(i + 1, j, k).x)*(T)0.5;
			T v_vel = (velocity_array(i, j - 1, k).y + velocity_array(i, j + 1, k).y)*(T)0.5;
			T w_vel = (velocity_array(i, j, k - 1).z + velocity_array(i, j, k + 1).z)*(T)0.5;

			if (!use_rk2)
			{
				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_m_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				}
				else
				{
					if (v_vel > 0)
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_m_y + w_vel*rho_p_z);
						}
					}
					else
					{
						if (w_vel > 0)
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_m_z);
						}
						else
						{
							rho_array(i, j, k) = rho_array(i, j, k) - dt*(u_vel*rho_p_x + v_vel*rho_p_y + w_vel*rho_p_z);
						}
					}
				}
			}
		}
		multithreading.Sync(thread_id);
	}
};
