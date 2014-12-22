#pragma once

#include "LINEAR_SOLVER.h"
#include "LEVELSET_3D.h"
#include "CG_METHOD.h"
#include "PCG_METHOD.h"
#include "GAUSS_SEIDEL_METHOD.h"

class POISSON_SOLVER_3D
{
public: // Essential Data
	T									tolerance;
	T									sqr_tolerance;
	int									max_iteration;
	int									num_smoother_id;
	static int							smoother_id;				// Number of gss call per each frame, typically num_levels - 1

	bool								use_variable_density;
	bool								place_dirichlet_bc_at_face;	// Default is true as it prevents multithreading artifact under multilevel Poisson solver
	static bool							use_detailed_log;

	ARRAY<DYNAMIC_ARRAY<int>>			nb_bix_lists;				// Narrow band cell indices, for each thread

	GRID_STRUCTURE_3D*					grid_ghost;					// Note : To calculate bix of pressure field previously, this grid should have same size
	int									ghost_width;				// Pressure ghost cell width
	ARRAY<GRID_STRUCTURE_3D>*			partial_grids;				// For each thread

public: // Temporary, for statistics
	int									num_levels;					// Number of levels for mg, ml
	// For call of smoother
	static T*							euclidean_res_accums;
	static T*							infinity_res_accums;
	static int*							num_frames_ids;				// This is the number of frames, and it's current frame + 1

public: // For matrix-vector type solvers
	CSR_MATRIX<T>						A;
	VECTOR_ND<T>						x;
	VECTOR_ND<T>						b;

	// For Sherman-Morrison Formula
	CSR_MATRIX<T>						B;
	VECTOR_ND<T>						y_1;
	VECTOR_ND<T>						y_2;
	VECTOR_ND<T>						w;
	VECTOR_ND<T>						z;
	T									beta_u, beta_d, beta;

public: // pre-computed matrix-free Gauss-Seidel multiplication coefficients for direct (non matrix-vector) type solvers
	FIELD_STRUCTURE_3D<GS_COEFFICIENTS> gs_coeffients;

public: // Multithreading
	MULTITHREADING*						multithreading;

public: // Solver
	LINEAR_SOLVER*						linear_solver;

public: // Properties of Simulation
	bool								air_water_simulation, oil_water_simulation;

public: // Boundary Condition
	bool								Dirichlet_Boundary_Condition, Neumann_Boundary_Condition;

public: // Constructors and Destructor
	POISSON_SOLVER_3D(void)
		: tolerance((T)1e-4), sqr_tolerance(tolerance*tolerance), max_iteration(100), multithreading(0), linear_solver(0), num_levels(0), use_variable_density(false), place_dirichlet_bc_at_face(true), num_smoother_id(0), air_water_simulation(false), oil_water_simulation(false)
	{}

	~POISSON_SOLVER_3D(void)
	{
		DELETE_POINTER(linear_solver);
		DELETE_ARRAY(euclidean_res_accums);
		DELETE_ARRAY(infinity_res_accums);
		DELETE_ARRAY(num_frames_ids);
	}

public: // Initialization Functions
	void Initialize(MULTITHREADING* multithreading_input, const T& tolerance_input, const int& max_itr_input, GRID_STRUCTURE_3D* grid_ghost_input = 0, const int ghost_width_input = 0,  ARRAY<GRID_STRUCTURE_3D>* partial_grids_input = 0, int num_levels_input = 1, bool use_nb_gss_input = false)
	{
		num_levels = num_levels_input;
		ghost_width = ghost_width_input;

		// Initializing residual concerned members
		if(num_levels_input > 1)
		{
			// Reserving enough spaces for statistics
			num_smoother_id = num_levels_input*3;
			
			DELETE_ARRAY(euclidean_res_accums);
			DELETE_ARRAY(infinity_res_accums);
			DELETE_ARRAY(num_frames_ids);

			euclidean_res_accums = new T[num_smoother_id];
			infinity_res_accums  = new T[num_smoother_id];
			num_frames_ids		 = new int[num_smoother_id];

			for (int i = 0; i < num_smoother_id; i++)
			{
				euclidean_res_accums[i] = (T)0;
				infinity_res_accums[i]  = (T)0;
				num_frames_ids[i]		= 0;
			}
		}
		
		if(use_nb_gss_input)
		{
			grid_ghost = grid_ghost_input;
			partial_grids = partial_grids_input;

			// Initializing narrow band cell index list
			int num_thread = multithreading_input->num_threads;
			int max_bix = (int)((T)grid_ghost_input->ijk_res/(T)num_thread);

			nb_bix_lists.Initialize(num_thread);
			for (int thread = 0; thread < num_thread; thread++)
			{
				nb_bix_lists[thread].Initialize(max_bix,1);
			}
		}
		
		multithreading = multithreading_input;
		max_iteration = max_itr_input;

		InitializeLinearSolver(CG);			// Defalut solver is PCG

		SetTolerance(tolerance_input);
 
	}

	void InitializeLinearSolver(const POISSON_SOLVER_TYPE linear_solver_type)
	{
		DELETE_POINTER(linear_solver);

		switch(linear_solver_type)
		{
		case CG:
			linear_solver = new CG_METHOD();
			break;
		case PCG:
			linear_solver = new PCG_METHOD();
			break;
		case GS:
			linear_solver = new GAUSS_SEIDEL_METHOD();
			break;
		case NO_SOLVER:
		default:
			linear_solver = 0;
			break;
		}

		linear_solver->Initialize(multithreading, tolerance, max_iteration);
	}

	void InitializePressure(const int& thread_id, FIELD_STRUCTURE_3D<T>& pressure)
	{
		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			pressure(i, j, k) = (T)0;
		}
		END_GRID_ITERATION_3D;
	}

public: // Member Functions
	void Solve(const int& thread_id, FIELD_STRUCTURE_3D<T>& pressure, FIELD_STRUCTURE_3D<T>& density, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div, const FIELD_STRUCTURE_3D<T>& variable, const LEVELSET_3D& levelset,  const LEVELSET_3D& boundary_levelset, const FIELD_STRUCTURE_3D<T>& jc_on_solution, const FIELD_STRUCTURE_3D<T>& jc_on_derivative)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve");

		BuildLinearSystemNodeJumpConditionVaribleCoefficient(thread_id, A, x, b, pressure, variable, bc, div, levelset, boundary_levelset, jc_on_solution, jc_on_derivative);
		
		T sum_for_b(0);
		
		if (Neumann_Boundary_Condition)
		{
			BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
			{
				if (bc(i, j, k) > -1)
				{
					sum_for_b += b[bc(i, j, k)];
				}
			}
			END_GRID_ITERATION_SUM(sum_for_b);

			T one_over_b_dimension = 1/(T)b.num_dimension;
			
			BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
			{
				if (bc(i, j, k) > -1)
				{
					b[bc(i, j, k)] -= one_over_b_dimension*sum_for_b;	
				}
			}
			END_GRID_ITERATION_3D;
		}

		linear_solver->Solve(thread_id, A, x, b, bc);
		VectorToGrid(thread_id, x, pressure, bc);

		LOG::End(thread_id);
	}

	void SolveForViscosity(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve For Semi-Implicit Method");

		BuildLinearSystemNodeForSemiImplicitViscosity(thread_id, A, x, b, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6, velocity, bc, explicit_term);
		
		linear_solver->Solve(thread_id, A, x, b, bc);
		VectorToGrid(thread_id, x, velocity, bc);

		LOG::End(thread_id);
	}

	void SolveForViscosity(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity, FIELD_STRUCTURE_3D<T>& velocity_ghost, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term, bool is_vertical = true)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve For Semi-Implicit Method");

		BuildLinearSystemNodeForSemiImplicitViscosity(thread_id, A, x, b, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6, velocity_ghost, bc, explicit_term);
		
		linear_solver->Solve(thread_id, A, x, b, bc);
		VectorToGrid(thread_id, x, velocity, bc);
		
		// Need to fix this when you simulate vertical pipe
		if (is_vertical)
		{
			if (velocity.is_y_component)
			{
				BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
				{
					velocity(i, velocity.grid.j_end, k) = velocity(i, velocity.grid.j_start + 1, k);
				}
				END_GRID_ITERATION_3D;
			}  
		}
		else
		{
			if (velocity.is_x_component)
			{
				BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
				{
					velocity(velocity.grid.i_end, j, k) = velocity(velocity.grid.i_start + 1, j, k);
				}
				END_GRID_ITERATION_3D;
			}  
		}
		LOG::End(thread_id);
	}

	void SolveForViscosity(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity, FIELD_STRUCTURE_3D<T>& velocity_ghost, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term, const FIELD_STRUCTURE_3D<T>& boundary_levelset, bool is_vertical = true)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve For Semi-Implicit Method");

		BuildLinearSystemNodeForSemiImplicitViscosity(thread_id, A, x, b, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6, velocity_ghost, bc, explicit_term, boundary_levelset);
		
		linear_solver->Solve(thread_id, A, x, b, bc);
		VectorToGrid(thread_id, x, velocity, bc);
		
		// Need to fix this when you simulate vertical pipe
		if (is_vertical)
		{
			if (velocity.is_y_component)
			{
				BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
				{
					velocity(i, velocity.grid.j_end, k) = velocity(i, velocity.grid.j_start + 1, k);
				}
				END_GRID_ITERATION_3D;
			}  
		}
		else
		{
			if (velocity.is_x_component)
			{
				BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
				{
					velocity(velocity.grid.i_end, j, k) = velocity(velocity.grid.i_start + 1, j, k);
				}
				END_GRID_ITERATION_3D;
			}  
		}
		LOG::End(thread_id);
	}

	void SolveForViscosity(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity, FIELD_STRUCTURE_3D<T>& velocity_ghost, FIELD_STRUCTURE_3D<T>& velocity_ghost_x, FIELD_STRUCTURE_3D<T>& velocity_ghost_y, FIELD_STRUCTURE_3D<T>& velocity_ghost_z, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term, bool is_vertical = true)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve For Semi-Implicit Method");

		BuildLinearSystemNodeForSemiImplicitViscosity(thread_id, A, x, b, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6, velocity_ghost, velocity_ghost_x, velocity_ghost_y, velocity_ghost_z, bc, explicit_term);

		linear_solver->Solve(thread_id, A, x, b, bc);
		VectorToGrid(thread_id, x, velocity, bc);
		
		// Need to fix this when you simulate vertical pipe
		if (is_vertical)
		{
			if (velocity.is_y_component)
			{
				BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
				{
					velocity(i, velocity.grid.j_end, k) = velocity(i, velocity.grid.j_start + 1, k);
				}
				END_GRID_ITERATION_3D;
			}  
		}
		else
		{
			if (velocity.is_x_component)
			{
				BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
				{
					velocity(velocity.grid.i_end, j, k) = velocity(velocity.grid.i_start + 1, j, k);
				}
				END_GRID_ITERATION_3D;
			}  
		}

		LOG::End(thread_id);
	}

	void SolveForViscosityYPeriodic(const int& thread_id, FIELD_STRUCTURE_3D<T>& velocity, FIELD_STRUCTURE_3D<T>& velocity_ghost, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term, bool is_vertical = true)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve For Semi-Implicit Method");

		BuildLinearSystemNodeForSemiImplicitViscosityPeriodic(thread_id, A, y_1, b, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6, velocity_ghost, bc, explicit_term);

		linear_solver->Solve(thread_id, A, y_1, b, bc);
		
		BEGIN_HEAD_THREAD_WORK
		{
			w.Initialize(b.num_dimension, true);
			z.Initialize(b.num_dimension, true);
			x.Initialize(b.num_dimension, true);
		}
		END_HEAD_THREAD_WORK;

		BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}
			if (bc(i, j - 1, k) < 0)
			{
				w[bc(i, j, k)] = -1;
			}
			else if (bc(i, j + 1, k) < 0)
			{
				w[bc(i, j, k)] = -1;
			}
			else
			{
				w[bc(i, j, k)] = 0;
			}
		}
		END_GRID_ITERATION_3D;

		BuildLinearSystemNodeForSemiImplicitViscosityPeriodic(thread_id, B, y_2, w, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6, velocity_ghost, bc);
		
		linear_solver->Solve(thread_id, B, y_2, w, bc);

		BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}
			if (bc(i, j - 1, k) < 0)
			{
				z[bc(i, bc.grid.j_end - (j + 1 - bc.grid.j_start), k)] = -coef_4(i, j, k)*bc.one_over_dy2;
			}
			else if (bc(i, j + 1, k) < 0)
			{
				z[bc(i, bc.grid.j_start + (j + 1 - bc.grid.j_end), k)] = -coef_3(i, j, k)*bc.one_over_dy2;
			}
			else
			{
				z[bc(i, j, k)] = 0;
			}
		}
		END_GRID_ITERATION_3D;

		BEGIN_HEAD_THREAD_WORK
		{
			beta_d = 0;
			beta_u = 0;
			beta   = 0;
		}
		END_HEAD_THREAD_WORK;

		DotProduct(multithreading, thread_id, z, y_1, beta_u);
		DotProduct(multithreading, thread_id, z, y_2, beta_d);

		beta = beta_u/((T)1 - beta_d);

		const int start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

		for (int i = start_ix; i <= end_ix; i++)
		{
			x[i] = y_1[i] + beta*y_2[i];
		}
		multithreading->Sync(thread_id);

		VectorToGrid(thread_id, x, velocity, bc);

		//// Need to fix this when you simulate vertical pipe
		//if (is_vertical)
		//{
		//	/*if (velocity.is_x_component)
		//	{
		//		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		//		{
		//			velocity(i, velocity.grid.j_end, k) = velocity(i, velocity.grid.j_start, k);
		//		}
		//		END_GRID_ITERATION_3D;
		//	}*/
		//	if (velocity.is_y_component)
		//	{
		//		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		//		{
		//			velocity(i, velocity.grid.j_end, k) = velocity(i, velocity.grid.j_start + 1, k);
		//		}
		//		END_GRID_ITERATION_3D;
		//	}  
		//}
		//else
		//{
		//	if (velocity.is_x_component)
		//	{
		//		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		//		{
		//			velocity(velocity.grid.i_end, j, k) = velocity(velocity.grid.i_start + 1, j, k);
		//		}
		//		END_GRID_ITERATION_3D;
		//	}  
		//}

		LOG::End(thread_id);
	}

	void Solve2ndOrder(const int& thread_id, FIELD_STRUCTURE_3D<T>& pressure, FIELD_STRUCTURE_3D<T>& density, FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div, const LEVELSET_3D& levelset)
	{
		assert(linear_solver != 0);

		LOG::Begin(thread_id, "Poisson Solve");

		BuildLinearSystem2ndOrderFreesurface(thread_id, A, x, b, pressure, density, bc, div, levelset);
		linear_solver->Solve(thread_id, A, x, b, bc);
		VectorToGrid(thread_id, x, pressure, bc);

		LOG::End(thread_id);
	}

	void BuildLinearSystem(const int& thread_id, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<T>& density, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div)
	{
		if (place_dirichlet_bc_at_face)
		{
			BuildLinearSystemFaceDirichlet(thread_id, A, x, b, pressure, bc, div);
		}
		else
		{
			BuildLinearSystemNodeDirichlet(thread_id, A, x, b, pressure, bc, div);
		}
	}

	int AssignSequentialindicesToFullCells(const int& thread_id, const FIELD_STRUCTURE_3D<int>& bc)
	{
		// Count number of full cells in this thread domain
		int full_ix(0);

		BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
		{
			if(bc(i, j, k) > -1) ++full_ix;
		}
		END_GRID_ITERATION_3D;

		int start_full_ix, end_full_ix;
		multithreading->SyncDomainIndices1D(thread_id, full_ix, start_full_ix, end_full_ix);

		BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
		{
			if(bc(i, j, k) > -1) bc(i, j, k) = start_full_ix++;
		}
		END_GRID_ITERATION_3D;

		assert(start_full_ix - 1 == end_full_ix);

		multithreading->SyncSum(thread_id, full_ix);

		return full_ix;
	}

	int CountNonZeroElements(const int& thread_id, const FIELD_STRUCTURE_3D<int>& bc)
	{
		int nnz(0);

		BEGIN_GRID_ITERATION_3D(bc.partial_grids[thread_id])
		{
			if(bc(i, j, k) > -1)
			{
				nnz++;
				if(bc(i+1, j, k) > -1) nnz++;
				if(bc(i-1, j, k) > -1) nnz++;
				if(bc(i, j+1, k) > -1) nnz++;
				if(bc(i, j-1, k) > -1) nnz++;
				if(bc(i, j, k+1) > -1) nnz++;
				if(bc(i, j, k-1) > -1) nnz++;
				// For Periodic Boundary Condition
				if(bc(i+1, j, k) == BC_PER) nnz++;
				if(bc(i-1, j, k) == BC_PER) nnz++;
				if(bc(i, j+1, k) == BC_PER) nnz++;
				if(bc(i, j-1, k) == BC_PER) nnz++;
				if(bc(i, j, k+1) == BC_PER) nnz++;
				if(bc(i, j, k-1) == BC_PER) nnz++;
			}
		}
		END_GRID_ITERATION_SUM(nnz);

		return nnz;
	}

	void BuildLinearSystemNodeDirichlet(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);

		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Make A matrix, x vector, and b vector
		const T dxdx = POW2(pressure.grid.dx), inv_dxdx = (T)1/dxdx;

		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			if(bc(i, j, k) < 0) continue;

			int coef_ijk = 0;		// For optimization, inv_dxdx is multiplied at the end
			
			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -inv_dxdx);
			}
			if(bc(i - 1, j, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -inv_dxdx);
			}
			if(bc(i, j + 1, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -inv_dxdx);
			}
			if(bc(i, j - 1, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -inv_dxdx);
			}
			if(bc(i, j, k + 1) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -inv_dxdx);
			}
			if(bc(i, j, k - 1) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -inv_dxdx);
			}

			if(bc(i + 1, j , k) == BC_DIR) coef_ijk += 1;
			if(bc(i - 1, j , k) == BC_DIR) coef_ijk += 1;
			if(bc(i, j + 1 , k) == BC_DIR) coef_ijk += 1;
			if(bc(i, j - 1 , k) == BC_DIR) coef_ijk += 1;
			if(bc(i, j , k + 1) == BC_DIR) coef_ijk += 1;
			if(bc(i, j , k - 1) == BC_DIR) coef_ijk += 1;

			if(coef_ijk == 0) coef_ijk = 1;				// No nearby FULL and DIR cells, all neumann condition (null space)

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), inv_dxdx*(T)coef_ijk);

			b_vector[bc(i, j, k)] = -div(i, j, k);		// Note: We solver -Ax = -b;
			x_vector[bc(i, j, k)] = pressure(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void BuildLinearSystemFaceDirichlet(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK
		
		// Make A matrix, x vector, and b vector
		const T dxdx = POW2(pressure.grid.dx), inv_dxdx = (T)1/dxdx;

		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			if(bc(i, j, k) < 0) continue;

			int coef_ijk = 0;		// For optimization, inv_dxdx is multiplied at the end
			
			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -inv_dxdx);
			}
			if(bc(i - 1, j, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -inv_dxdx);
			}
			if(bc(i, j + 1, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -inv_dxdx);
			}
			if(bc(i, j - 1, k) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -inv_dxdx);
			}
			if(bc(i, j, k + 1) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -inv_dxdx);
			}
			if(bc(i, j, k - 1) > -1)
			{
				coef_ijk += 1;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -inv_dxdx);
			}

			if(bc(i + 1, j , k) == BC_DIR) coef_ijk += 2;
			if(bc(i - 1, j , k) == BC_DIR) coef_ijk += 2;
			if(bc(i, j + 1 , k) == BC_DIR) coef_ijk += 2;
			if(bc(i, j - 1 , k) == BC_DIR) coef_ijk += 2;
			if(bc(i, j , k + 1) == BC_DIR) coef_ijk += 2;
			if(bc(i, j , k - 1) == BC_DIR) coef_ijk += 2;

			if(coef_ijk == 0) coef_ijk = 1;				// No nearby FULL and DIR cells, all neumann condition (null space)

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), (T)coef_ijk*inv_dxdx);

			b_vector[bc(i, j, k)] = -div(i, j, k);		// Note: We solver -Ax = -b;
			x_vector[bc(i, j, k)] = pressure(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void BuildLinearSystemVariableDensity(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<T>& prejection_density, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);

		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Make A matrix, x vector, and b vector
		const T dxdx = POW2(pressure.grid.dx), inv_dxdx = (T)1/dxdx;

		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			if(bc(i, j, k) < 0) continue;

			const T density_ijk = prejection_density(i, j, k);

			T coef_ijk = 0;
			
			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				const T coef = (T)2/(density_ijk + prejection_density(i + 1, j, k));
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef);
			}
			if(bc(i - 1, j, k) > -1)
			{
				const T coef = (T)2/(density_ijk + prejection_density(i - 1, j, k));
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef);
			}
			if(bc(i, j + 1, k) > -1)
			{
				const T coef = (T)2/(density_ijk + prejection_density(i, j + 1, k));
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef);
			}
			if(bc(i, j - 1, k) > -1)
			{
				const T coef = (T)2/(density_ijk + prejection_density(i, j - 1, k));
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef);
			}
			if(bc(i, j, k + 1) > -1)
			{
				const T coef = (T)2/(density_ijk + prejection_density(i , j, k + 1));
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef);
			}
			if(bc(i, j, k - 1) > -1)
			{
				const T coef = (T)2/(density_ijk + prejection_density(i , j, k - 1));
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef);
			}

			if(bc(i + 1, j , k) == BC_DIR) coef_ijk += (T)2/(density_ijk + prejection_density(i + 1, j, k));
			if(bc(i - 1, j , k) == BC_DIR) coef_ijk += (T)2/(density_ijk + prejection_density(i - 1, j, k));
			if(bc(i, j + 1 , k) == BC_DIR) coef_ijk += (T)2/(density_ijk + prejection_density(i, j + 1, k));
			if(bc(i, j - 1 , k) == BC_DIR) coef_ijk += (T)2/(density_ijk + prejection_density(i, j - 1, k));
			if(bc(i, j , k + 1) == BC_DIR) coef_ijk += (T)2/(density_ijk + prejection_density(i, j, k + 1));
			if(bc(i, j , k - 1) == BC_DIR) coef_ijk += (T)2/(density_ijk + prejection_density(i, j, k - 1));

			if(coef_ijk == (T)0) coef_ijk = (T)1;				// No nearby FULL and DIR cells, all neumann condition (null space)

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), coef_ijk);

			b_vector[bc(i, j, k)] = -div(i, j, k)*dxdx;		// Note: We solver -Ax = -b;
			x_vector[bc(i, j, k)] = pressure(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	// After studying boundary condition set up, review this one more time
	void BuildLinearSystem2ndOrderFreesurface(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<T>& density, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div, const LEVELSET_3D& levelset)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);

		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.start_ix[thread_id - 1] + multithreading->sync_value_int[thread_id - 1];
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = A_matrix.prev_row_array[thread_id - 1] + multithreading->sync_value_int_array[thread_id - 1];
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Make A matrix, x vector, and b vector
		const T dxdx = POW2(pressure.grid.dx);

		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			if(bc(i, j, k) < 0) continue;

			T coef_ijk = 0;
			
			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				const T coef = (T)1;
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef);
			}
			if(bc(i - 1, j, k) > -1)
			{
				const T coef = (T)1;
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef);
			}
			if(bc(i, j + 1, k) > -1)
			{
				const T coef = (T)1;
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef);
			}
			if(bc(i, j - 1, k) > -1)
			{
				const T coef = (T)1;
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef);
			}
			if(bc(i, j, k + 1) > -1)
			{
				const T coef = (T)1;
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef);
			}
			if(bc(i, j, k - 1) > -1)
			{
				const T coef = (T)1;
				coef_ijk += coef;
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef);
			}

			const T phi_ijk = MAX(ABS(levelset.phi(i, j, k)), (T)1e-4);

			if(bc(i + 1, j , k) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i + 1, j, k)), (T)1e-4);
				coef_ijk += phi/phi_ijk + (T)1;
			}
			if(bc(i - 1, j , k) == BC_DIR) 
			{
				const T phi = MAX(ABS(levelset.phi(i - 1, j, k)), (T)1e-4);
				coef_ijk += phi/phi_ijk + (T)1;
			}
			if(bc(i, j + 1 , k) == BC_DIR) 
			{
				const T phi = MAX(ABS(levelset.phi(i, j + 1, k)), (T)1e-4);
				coef_ijk += phi/phi_ijk + (T)1;
			}
			if(bc(i, j - 1 , k) == BC_DIR) 
			{
				const T phi = MAX(ABS(levelset.phi(i, j - 1, k)), (T)1e-4);
				coef_ijk += phi/phi_ijk + (T)1;
			}
			if(bc(i, j , k + 1) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i, j, k + 1)), (T)1e-4);
				coef_ijk += phi/phi_ijk + (T)1;
			}
			if(bc(i, j , k - 1) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i, j, k - 1)), (T)1e-4);
				coef_ijk += phi/phi_ijk + (T)1;
			}

			if(coef_ijk == (T)0) coef_ijk = (T)1;				// No nearby FULL and DIR cells, all neumann condition (null space)

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), coef_ijk);

			b_vector[bc(i, j, k)] = -div(i, j, k)*dxdx;		// Note: We solver -Ax = -b;
			x_vector[bc(i, j, k)] = pressure(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}
	
	void BuildLinearSystemNodeJumpConditionVaribleCoefficient(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<T>& variable, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div, const LEVELSET_3D& levelset, const LEVELSET_3D& boundary_levelset, const FIELD_STRUCTURE_3D<T>& jc_on_solution, const FIELD_STRUCTURE_3D<T>& jc_on_derivative)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speed up variables 
		int i_start(pressure.i_start), i_end(pressure.i_end), j_start(pressure.j_start), j_end(pressure.j_end);
		T dx(pressure.dx), dy(pressure.dy), one_over_dx(pressure.one_over_dx), one_over_dy(pressure.one_over_dy), one_over_dz(pressure.one_over_dz), one_over_dx2(pressure.one_over_dx2), one_over_dy2(pressure.one_over_dy2), one_over_dz2(pressure.one_over_dz2);

		const T dxdx = POW2(pressure.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(pressure.grid.dy), inv_dydy = (T)1/dydy, dzdz = POW2(pressure.grid.dz), inv_dzdz = (T)1/dzdz;
		
		int i_start_for_domain(0), j_start_for_domain(0), k_start_for_domain(0);

		if (air_water_simulation)
		{
			BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
			{
				if (bc(i, j, k) < 0)
				{
					continue;
				}

				T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
				// Define betas
				T beta_l, beta_r, beta_b, beta_t, beta_d, beta_u;
				T levelset_c = levelset(i, j, k), levelset_l = levelset(i - 1, j, k), levelset_r = levelset(i + 1, j, k), levelset_b = levelset(i, j - 1, k), levelset_t = levelset(i, j + 1, k), levelset_d = levelset(i, j, k - 1), levelset_u = levelset(i, j, k + 1);
	
				beta_l = variable(i, j, k)*variable(i - 1, j, k)*(abs(levelset_l) + abs(levelset_c))/(variable(i, j, k)*abs(levelset_l) + variable(i - 1, j, k)*abs(levelset_c));
				beta_r = variable(i, j, k)*variable(i + 1, j, k)*(abs(levelset_r) + abs(levelset_c))/(variable(i + 1, j, k)*abs(levelset_c) + variable(i, j, k)*abs(levelset_r));
				beta_b = variable(i, j, k)*variable(i, j - 1, k)*(abs(levelset_b) + abs(levelset_c))/(variable(i, j, k)*abs(levelset_b) + variable(i, j - 1, k)*abs(levelset_c));
				beta_t = variable(i, j, k)*variable(i, j + 1, k)*(abs(levelset_t) + abs(levelset_c))/(variable(i, j + 1, k)*abs(levelset_c) + variable(i, j, k)*abs(levelset_t));
				beta_d = variable(i, j, k)*variable(i, j, k - 1)*(abs(levelset_d) + abs(levelset_c))/(variable(i, j, k)*abs(levelset_d) + variable(i, j, k - 1)*abs(levelset_c));
				beta_u = variable(i, j, k)*variable(i, j, k + 1)*(abs(levelset_u) + abs(levelset_c))/(variable(i, j, k + 1)*abs(levelset_c) + variable(i, j, k)*abs(levelset_u));

				// If neighbor is full cell
				if (bc(i - 1, j, k) > -1)
				{
					coef_ijk += beta_l;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -beta_l*inv_dxdx);
				}
				if (bc(i + 1, j, k) > -1)
				{
					coef_ijk += beta_r;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -beta_r*inv_dxdx);
				}
				if (bc(i, j - 1, k) > -1)
				{
					coef_ijk += beta_b;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -beta_b*inv_dydy);
				}
				if (bc(i, j + 1, k) > -1)
				{
					coef_ijk += beta_t;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -beta_t*inv_dydy);
				}
				if (bc(i, j, k - 1) > -1)
				{
					coef_ijk += beta_d;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -beta_d*inv_dzdz);
				}
				if (bc(i, j, k + 1) > -1)
				{
					coef_ijk += beta_u;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -beta_u*inv_dzdz);
				}
				
				if (bc(i - 1, j, k) == BC_DIR)
				{
					coef_ijk += beta_l;
				}
				if (bc(i + 1, j, k) == BC_DIR)
				{
					coef_ijk += beta_r;
				}
				if (bc(i, j - 1, k) == BC_DIR)
				{
					coef_ijk += beta_b;
				}
				if (bc(i, j + 1, k) == BC_DIR)
				{
					coef_ijk += beta_t;
				}
				if (bc(i, j, k - 1) == BC_DIR)
				{
					coef_ijk += beta_d;
				}
				if (bc(i, j, k + 1) == BC_DIR)
				{
					coef_ijk += beta_u;
				}
				
				// For saving starting index of the given domain
				if (bc(i, j, k) == 0)
				{
					i_start_for_domain = i;
					j_start_for_domain = j;
					k_start_for_domain = k;
				}
				
				// Tuning the given matrix into the invertiable matrix -- "On deflation and singular symmetric positive semi-definite matrices -- J.M.Tang, C.Vuik"
				if ((thread_id == 0) && (i == i_start_for_domain) && (j == j_start_for_domain) && (k == k_start_for_domain))
				{
					if (bc(i - 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i + 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j - 1, k) == BC_NEUM)
					{
						coef_ijk += beta_b;
					}
					if (bc(i, j + 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j, k - 1) == BC_NEUM)
					{
						coef_ijk += beta_d;
					}
					if (bc(i, j, k + 1) == BC_NEUM)
					{
						coef_ijk += beta_u;
					}
				}
				else
				{
					if (bc(i - 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i + 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j - 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j + 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j, k - 1) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j, k + 1) == BC_NEUM)
					{
						coef_ijk += 0;
					}
				}
							
				if (coef_ijk == 0)
				{
					coef_ijk = 1;
				}

				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), inv_dxdx*(T)coef_ijk);
				
				// Divide given region into different region - Boundary Capturing method for Poisson Equation on irregular domain
				T F_L, F_R, F_B ,F_T, F_D, F_U;
					
				// Left arm Stencil
				T subcell_l = abs(levelset(i - 1, j, k))/(abs(levelset(i, j, k)) + abs(levelset(i - 1, j, k)));
				T a_l = (jc_on_solution(i, j, k)*abs(levelset(i - 1, j, k)) + jc_on_solution(i - 1, j, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i - 1, j, k)));
				T b_l = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).x*abs(levelset(i - 1, j, k)) + jc_on_derivative(i - 1, j, k)*levelset.normal(i - 1, j, k).x*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i - 1, j, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i - 1, j, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i - 1, j, k) > 0)))
				{
					F_L = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i - 1, j, k) > 0))
				{
					F_L = one_over_dx2*a_l*beta_l - one_over_dx*beta_l*b_l*subcell_l/variable(i - 1, j, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i - 1, j, k) <= 0))
				{
					F_L = -one_over_dx2*a_l*beta_l + one_over_dx*beta_l*b_l*subcell_l/variable(i - 1, j, k);
				}
	
				// Right arm Stencil
				T subcell_r = abs(levelset(i + 1, j, k))/(abs(levelset(i, j, k)) + abs(levelset(i + 1, j, k)));
				T a_r = (jc_on_solution(i, j, k)*abs(levelset(i + 1, j, k)) + jc_on_solution(i + 1, j, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i + 1, j, k)));
				T b_r = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).x*abs(levelset(i + 1, j, k)) + jc_on_derivative(i + 1, j, k)*levelset.normal(i + 1, j, k).x*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i + 1, j, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i + 1, j, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i + 1, j, k) > 0)))
				{
					F_R = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i + 1, j, k) > 0))
				{
					F_R = one_over_dx2*a_r*beta_r + one_over_dx*beta_r*b_r*subcell_r/variable(i + 1, j, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i + 1, j, k) <= 0))
				{
					F_R = -one_over_dx2*a_r*beta_r - one_over_dx*beta_r*b_r*subcell_r/variable(i + 1, j, k);
				}
	
				// Bottom arm Stencil
				T subcell_b = abs(levelset(i, j - 1, k))/(abs(levelset(i, j, k)) + abs(levelset(i, j - 1, k)));
				T a_b = (jc_on_solution(i, j, k)*abs(levelset(i, j - 1, k)) + jc_on_solution(i, j - 1, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j - 1, k)));
				T b_b = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).y*abs(levelset(i, j - 1, k)) + jc_on_derivative(i, j - 1, k)*levelset.normal(i, j - 1, k).y*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j - 1, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i, j - 1, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j - 1, k) > 0)))
				{
					F_B = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j - 1, k) > 0))
				{
					F_B = one_over_dy2*a_b*beta_b - one_over_dy*beta_b*b_b*subcell_b/variable(i, j - 1, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j - 1, k) <= 0))
				{
					F_B = -one_over_dy2*a_b*beta_b + one_over_dy*beta_b*b_b*subcell_b/variable(i, j - 1, k);
				}
	
				// Top arm Stencil
				T subcell_t = abs(levelset(i, j + 1, k))/(abs(levelset(i, j, k)) + abs(levelset(i, j + 1, k)));
				T a_t = (jc_on_solution(i, j, k)*abs(levelset(i, j + 1, k)) + jc_on_solution(i, j + 1, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j + 1, k)));
				T b_t = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).y*abs(levelset(i, j + 1, k)) + jc_on_derivative(i, j + 1, k)*levelset.normal(i, j + 1, k).y*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j + 1, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i, j + 1, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j + 1, k) > 0)))
				{
					F_T = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j + 1, k) > 0))
				{
					F_T = one_over_dy2*a_t*beta_t + one_over_dy*beta_t*b_t*subcell_t/variable(i, j + 1, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j + 1, k) <= 0))
				{
					F_T = -one_over_dy2*a_t*beta_t - one_over_dy*beta_t*b_t*subcell_t/variable(i, j + 1, k);
				}
				
				// Downward arm Stencil
				T subcell_d = abs(levelset(i, j, k - 1))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k - 1)));
				T a_d = (jc_on_solution(i, j, k)*abs(levelset(i, j, k - 1)) + jc_on_solution(i, j, k - 1)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k - 1)));
				T b_d = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).z*abs(levelset(i, j, k - 1)) + jc_on_derivative(i, j, k - 1)*levelset.normal(i, j, k - 1).z*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k - 1)));
	
				if (((levelset(i, j, k) <= 0) && (levelset(i, j, k - 1) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j, k - 1) > 0)))
				{
					F_D = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j, k - 1) > 0))
				{
					F_D = one_over_dz2*a_d*beta_d - one_over_dz*beta_d*b_d*subcell_d/variable(i, j, k - 1);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j, k - 1) <= 0))
				{
					F_D = -one_over_dz2*a_d*beta_d + one_over_dz*beta_d*b_d*subcell_d/variable(i, j, k - 1);
				}
	
				// Upward arm Stencil
				T subcell_u = abs(levelset(i, j, k + 1))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k + 1)));
				T a_u = (jc_on_solution(i, j, k)*abs(levelset(i, j, k + 1)) + jc_on_solution(i, j, k + 1)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k + 1)));
				T b_u = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).z*abs(levelset(i, j, k + 1)) + jc_on_derivative(i, j, k + 1)*levelset.normal(i, j, k + 1).z*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k + 1)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i, j, k + 1) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j, k + 1) > 0)))
				{
					F_U = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j, k + 1) > 0))
				{
					F_U = one_over_dz2*a_u*beta_u + one_over_dz*beta_u*b_u*subcell_u/variable(i, j, k + 1);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j, k + 1) <= 0))
				{
					F_U = -one_over_dz2*a_u*beta_u - one_over_dz*beta_u*b_u*subcell_u/variable(i, j, k + 1);
				}
	
				T F_X = F_L + F_R, F_Y = F_B + F_T, F_Z = F_D + F_U;
	
				b_vector[bc(i, j, k)] = -div(i, j, k) - F_X - F_Y - F_Z;
				
				if (bc(i - 1, j, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_l*inv_dxdx*pressure(i - 1, j, k);
				}
				if (bc(i + 1, j, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_r*inv_dxdx*pressure(i + 1, j, k);
				}
				if (bc(i, j - 1, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_b*inv_dydy*pressure(i, j - 1, k);
				}
				if (bc(i, j + 1, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_t*inv_dydy*pressure(i, j + 1, k);
				}
				if (bc(i, j, k - 1) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_d*inv_dzdz*pressure(i, j, k - 1);
				}		
				if (bc(i, j, k + 1) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_u*inv_dzdz*pressure(i, j, k + 1);
				}		
	
				x_vector[bc(i, j, k)] = pressure(i, j, k);
			}
			END_GRID_ITERATION_3D;
		}

		if (oil_water_simulation)
		{
			BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
			{
				if (bc(i, j, k) < 0)
				{
					continue;
				}

				T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
				// Define betas
				T beta_l, beta_r, beta_b, beta_t, beta_d, beta_u;
				T levelset_c = levelset(i, j, k), levelset_l = levelset(i - 1, j, k), levelset_r = levelset(i + 1, j, k), levelset_b = levelset(i, j - 1, k), levelset_t = levelset(i, j + 1, k), levelset_d = levelset(i, j, k - 1), levelset_u = levelset(i, j, k + 1);
	
				beta_l = variable(i, j, k)*variable(i - 1, j, k)*(abs(levelset_l) + abs(levelset_c))/(variable(i, j, k)*abs(levelset_l) + variable(i - 1, j, k)*abs(levelset_c));
				beta_r = variable(i, j, k)*variable(i + 1, j, k)*(abs(levelset_r) + abs(levelset_c))/(variable(i + 1, j, k)*abs(levelset_c) + variable(i, j, k)*abs(levelset_r));
				beta_b = variable(i, j, k)*variable(i, j - 1, k)*(abs(levelset_b) + abs(levelset_c))/(variable(i, j, k)*abs(levelset_b) + variable(i, j - 1, k)*abs(levelset_c));
				beta_t = variable(i, j, k)*variable(i, j + 1, k)*(abs(levelset_t) + abs(levelset_c))/(variable(i, j + 1, k)*abs(levelset_c) + variable(i, j, k)*abs(levelset_t));
				beta_d = variable(i, j, k)*variable(i, j, k - 1)*(abs(levelset_d) + abs(levelset_c))/(variable(i, j, k)*abs(levelset_d) + variable(i, j, k - 1)*abs(levelset_c));
				beta_u = variable(i, j, k)*variable(i, j, k + 1)*(abs(levelset_u) + abs(levelset_c))/(variable(i, j, k + 1)*abs(levelset_c) + variable(i, j, k)*abs(levelset_u));

				// If neighbor is full cell
				if (bc(i - 1, j, k) > -1)
				{
					coef_ijk += beta_l;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -beta_l*inv_dxdx);
				}
				if (bc(i + 1, j, k) > -1)
				{
					coef_ijk += beta_r;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -beta_r*inv_dxdx);
				}
				if (bc(i, j - 1, k) > -1)
				{
					coef_ijk += beta_b;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -beta_b*inv_dydy);
				}
				if (bc(i, j + 1, k) > -1)
				{
					coef_ijk += beta_t;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -beta_t*inv_dydy);
				}
				if (bc(i, j, k - 1) > -1)
				{
					coef_ijk += beta_d;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -beta_d*inv_dzdz);
				}
				if (bc(i, j, k + 1) > -1)
				{
					coef_ijk += beta_u;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -beta_u*inv_dzdz);
				}
				
				// Periodic Boundary Condition
				if (bc(i, j - 1, k) == BC_PER)
				{
					coef_ijk += beta_b;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j_end + (j - 1 - j_start), k), -beta_b*inv_dydy);
				}
				if (bc(i, j + 1, k) == BC_PER)
				{
					coef_ijk += beta_t;
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j_start + (j + 1 - j_end), k), -beta_t*inv_dydy);
				}

				if (bc(i - 1, j, k) == BC_DIR)
				{
					coef_ijk += beta_l;
				}
				if (bc(i + 1, j, k) == BC_DIR)
				{
					coef_ijk += beta_r;
				}
				if (bc(i, j - 1, k) == BC_DIR)
				{
					coef_ijk += beta_b;
				}
				if (bc(i, j + 1, k) == BC_DIR)
				{
					coef_ijk += beta_t;
				}
				if (bc(i, j, k - 1) == BC_DIR)
				{
					coef_ijk += beta_d;
				}
				if (bc(i, j, k + 1) == BC_DIR)
				{
					coef_ijk += beta_u;
				}
				
				// For saving starting index of the given domain
				if (bc(i, j, k) == 0)
				{
					i_start_for_domain = i;
					j_start_for_domain = j;
					k_start_for_domain = k;
				}
				
				// Tuning the given matrix into the invertiable matrix -- "On deflation and singular symmetric positive semi-definite matrices -- J.M.Tang, C.Vuik"
				if ((thread_id == 0) && (i == i_start_for_domain) && (j == j_start_for_domain) && (k == k_start_for_domain))
				{
					if (bc(i - 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i + 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j - 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j + 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j, k - 1) == BC_NEUM)
					{
						coef_ijk += beta_d;
						//coef_ijk += 0;
					}
					if (bc(i, j, k + 1) == BC_NEUM)
					{
						coef_ijk += 0;
					}
				}
				else
				{
					if (bc(i - 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i + 1, j, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j - 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j + 1, k) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j, k - 1) == BC_NEUM)
					{
						coef_ijk += 0;
					}
					if (bc(i, j, k + 1) == BC_NEUM)
					{
						coef_ijk += 0;
					}
				}

				// For biCG
				//if (bc(i - 1, j, k) == BC_NEUM)
				//{
				//	T theta = acos(-boundary_levelset.normal(i - 1, j, k).x);
				//	T tan_theta = tan(theta);
				//	
				//	// For the floating point error
				//	T epsilon = 10e-7;
				//	
				//	if (tan_theta > 0 && tan_theta <= 1)
				//	{
				//		if (boundary_levelset.normal(i - 1, j, k).z > 0)
				//		{
				//			coef_ijk += beta_l*tan(theta);
				//			A_matrix(bc(i, j, k), bc(i, j, k - 1)) -= beta_l*inv_dxdx*tan(theta);
				//		}
				//		else
				//		{
				//			coef_ijk += beta_l*tan(theta);
				//			A_matrix(bc(i, j, k), bc(i, j, k + 1)) -= beta_l*inv_dxdx*tan(theta);
				//		}
				//	}
				//	else if (tan_theta > 1 + epsilon)
				//	{
				//		if (boundary_levelset.normal(i - 1, j, k).z > 0)
				//		{
				//			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k - 1), beta_r*inv_dxdx*(1 - 1/tan(theta)));
				//			A_matrix(bc(i, j, k), bc(i, j, k - 1)) -= beta_r*inv_dxdx*1/tan(theta);
				//		}
				//		else
				//		{
				//			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k - 1), beta_r*inv_dxdx*(1 - 1/tan(theta)));
				//			A_matrix(bc(i, j, k), bc(i, j, k + 1)) -= beta_r*inv_dxdx*1/tan(theta);
				//		}
				//	}
				//	else
				//	{
				//		coef_ijk += 0;
				//	}
				//}
				//
				//if (bc(i + 1, j, k) == BC_NEUM)
				//{
				//	T theta = acos(boundary_levelset.normal(i + 1, j, k).x);
				//	T tan_theta = tan(theta);

				//	// For the floating point error
				//	T epsilon = 10e-7;

				//	if (tan_theta > 0 && tan_theta <= 1)
				//	{
				//		if (boundary_levelset.normal(i + 1, j, k).z > (T)0)
				//		{
				//			coef_ijk += beta_r*inv_dxdx*tan(theta);
				//			A_matrix(bc(i, j, k), bc(i, j, k - 1)) -= beta_r*inv_dxdx*tan(theta);
				//		}
				//		else
				//		{
				//			coef_ijk += beta_r*inv_dxdx*tan(theta);
				//			A_matrix(bc(i, j, k), bc(i, j, k + 1)) -= beta_r*inv_dxdx*tan(theta);
				//		}
				//	}
				//	else if (tan_theta > 1 + epsilon)
				//	{
				//		if (boundary_levelset.normal(i + 1, j, k).z > (T)0)
				//		{
				//			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k - 1), -beta_r*inv_dxdx*(1 - 1/tan(theta)));
				//			A_matrix(bc(i, j, k), bc(i, j, k - 1)) -= beta_r*inv_dxdx*1/tan(theta);
				//		}
				//		else
				//		{
				//			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k + 1), -beta_r*inv_dxdx*(1 - 1/tan(theta)));
				//			A_matrix(bc(i, j, k), bc(i, j, k + 1)) -= beta_r*inv_dxdx*1/tan(theta);
				//		}
				//	}
				//	else
				//	{
				//		coef_ijk += 0;
				//	}
				//}
				//if (bc(i, j - 1, k) == BC_NEUM)
				//{
				//	coef_ijk += 0;
				//}
				//if (bc(i, j + 1, k) == BC_NEUM)
				//{
				//	coef_ijk += 0;
				//}
				//if (bc(i, j, k - 1) == BC_NEUM)
				//{
				//	T theta = acos(-boundary_levelset.normal(i, j, k - 1).x);
				//	T epsilon = 10e-7;

				//	if (theta > 0 && theta <= PI/4)
				//	{
				//		A_matrix(bc(i, j, k), bc(i + 1, j, k)) -= beta_b*inv_dzdz*(tan(theta));
				//		A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k - 1), -beta_b*inv_dzdz*(1 - tan(theta)));
				//	}
				//	else if (theta > PI/4 + epsilon && theta < PI/2)
				//	{
				//		coef_ijk += beta_b*inv_dzdz*1/tan(theta);
				//		A_matrix(bc(i, j, k), bc(i + 1, j, k)) -= beta_b*inv_dzdz*(1 - 1/tan(theta));
				//	}
				//	else if (theta > PI/2 + epsilon && theta < 3*PI/4)
				//	{
				//		coef_ijk += beta_b*inv_dzdz*1/tan(PI-theta);
				//		A_matrix(bc(i, j, k), bc(i - 1, j, k)) -= beta_b*inv_dzdz*(1 - 1/tan(PI-theta));
				//	}
				//	else if (theta >= 3*PI/4 && theta < PI)
				//	{
				//		A_matrix(bc(i, j, k), bc(i - 1, j, k)) -= beta_b*inv_dzdz*(tan(PI-theta));
				//		A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k - 1), -beta_b*inv_dzdz*(1 - tan(PI-theta)));
				//	}
				//	else
				//	{
				//		coef_ijk += 0;
				//	}
				//}
				//if (bc(i, j, k + 1) == BC_NEUM)
				//{
				//	T theta = acos(boundary_levelset.normal(i, j, k + 1).x);	
				//	
				//	// For floting point error
				//	T epsilon = 10e-7;

				//	if (theta > 0 && theta <= PI/4)
				//	{
				//		A_matrix(bc(i, j, k), bc(i - 1, j, k)) -= beta_u*inv_dzdz*tan(theta);
				//		A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k - 1), -beta_u*inv_dzdz*(1 - tan(theta)));
				//	}	
				//	else if (theta > PI/4 + epsilon && theta < PI/2)
				//	{
				//		coef_ijk += beta_u*inv_dzdz*1/tan(theta);
				//		A_matrix(bc(i, j, k), bc(i - 1, j, k)) -= beta_u*inv_dzdz*(1/tan(theta));
				//	}
				//	else if (theta > PI/2 + epsilon && theta < 3*PI/4)
				//	{
				//		coef_ijk += beta_u*inv_dzdz*1/tan(PI-theta);
				//		A_matrix(bc(i, j, k), bc(i + 1, j, k)) -= beta_u*inv_dzdz*(1/tan(PI-theta));
				//	}
				//	else if (theta >= 3*PI/4 && theta < PI)
				//	{
				//		A_matrix(bc(i, j, k), bc(i + 1, j, k)) -= beta_u*inv_dzdz*tan(PI-theta);
				//		A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k - 1), -beta_u*inv_dzdz*(1 - tan(PI-theta)));
				//	}
				//	else
				//	{
				//		coef_ijk += 0;
				//	}					
				//}
				
				if (coef_ijk == 0)
				{
					coef_ijk = 1;
				}
				
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), inv_dxdx*(T)coef_ijk);
				
				// Divide given region into different region - Boundary Capturing method for Poisson Equation on irregular domain
				T F_L, F_R, F_B ,F_T, F_D, F_U;
					
				// Left arm Stencil
				T subcell_l = abs(levelset(i - 1, j, k))/(abs(levelset(i, j, k)) + abs(levelset(i - 1, j, k)));
				T a_l = (jc_on_solution(i, j, k)*abs(levelset(i - 1, j, k)) + jc_on_solution(i - 1, j, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i - 1, j, k)));
				T b_l = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).x*abs(levelset(i - 1, j, k)) + jc_on_derivative(i - 1, j, k)*levelset.normal(i - 1, j, k).x*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i - 1, j, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i - 1, j, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i - 1, j, k) > 0)))
				{
					F_L = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i - 1, j, k) > 0))
				{
					F_L = one_over_dx2*a_l*beta_l - one_over_dx*beta_l*b_l*subcell_l/variable(i - 1, j, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i - 1, j, k) <= 0))
				{
					F_L = -one_over_dx2*a_l*beta_l + one_over_dx*beta_l*b_l*subcell_l/variable(i - 1, j, k);
				}
	
				// Right arm Stencil
				T subcell_r = abs(levelset(i + 1, j, k))/(abs(levelset(i, j, k)) + abs(levelset(i + 1, j, k)));
				T a_r = (jc_on_solution(i, j, k)*abs(levelset(i + 1, j, k)) + jc_on_solution(i + 1, j, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i + 1, j, k)));
				T b_r = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).x*abs(levelset(i + 1, j, k)) + jc_on_derivative(i + 1, j, k)*levelset.normal(i + 1, j, k).x*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i + 1, j, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i + 1, j, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i + 1, j, k) > 0)))
				{
					F_R = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i + 1, j, k) > 0))
				{
					F_R = one_over_dx2*a_r*beta_r + one_over_dx*beta_r*b_r*subcell_r/variable(i + 1, j, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i + 1, j, k) <= 0))
				{
					F_R = -one_over_dx2*a_r*beta_r - one_over_dx*beta_r*b_r*subcell_r/variable(i + 1, j, k);
				}
	
				// Bottom arm Stencil
				T subcell_b = abs(levelset(i, j - 1, k))/(abs(levelset(i, j, k)) + abs(levelset(i, j - 1, k)));
				T a_b = (jc_on_solution(i, j, k)*abs(levelset(i, j - 1, k)) + jc_on_solution(i, j - 1, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j - 1, k)));
				T b_b = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).y*abs(levelset(i, j - 1, k)) + jc_on_derivative(i, j - 1, k)*levelset.normal(i, j - 1, k).y*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j - 1, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i, j - 1, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j - 1, k) > 0)))
				{
					F_B = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j - 1, k) > 0))
				{
					F_B = one_over_dy2*a_b*beta_b - one_over_dy*beta_b*b_b*subcell_b/variable(i, j - 1, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j - 1, k) <= 0))
				{
					F_B = -one_over_dy2*a_b*beta_b + one_over_dy*beta_b*b_b*subcell_b/variable(i, j - 1, k);
				}
	
				// Top arm Stencil
				T subcell_t = abs(levelset(i, j + 1, k))/(abs(levelset(i, j, k)) + abs(levelset(i, j + 1, k)));
				T a_t = (jc_on_solution(i, j, k)*abs(levelset(i, j + 1, k)) + jc_on_solution(i, j + 1, k)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j + 1, k)));
				T b_t = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).y*abs(levelset(i, j + 1, k)) + jc_on_derivative(i, j + 1, k)*levelset.normal(i, j + 1, k).y*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j + 1, k)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i, j + 1, k) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j + 1, k) > 0)))
				{
					F_T = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j + 1, k) > 0))
				{
					F_T = one_over_dy2*a_t*beta_t + one_over_dy*beta_t*b_t*subcell_t/variable(i, j + 1, k);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j + 1, k) <= 0))
				{
					F_T = -one_over_dy2*a_t*beta_t - one_over_dy*beta_t*b_t*subcell_t/variable(i, j + 1, k);
				}
				
				// Downward arm Stencil
				T subcell_d = abs(levelset(i, j, k - 1))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k - 1)));
				T a_d = (jc_on_solution(i, j, k)*abs(levelset(i, j, k - 1)) + jc_on_solution(i, j, k - 1)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k - 1)));
				T b_d = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).z*abs(levelset(i, j, k - 1)) + jc_on_derivative(i, j, k - 1)*levelset.normal(i, j, k - 1).z*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k - 1)));
	
				if (((levelset(i, j, k) <= 0) && (levelset(i, j, k - 1) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j, k - 1) > 0)))
				{
					F_D = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j, k - 1) > 0))
				{
					F_D = one_over_dz2*a_d*beta_d - one_over_dz*beta_d*b_d*subcell_d/variable(i, j, k - 1);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j, k - 1) <= 0))
				{
					F_D = -one_over_dz2*a_d*beta_d + one_over_dz*beta_d*b_d*subcell_d/variable(i, j, k - 1);
				}
	
				// Upward arm Stencil
				T subcell_u = abs(levelset(i, j, k + 1))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k + 1)));
				T a_u = (jc_on_solution(i, j, k)*abs(levelset(i, j, k + 1)) + jc_on_solution(i, j, k + 1)*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k + 1)));
				T b_u = (jc_on_derivative(i, j, k)*levelset.normal(i, j, k).z*abs(levelset(i, j, k + 1)) + jc_on_derivative(i, j, k + 1)*levelset.normal(i, j, k + 1).z*abs(levelset(i, j, k)))/(abs(levelset(i, j, k)) + abs(levelset(i, j, k + 1)));
				
				if (((levelset(i, j, k) <= 0) && (levelset(i, j, k + 1) <= 0)) || ((levelset(i, j, k) > 0) && (levelset(i, j, k + 1) > 0)))
				{
					F_U = 0;
				}
				else if ((levelset(i, j, k) <= 0) && (levelset(i, j, k + 1) > 0))
				{
					F_U = one_over_dz2*a_u*beta_u + one_over_dz*beta_u*b_u*subcell_u/variable(i, j, k + 1);
				}
				else if ((levelset(i, j, k) > 0) && (levelset(i, j, k + 1) <= 0))
				{
					F_U = -one_over_dz2*a_u*beta_u - one_over_dz*beta_u*b_u*subcell_u/variable(i, j, k + 1);
				}
	
				T F_X = F_L + F_R, F_Y = F_B + F_T, F_Z = F_D + F_U;
	
				b_vector[bc(i, j, k)] = -div(i, j, k) - F_X - F_Y - F_Z;
				
				if (bc(i - 1, j, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_l*inv_dxdx*pressure(i - 1, j, k);
				}
				if (bc(i + 1, j, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_r*inv_dxdx*pressure(i + 1, j, k);
				}
				if (bc(i, j - 1, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_b*inv_dydy*pressure(i, j - 1, k);
				}
				if (bc(i, j + 1, k) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_t*inv_dydy*pressure(i, j + 1, k);
				}
				if (bc(i, j, k - 1) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_d*inv_dzdz*pressure(i, j, k - 1);
				}		
				if (bc(i, j, k + 1) == BC_DIR)
				{
					b_vector[bc(i, j, k)] += beta_u*inv_dzdz*pressure(i, j, k + 1);
				}		
	
				x_vector[bc(i, j, k)] = pressure(i, j, k);
			}
			END_GRID_ITERATION_3D;
		}
	}
	
	void BuildLinearSystemNodeForSemiImplicitViscosity(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, const FIELD_STRUCTURE_3D<T>& velocity, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speed up variables 
		int i_res(velocity.grid.i_res), j_res(velocity.grid.j_res), k_res(velocity.grid.k_res), i_start(velocity.i_start), i_end(velocity.i_end), j_start(velocity.j_start), j_end(velocity.j_end), k_start(velocity.k_start), k_end(velocity.k_end);
		T dx(velocity.dx), dy(velocity.dy), one_over_dx(velocity.one_over_dx), one_over_dy(velocity.one_over_dy), one_over_dz(velocity.one_over_dz), one_over_dx2(velocity.one_over_dx2), one_over_dy2(velocity.one_over_dy2), one_over_dz2(velocity.one_over_dz2);
		T x_min(velocity.grid.x_min), y_min(velocity.grid.y_min), z_min(velocity.grid.z_min), x_max(velocity.grid.x_max), y_max(velocity.grid.y_max), z_max(velocity.grid.z_max);

		const T dxdx = POW2(velocity.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(velocity.grid.dy), inv_dydy = (T)1/dydy, dzdz = POW2(velocity.grid.dz), inv_dzdz = (T)1/dzdz;
		
		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j, k) > -1)
			{
				coef_ijk += coef_2(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef_2(i, j, k)*inv_dxdx);
			}
			if (bc(i + 1, j, k) > -1)
			{
				coef_ijk += coef_1(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef_1(i, j, k)*inv_dxdx);
			}
			if (bc(i, j - 1, k) > -1)
			{
				coef_ijk += coef_4(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef_4(i, j, k)*inv_dydy);
			}
			if (bc(i, j + 1, k) > -1)
			{
				coef_ijk += coef_3(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef_3(i, j, k)*inv_dydy);
			}
			if (bc(i, j, k - 1) > -1)
			{
				coef_ijk += coef_6(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef_6(i, j, k)*inv_dzdz);
			}
			if (bc(i, j, k + 1) > -1)
			{
				coef_ijk += coef_5(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef_5(i, j, k)*inv_dzdz);
			}

			if (bc(i - 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_2(i, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_1(i, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += coef_4(i, j, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += coef_3(i, j, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				coef_ijk += coef_6(i, j, k);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				coef_ijk += coef_5(i, j, k);
			}
			
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), (T)1 + inv_dxdx*(T)coef_ijk);
			
			b_vector[bc(i, j, k)] = explicit_term(i, j, k);

			if (bc(i - 1, j, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_2(i, j, k)*inv_dxdx*velocity(i - 1, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_1(i, j, k)*inv_dxdx*velocity(i + 1, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_4(i, j, k)*inv_dydy*velocity(i, j - 1, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_3(i, j, k)*inv_dydy*velocity(i, j + 1, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_6(i, j, k)*inv_dzdz*velocity(i, j, k - 1);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_5(i, j, k)*inv_dzdz*velocity(i, j, k + 1);
			}		
					
			x_vector[bc(i, j, k)] = velocity(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void BuildLinearSystemNodeForSemiImplicitViscosity(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, const FIELD_STRUCTURE_3D<T>& velocity, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term, const FIELD_STRUCTURE_3D<T>& boundary_levelset)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speed up variables 
		int i_res(velocity.grid.i_res), j_res(velocity.grid.j_res), k_res(velocity.grid.k_res), i_start(velocity.i_start), i_end(velocity.i_end), j_start(velocity.j_start), j_end(velocity.j_end), k_start(velocity.k_start), k_end(velocity.k_end);
		T dx(velocity.dx), dy(velocity.dy), one_over_dx(velocity.one_over_dx), one_over_dy(velocity.one_over_dy), one_over_dz(velocity.one_over_dz), one_over_dx2(velocity.one_over_dx2), one_over_dy2(velocity.one_over_dy2), one_over_dz2(velocity.one_over_dz2);
		T x_min(velocity.grid.x_min), y_min(velocity.grid.y_min), z_min(velocity.grid.z_min), x_max(velocity.grid.x_max), y_max(velocity.grid.y_max), z_max(velocity.grid.z_max);

		const T dxdx = POW2(velocity.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(velocity.grid.dy), inv_dydy = (T)1/dydy, dzdz = POW2(velocity.grid.dz), inv_dzdz = (T)1/dzdz;
		
		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j, k) > -1)
			{
				coef_ijk += coef_2(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef_2(i, j, k)*inv_dxdx);
			}
			if (bc(i + 1, j, k) > -1)
			{
				coef_ijk += coef_1(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef_1(i, j, k)*inv_dxdx);
			}
			if (bc(i, j - 1, k) > -1)
			{
				coef_ijk += coef_4(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef_4(i, j, k)*inv_dydy);
			}
			if (bc(i, j + 1, k) > -1)
			{
				coef_ijk += coef_3(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef_3(i, j, k)*inv_dydy);
			}
			if (bc(i, j, k - 1) > -1)
			{
				coef_ijk += coef_6(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef_6(i, j, k)*inv_dzdz);
			}
			if (bc(i, j, k + 1) > -1)
			{
				coef_ijk += coef_5(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef_5(i, j, k)*inv_dzdz);
			}

			if (bc(i - 1, j, k) == BC_DIR)
			{
				T theta = abs(boundary_levelset(i, j, k))/(abs(boundary_levelset(i, j, k)) + abs(boundary_levelset(i - 1, j, k)));
				coef_ijk += (theta - 1)/theta*coef_2(i, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				T theta = abs(boundary_levelset(i, j, k))/(abs(boundary_levelset(i, j, k)) + abs(boundary_levelset(i + 1, j, k)));
				coef_ijk += (theta - 1)/theta*coef_1(i, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += coef_4(i, j, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += coef_3(i, j, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				T theta = abs(boundary_levelset(i, j, k))/(abs(boundary_levelset(i, j, k)) + abs(boundary_levelset(i, j, k - 1)));
				coef_ijk += (theta - 1)/theta*coef_6(i, j, k);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				T theta = abs(boundary_levelset(i, j, k))/(abs(boundary_levelset(i, j, k)) + abs(boundary_levelset(i, j, k + 1)));
				coef_ijk += (theta - 1)/theta*coef_5(i, j, k);
			}
			
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), (T)1 + inv_dxdx*(T)coef_ijk);
			
			b_vector[bc(i, j, k)] = explicit_term(i, j, k);

			if (bc(i - 1, j, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += 0;//coef_2(i, j, k)*inv_dxdx*velocity(i - 1, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += 0;//coef_1(i, j, k)*inv_dxdx*velocity(i + 1, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_4(i, j, k)*inv_dydy*velocity(i, j - 1, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_3(i, j, k)*inv_dydy*velocity(i, j + 1, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += 0;//coef_6(i, j, k)*inv_dzdz*velocity(i, j, k - 1);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += 0;//coef_5(i, j, k)*inv_dzdz*velocity(i, j, k + 1);
			}		
					
			x_vector[bc(i, j, k)] = velocity(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void BuildLinearSystemNodeForSemiImplicitViscosity(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, const FIELD_STRUCTURE_3D<T>& velocity, const FIELD_STRUCTURE_3D<T>& velocity_x, const FIELD_STRUCTURE_3D<T>& velocity_y, const FIELD_STRUCTURE_3D<T>& velocity_z, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speed up variables 
		int i_res(velocity.grid.i_res), j_res(velocity.grid.j_res), k_res(velocity.grid.k_res), i_start(velocity.i_start), i_end(velocity.i_end), j_start(velocity.j_start), j_end(velocity.j_end), k_start(velocity.k_start), k_end(velocity.k_end);
		T dx(velocity.dx), dy(velocity.dy), one_over_dx(velocity.one_over_dx), one_over_dy(velocity.one_over_dy), one_over_dz(velocity.one_over_dz), one_over_dx2(velocity.one_over_dx2), one_over_dy2(velocity.one_over_dy2), one_over_dz2(velocity.one_over_dz2);
		T x_min(velocity.grid.x_min), y_min(velocity.grid.y_min), z_min(velocity.grid.z_min), x_max(velocity.grid.x_max), y_max(velocity.grid.y_max), z_max(velocity.grid.z_max);

		const T dxdx = POW2(velocity.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(velocity.grid.dy), inv_dydy = (T)1/dydy, dzdz = POW2(velocity.grid.dz), inv_dzdz = (T)1/dzdz;
		
		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j, k) > -1)
			{
				coef_ijk += coef_2(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef_2(i, j, k)*inv_dxdx);
			}
			if (bc(i + 1, j, k) > -1)
			{
				coef_ijk += coef_1(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef_1(i, j, k)*inv_dxdx);
			}
			if (bc(i, j - 1, k) > -1)
			{
				coef_ijk += coef_4(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef_4(i, j, k)*inv_dydy);
			}
			if (bc(i, j + 1, k) > -1)
			{
				coef_ijk += coef_3(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef_3(i, j, k)*inv_dydy);
			}
			if (bc(i, j, k - 1) > -1)
			{
				coef_ijk += coef_6(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef_6(i, j, k)*inv_dzdz);
			}
			if (bc(i, j, k + 1) > -1)
			{
				coef_ijk += coef_5(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef_5(i, j, k)*inv_dzdz);
			}

			// Periodic Boundary Condition 
			if (bc(i - 1, j, k) == BC_PER)
			{
				coef_ijk += coef_2(i, j, k);
				if (velocity.is_x_component)
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i_end + (i - 1 - (i_start + 1)), j, k), -coef_2(i, j, k)*inv_dxdx);
				}
				else
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i_end + (i - 1 - i_start), j, k), -coef_2(i, j, k)*inv_dxdx);
				}
			}
			if (bc(i + 1, j, k) == BC_PER)
			{
				coef_ijk += coef_1(i, j, k);
				if (velocity.is_x_component)
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i_start + (i + 1 - (i_end - 1)), j, k), -coef_1(i, j, k)*inv_dxdx);
				}
				else
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i_start + (i + 1 - i_end), j, k), -coef_1(i, j, k)*inv_dxdx);
				}
			}
			if (bc(i, j - 1, k) == BC_PER)
			{
				coef_ijk += coef_4(i, j, k);
				if (velocity.is_y_component)
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j_end + (j - 1 - (j_start + 1)), k), -coef_4(i, j, k)*inv_dydy);
				}
				else
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j_end + (j - 1 - j_start), k), -coef_4(i, j, k)*inv_dydy);
				}
			}
			if (bc(i, j + 1, k) == BC_PER)
			{
				coef_ijk += coef_3(i, j, k);
				if (velocity.is_y_component)
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j_start + (j + 1 - (j_end - 1)), k), -coef_3(i, j, k)*inv_dydy);
				}
				else
				{
					A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j_start + (j + 1 - j_end), k), -coef_3(i, j, k)*inv_dydy);
				}
			}

			if (bc(i - 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_2(i, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_1(i, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += coef_4(i, j, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += coef_3(i, j, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				coef_ijk += coef_6(i, j, k);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				coef_ijk += coef_5(i, j, k);
			}
						
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), (T)1 + inv_dxdx*(T)coef_ijk);
			
			b_vector[bc(i, j, k)] = explicit_term(i, j, k);

			if (bc(i - 1, j, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_2(i, j, k)*inv_dxdx*velocity_x(i - 1, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_1(i, j, k)*inv_dxdx*velocity_x(i + 1, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_4(i, j, k)*inv_dydy*velocity_y(i, j - 1, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_3(i, j, k)*inv_dydy*velocity_y(i, j + 1, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_6(i, j, k)*inv_dzdz*velocity_z(i, j, k - 1);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				b_vector[bc(i, j, k)] += coef_5(i, j, k)*inv_dzdz*velocity_z(i, j, k + 1);
			}		
					
			x_vector[bc(i, j, k)] = velocity(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void BuildLinearSystemNodeForSemiImplicitViscosityPeriodic(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, const FIELD_STRUCTURE_3D<T>& velocity, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& explicit_term)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speed up variables 
		int i_res(velocity.grid.i_res), j_res(velocity.grid.j_res), k_res(velocity.grid.k_res), i_start(velocity.i_start), i_end(velocity.i_end), j_start(velocity.j_start), j_end(velocity.j_end), k_start(velocity.k_start), k_end(velocity.k_end);
		T dx(velocity.dx), dy(velocity.dy), one_over_dx(velocity.one_over_dx), one_over_dy(velocity.one_over_dy), one_over_dz(velocity.one_over_dz), one_over_dx2(velocity.one_over_dx2), one_over_dy2(velocity.one_over_dy2), one_over_dz2(velocity.one_over_dz2);
		T x_min(velocity.grid.x_min), y_min(velocity.grid.y_min), z_min(velocity.grid.z_min), x_max(velocity.grid.x_max), y_max(velocity.grid.y_max), z_max(velocity.grid.z_max);

		const T dxdx = POW2(velocity.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(velocity.grid.dy), inv_dydy = (T)1/dydy, dzdz = POW2(velocity.grid.dz), inv_dzdz = (T)1/dzdz;
		
		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j, k) > -1)
			{
				coef_ijk += coef_2(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef_2(i, j, k)*inv_dxdx);
			}
			if (bc(i + 1, j, k) > -1)
			{
				coef_ijk += coef_1(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef_1(i, j, k)*inv_dxdx);
			}
			if (bc(i, j - 1, k) > -1)
			{
				coef_ijk += coef_4(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef_4(i, j, k)*inv_dydy);
			}
			if (bc(i, j + 1, k) > -1)
			{
				coef_ijk += coef_3(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef_3(i, j, k)*inv_dydy);
			}
			if (bc(i, j, k - 1) > -1)
			{
				coef_ijk += coef_6(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef_6(i, j, k)*inv_dzdz);
			}
			if (bc(i, j, k + 1) > -1)
			{
				coef_ijk += coef_5(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef_5(i, j, k)*inv_dzdz);
			}

			if (bc(i - 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_2(i, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_1(i, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += coef_4(i, j, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += coef_3(i, j, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				coef_ijk += coef_6(i, j, k);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				coef_ijk += coef_5(i, j, k);
			}
						
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), (T)1 + inv_dxdx*(T)coef_ijk);
			
			b_vector[bc(i, j, k)] = explicit_term(i, j, k);

			x_vector[bc(i, j, k)] = velocity(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void BuildLinearSystemNodeForSemiImplicitViscosityPeriodic(const int& thread_id, CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_3D<T>& coef_1, const FIELD_STRUCTURE_3D<T>& coef_2, const FIELD_STRUCTURE_3D<T>& coef_3, const FIELD_STRUCTURE_3D<T>& coef_4, const FIELD_STRUCTURE_3D<T>& coef_5, const FIELD_STRUCTURE_3D<T>& coef_6, const FIELD_STRUCTURE_3D<T>& velocity, const FIELD_STRUCTURE_3D<int>& bc)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(thread_id, bc);
		const int nnz = CountNonZeroElements(thread_id, bc);
		
		// Initialize A matrix, x vector, and b vector
		HEAD_THREAD_WORK(A_matrix.Initialize(multithreading, num_all_full_cells, nnz));
		HEAD_THREAD_WORK(x_vector.Initialize(num_all_full_cells, true));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			A_matrix.start_ix[0] = 0;
			A_matrix.end_ix[0] = multithreading->sync_value_int[0] - 1;
			A_matrix.prev_row_array[0] = -1;
			A_matrix.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				A_matrix.start_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + 1;
				A_matrix.end_ix[thread_id] = A_matrix.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				A_matrix.prev_row_array[thread_id] = -1;
				A_matrix.values_ix_array[thread_id] = A_matrix.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speed up variables 
		int i_res(velocity.grid.i_res), j_res(velocity.grid.j_res), k_res(velocity.grid.k_res), i_start(velocity.i_start), i_end(velocity.i_end), j_start(velocity.j_start), j_end(velocity.j_end), k_start(velocity.k_start), k_end(velocity.k_end);
		T dx(velocity.dx), dy(velocity.dy), one_over_dx(velocity.one_over_dx), one_over_dy(velocity.one_over_dy), one_over_dz(velocity.one_over_dz), one_over_dx2(velocity.one_over_dx2), one_over_dy2(velocity.one_over_dy2), one_over_dz2(velocity.one_over_dz2);
		T x_min(velocity.grid.x_min), y_min(velocity.grid.y_min), z_min(velocity.grid.z_min), x_max(velocity.grid.x_max), y_max(velocity.grid.y_max), z_max(velocity.grid.z_max);

		const T dxdx = POW2(velocity.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(velocity.grid.dy), inv_dydy = (T)1/dydy, dzdz = POW2(velocity.grid.dz), inv_dzdz = (T)1/dzdz;
		
		BEGIN_GRID_ITERATION_3D(velocity.partial_grids[thread_id])
		{
			if (bc(i, j, k) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j, k) > -1)
			{
				coef_ijk += coef_2(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i - 1, j, k), -coef_2(i, j, k)*inv_dxdx);
			}
			if (bc(i + 1, j, k) > -1)
			{
				coef_ijk += coef_1(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i + 1, j, k), -coef_1(i, j, k)*inv_dxdx);
			}
			if (bc(i, j - 1, k) > -1)
			{
				coef_ijk += coef_4(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j - 1, k), -coef_4(i, j, k)*inv_dydy);
			}
			if (bc(i, j + 1, k) > -1)
			{
				coef_ijk += coef_3(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j + 1, k), -coef_3(i, j, k)*inv_dydy);
			}
			if (bc(i, j, k - 1) > -1)
			{
				coef_ijk += coef_6(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k - 1), -coef_6(i, j, k)*inv_dzdz);
			}
			if (bc(i, j, k + 1) > -1)
			{
				coef_ijk += coef_5(i, j, k);
				A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k + 1), -coef_5(i, j, k)*inv_dzdz);
			}

			if (bc(i - 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_2(i, j, k);
			}
			if (bc(i + 1, j, k) == BC_DIR)
			{
				coef_ijk += coef_1(i, j, k);
			}
			if (bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += coef_4(i, j, k);
			}
			if (bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += coef_3(i, j, k);
			}
			if (bc(i, j, k - 1) == BC_DIR)
			{
				coef_ijk += coef_6(i, j, k);
			}
			if (bc(i, j, k + 1) == BC_DIR)
			{
				coef_ijk += coef_5(i, j, k);
			}
						
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(thread_id, bc(i, j, k), bc(i, j, k), (T)1 + inv_dxdx*(T)coef_ijk);
			
			x_vector[bc(i, j, k)] = velocity(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void VectorToGrid(const int& thread_id, FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<int>& bc)
	{
		VectorToGrid(thread_id, x, pressure, bc);
	}

	void VectorToGrid(const int& thread_id, const VECTOR_ND<T>& x, FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<int>& bc)
	{
		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			if (pressure.fixed(i, j, k) == true)
			{
				continue;
			}
			else
			{
				if (bc(i, j, k) < 0)
				{
					continue;
				}
				pressure(i, j, k) = x[bc(i, j, k)];
			}
		}
		END_GRID_ITERATION_3D;
	}

	void GridToVector(const int& thread_id, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<T>& div, const FIELD_STRUCTURE_3D<int>& bc)
	{
		GridToVector(thread_id, pressure, div, bc, x, b);
	}

	void GridToVector(const int& thread_id, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<T>& div, const FIELD_STRUCTURE_3D<int>& bc, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector)
	{
		int gridtovectorindex;

		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			if(bc(i, j, k) < 0)
				continue;

			gridtovectorindex = bc(i, j, k);

			b_vector[gridtovectorindex] = -div(i, j, k);

			x_vector[gridtovectorindex] = pressure(i, j, k);
		}
		END_GRID_ITERATION_3D;
	}

	void PrecomputeGSCoefficients(const int& thread_id, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& density)
	{
		if (place_dirichlet_bc_at_face)
		{
			PrecomputeGSCoefficientsFaceDirichlet(thread_id, bc, density);
		}
		else
		{
			PrecomputeGSCoefficientsNodeDirichlet(thread_id, bc, density);
		}
	}

	void PrecomputeGSCoefficientsNodeDirichlet(const int& thread_id, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& density)
	{
		// Speedup constants and references
		const T dx = bc.grid.dx, inv_dx = (T)1/dx, half_dx = (T)0.5*dx;
		const T dxdx = dx*dx, inv_dxdx = (T)1/dxdx;
		const ARRAY_3D<int>& bc_arr(bc.array_for_this);
		const int& i_res(bc_arr.i_res), ij_res(bc_arr.ij_res), ijk_res(bc_arr.ijk_res);

		T parrsum(0), den_c(0), inv_dx_inv_den_ave(inv_dx), coefsum(0);
		const int dixiev[3] = {1, i_res, ij_res};				// Dimensional index deviation

		BEGIN_GRID_ITERATION_3D(gs_coeffients.partial_grids[thread_id])
		{
			GS_COEFFICIENTS& gscoef(gs_coeffients(i, j, k));
			
			gscoef.Initialize();

			if(bc_arr(i, j, k) < 0) continue;

			const T density_center = (T)1;					// const T density_center = density(i, j, k)

			T coef_ijk = 0;

			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.plus[0] = coef;
			}

			if(bc(i - 1, j, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.minus[0] = coef;
			}

			if(bc(i, j + 1, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.plus[1] = coef;
			}

			if(bc(i, j - 1, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.minus[1] = coef;
			}

			if(bc(i, j, k + 1) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.plus[2] = coef;
			}

			if(bc(i, j, k - 1) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.minus[2] = coef;
			}

			// If neighbor is Dirichlet bc (+=2 means bc is placed between cells)
			// We assume density of Dirichlet bc cells is density_ijk
			if(bc(i + 1, j, k) == BC_DIR)
			{
				coef_ijk += (T)1*inv_dxdx;
			}

			if(bc(i - 1, j, k) == BC_DIR)
			{
				coef_ijk += (T)1*inv_dxdx;
			}

			if(bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += (T)1*inv_dxdx;
			}

			if(bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += (T)1*inv_dxdx;
			}

			if(bc(i, j, k + 1) == BC_DIR)
			{
				coef_ijk += (T)1*inv_dxdx;
			}

			if(bc(i, j, k - 1) == BC_DIR)
			{
				coef_ijk += (T)1*inv_dxdx;
			}

			gscoef.center = -coef_ijk;

			if (gscoef.center != (T)0)
			{
				gscoef.inv_center = (T)1/-coef_ijk;
			}
			else
			{
				gscoef.inv_center = (T)0;
			}
		}
		END_GRID_ITERATION_3D;
	}
	
	void PrecomputeGSCoefficientsFaceDirichlet(const int& thread_id, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& density)
	{
		// Speedup constants and references
		const T dx = bc.grid.dx, inv_dx = (T)1/dx, half_dx = (T)0.5*dx;
		const T dxdx = dx*dx, inv_dxdx = (T)1/dxdx;
		const ARRAY_3D<int>& bc_arr(bc.array_for_this);
		const int& i_res(bc_arr.i_res), ij_res(bc_arr.ij_res), ijk_res(bc_arr.ijk_res);

		T parrsum(0), den_c(0), inv_dx_inv_den_ave(inv_dx), coefsum(0);
		const int dixdev[3] = {1, i_res, ij_res};				// Dimensional index deviation

		BEGIN_GRID_ITERATION_3D(gs_coeffients.partial_grids[thread_id])
		{
			GS_COEFFICIENTS& gscoef(gs_coeffients(i, j, k));
			
			gscoef.Initialize();

			if(bc_arr(i, j, k) < 0) continue;

			const T density_center = (T)1;					// const T density_center = density(i, j, k)

			T coef_ijk = 0;

			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.plus[0] = coef;
			}

			if(bc(i - 1, j, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.minus[0] = coef;
			}

			if(bc(i, j + 1, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.plus[1] = coef;
			}

			if(bc(i, j - 1, k) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.minus[1] = coef;
			}

			if(bc(i, j, k + 1) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.plus[2] = coef;
			}

			if(bc(i, j, k - 1) > -1)
			{
				const T coef = inv_dxdx;
				coef_ijk += coef;
				gscoef.minus[2] = coef;
			}

			// If neighbor is Dirichlet bc (+=2 means bc is placed between cells)
			// We assume density of Dirichlet bc cells is density_ijk
			if(bc(i + 1, j, k) == BC_DIR)
			{
				coef_ijk += (T)2*inv_dxdx;
			}

			if(bc(i - 1, j, k) == BC_DIR)
			{
				coef_ijk += (T)2*inv_dxdx;
			}

			if(bc(i, j + 1, k) == BC_DIR)
			{
				coef_ijk += (T)2*inv_dxdx;
			}

			if(bc(i, j - 1, k) == BC_DIR)
			{
				coef_ijk += (T)2*inv_dxdx;
			}

			if(bc(i, j, k + 1) == BC_DIR)
			{
				coef_ijk += (T)2*inv_dxdx;
			}

			if(bc(i, j, k - 1) == BC_DIR)
			{
				coef_ijk += (T)2*inv_dxdx;
			}

			gscoef.center = -coef_ijk;

			if (gscoef.center != (T)0)
			{
				gscoef.inv_center = (T)1/-coef_ijk;
			}
			else
			{
				gscoef.inv_center = (T)0;
			}
		}
		END_GRID_ITERATION_3D;
	}

	void PrecomputeGSCoefficients2ndOrderFreesurface(const int& thread_id, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& density, const LEVELSET_3D& levelset)
	{
		// Speedup constants and references
		const T dx = bc.grid.dx, inv_dx = (T)1/dx, half_dx = (T)0.5*dx;
		const T dxdx = dx*dx, inv_dxdx = (T)1/dxdx;
		const ARRAY_3D<int>& bc_arr(bc.array_for_this);
		const int& i_res(bc_arr.i_res), ij_res(bc_arr.ij_res), ijk_res(bc_arr.ijk_res);

		T parrsum(0), den_c(0), inv_dx_inv_den_ave(inv_dx), coefsum(0);
		const int dixiev[3] = {1, i_res, ij_res};				// Dimensional index deviation

		BEGIN_GRID_ITERATION_3D(gs_coeffients.partial_grids[thread_id])
		{
			GS_COEFFICIENTS& gscoef(gs_coeffients(i, j, k));
			
			gscoef.Initialize();

			if(bc_arr(i, j, k) < 0) continue;

			const T density_center = (T)1;					// const T density_center = density(i, j, k)

			T coef_ijk = 0;

			// If neighbor is full cell
			if(bc(i + 1, j, k) > -1)
			{
				coef_ijk += inv_dxdx;
				gscoef.plus[0] = inv_dxdx;
			}

			if(bc(i - 1, j, k) > -1)
			{
				coef_ijk += inv_dxdx;
				gscoef.minus[0] = inv_dxdx;
			}

			if(bc(i, j + 1, k) > -1)
			{
				coef_ijk += inv_dxdx;
				gscoef.plus[1] = inv_dxdx;
			}

			if(bc(i, j - 1, k) > -1)
			{
				coef_ijk += inv_dxdx;
				gscoef.minus[1] = inv_dxdx;
			}

			if(bc(i, j, k + 1) > -1)
			{
				coef_ijk += inv_dxdx;
				gscoef.plus[2] = inv_dxdx;
			}

			if(bc(i, j, k - 1) > -1)
			{
				coef_ijk += inv_dxdx;
				gscoef.minus[2] = inv_dxdx;
			}

			const T phi_ijk = MAX(ABS(levelset.phi(i, j, k)), (T)1e-4);

			if(bc(i + 1, j, k) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i + 1, j, k)), (T)1e-4);
				coef_ijk += (phi_ijk + phi)/phi_ijk*inv_dxdx;
			}

			if(bc(i - 1, j, k) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i - 1, j, k)), (T)1e-4);
				coef_ijk += (phi_ijk + phi)/phi_ijk*inv_dxdx;
			}

			if(bc(i, j + 1, k) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i, j + 1, k)), (T)1e-4);
				coef_ijk += (phi_ijk + phi)/phi_ijk*inv_dxdx;
			}

			if(bc(i, j - 1, k) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i, j - 1, k)), (T)1e-4);
				coef_ijk += (phi_ijk + phi)/phi_ijk*inv_dxdx;
			}

			if(bc(i, j, k + 1) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i, j, k + 1)), (T)1e-4);
				coef_ijk += (phi_ijk + phi)/phi_ijk*inv_dxdx;
			}

			if(bc(i, j, k - 1) == BC_DIR)
			{
				const T phi = MAX(ABS(levelset.phi(i, j, k - 1)), (T)1e-4);
				coef_ijk += (phi_ijk + phi)/phi_ijk*inv_dxdx;
			}

			gscoef.center = -coef_ijk;

			if (gscoef.center != (T)0)
			{
				gscoef.inv_center = (T)1/-coef_ijk;
			}
			else
			{
				gscoef.inv_center = (T)0;
			}
		}
		END_GRID_ITERATION_3D;
	}

	void WeightedJacobiSmoothing(const int& thread_id, FIELD_STRUCTURE_3D<T>& pressure, FIELD_STRUCTURE_3D<T>& pressure_temp, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div, const int& max_itr, const T& weight)
	{
		ARRAY_3D<T>& p_arr(pressure.array_for_this);

		pressure_temp.CopyAllValuesGhostFrom(thread_id, pressure);

		int*				bc_val = bc.array_for_this.values;
		T*					p_val = pressure.array_for_this.values;
		T*					p_val_tmp = pressure_temp.array_for_this.values;
		T*					div_val = div.array_for_this.values;
		GS_COEFFICIENTS*	gscoef_val(gs_coeffients.array_for_this.values);

		GRID_STRUCTURE_3D&	grid(pressure.partial_grids[thread_id]);
		const int			i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start);
		const int			i_end(grid.i_end),	   j_end(grid.j_end),	  k_end(grid.k_end);
		const int&			i_res(pressure.array_for_this.i_res), ij_res(pressure.array_for_this.ij_res);

		T sqr_res_sum((T)0), max_abs_res((T)0);
		int num_full_cell(0);
		int itr;
		for (itr = 0; itr < max_itr; itr++)
		{
			pressure_temp.CopyAllValuesGhostFrom(thread_id, pressure);
			
			sqr_res_sum = (T)0;
			max_abs_res = (T)0;
			num_full_cell = 0;

			for (int k = k_start; k <= k_end; ++k)
			{
				for (int j = j_start; j < j_end; ++j)
				{
					for (int i = i_start; i < i_end; ++i)
					{
						const int bix(p_arr.Index1D(i, j, k));
						if(bc_val[bix] < 0) continue;

						const GS_COEFFICIENTS& gscoef(gscoef_val[bix]);

						const T res(div_val[bix] - (*(p_val_tmp + bix + 1))*gscoef.plus[0] - (*(p_val_tmp + bix + i_res))*gscoef.plus[1] - (*(p_val_tmp + bix + ij_res))*gscoef.plus[2]
												 - (*(p_val_tmp + bix - 1))*gscoef.minus[0] - (*(p_val_tmp + bix - i_res))*gscoef.minus[1] - (*(p_val_tmp + bix - ij_res))*gscoef.minus[2]
												 - (*(p_val_tmp + bix))*gscoef.center);
						*(p_val + bix) = *(p_val_tmp + bix) + weight*res*gscoef.inv_center;

						sqr_res_sum += POW2(res);
						max_abs_res = MAX(max_abs_res, ABS(res));
						num_full_cell++;
					}
				}
			}

			multithreading->SyncSum(thread_id, sqr_res_sum);
			multithreading->SyncMax(thread_id, max_abs_res);
			multithreading->SyncSum(thread_id, num_full_cell);

			if(sqr_res_sum/(T)num_full_cell <= sqr_tolerance) break;

			multithreading->Sync(thread_id);			// Note: This Sync() prevent crashes which happens when a certain fast thread reaches "sqr_res_sum = (T)0;" above before all threads finishes tolerance check above.
		}
		
		BEGIN_HEAD_THREAD_WORK
		{
			smoother_id++;
			if (smoother_id == num_smoother_id)
			{
				smoother_id = 0;
			}
		}
		END_HEAD_THREAD_WORK
	}

	void SetTolerance(const T& tolerance_input)
	{
		tolerance = tolerance_input;
		sqr_tolerance = tolerance*tolerance;

		linear_solver->SetTolerance(tolerance);
	}

	void CalculateResidual(const int& thread_id, const FIELD_STRUCTURE_3D<T>& residual, const FIELD_STRUCTURE_3D<T>& pressure, const FIELD_STRUCTURE_3D<int>& bc, const FIELD_STRUCTURE_3D<T>& div, const FIELD_STRUCTURE_3D<T>& den)
	{
		// Note: This method calculates residuals of fluid cells only.
		//		 Ghost values of pressure field should be filled before calling this method

		// Speedup constants and references
		const ARRAY_3D<T>&		res_arr(residual.array_for_this);
		const int*				bc_val = bc.array_for_this.values;
		const ARRAY_3D<T>&		parr(pressure.array_for_this);
		const T*				p_val = pressure.array_for_this.values;
		const T*				div_val = div.array_for_this.values;
		const GS_COEFFICIENTS*	gscoef_val(gs_coeffients.array_for_this.values);
		const int&				i_res(pressure.array_for_this.i_res), ij_res(pressure.array_for_this.ij_res);

		BEGIN_GRID_ITERATION_3D(pressure.partial_grids[thread_id])
		{
			const int bix(parr.Index1D(i, j, k));

			if(bc_val[bix] < 0)	continue;

			const GS_COEFFICIENTS& gscoef(gscoef_val[bix]);

			// This refers ghost values of pressure field.
			res_arr(i, j, k) = div_val[bix] - (*(p_val + bix + 1))*gscoef.plus[0] - (*(p_val + bix + i_res))*gscoef.plus[1] - (*(p_val + bix + ij_res))*gscoef.plus[2]
											- (*(p_val + bix - 1))*gscoef.minus[0] - (*(p_val + bix - i_res))*gscoef.minus[1] - (*(p_val + bix - ij_res))*gscoef.minus[2]
											- (*(p_val + bix))*gscoef.center;
		}
		END_GRID_ITERATION_3D;
	}
};
	