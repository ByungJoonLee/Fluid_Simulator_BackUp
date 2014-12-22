#pragma once

enum	POISSON_SOLVER_TYPE                                     {NO_SOLVER, MULTIGRID, CG, PCG, HYBRID, GS};
enum	SMOOTHER_TYPE											{NO_SMOOTHER, GS_SMOOTHER, NB_GS_SMOOTHER, WEIGHTED_JACOBI_SMOOTHER};

#define BC_FULL  0
#define BC_DIR  -1
#define BC_OBJ  -2
#define BC_NULL -3
#define BC_NEUM	-4
#define BC_PER  -5

#define PI														(T)3.14159265358979323846

#define ABS(a)													((a) > 0 ? (a) : -(a))

#define MIN(a, b)												((a) > (b) ? (b) : (a))
#define MIN3(a, b, c)											(MIN(MIN(a, b), (c)))
#define MIN4(a, b, c, d)										(MIN(MIN3(a, b, c), (d)))
#define MIN5(a, b, c, d, e)										(MIN(MIN4(a, b, c, d), (e)))
#define MIN7(a, b, c, d, e, f, g)								(MIN(MIN4(a, b, c, d), MIN3(e, f, g)))
#define MIN8(a, b, c, d, e, f, g, h)							(MIN(MIN7(a, b, c, d, e, f, g), (h)))

#define MAX(a, b)												((a) > (b) ? (a) : (b))
#define MAX3(a, b, c)											(MAX(MAX(a, b), (c)))
#define MAX4(a, b, c, d)										(MAX(MAX3(a, b, c), (d)))
#define MAX5(a, b, c, d, e)										(MAX(MAX4(a, b, c, d), (e)))
#define MAX7(a, b, c, d, e, f, g)								(MAX(MAX4(a, b, c, d), MAX3(e, f, g)))
#define MAX8(a, b, c, d, e, f, g, h)							(MAX(MAX7(a, b, c, d, e, f, g), (h)))

#define MIN_ABS(a, b)											(ABS(a) > ABS(b) ? (b) : (a))
#define MAX_ABS(a, b)											(ABS(a) > ABS(b) ? (a) : (b))

#define CLAMP(v, min, max)										((v) < (min) ? (min) : ((v) > (max) ? (max) : (v)))

#define SQUARE(a)												((a)*(a))
#define POW3(a)													((a)*(a)*(a))

#ifdef USE_FLOAT_T
	inline float POW2(const float& a)							{return a*a;}
	inline float POW4(const float& a)							{const float a2 = a*a; return a2*a2;}
	inline float POW5(const float& a)							{const float a2 = a*a; return a2*a2*a;}
	inline float POW6(const float& a)							{const float a3 = a*a*a; return a3*a3;}
	inline float POW7(const float& a)							{const float a3 = a*a*a; return a3*a3*a;}
	inline float POW8(const float& a)							{const float a2 = a*a; const float a4 = a2*a2; return a4*a4;}
#else
	inline double POW2(const double& a)							{return (a)*(a);}
	inline double POW4(const double& a)							{const double a2 = (a)*(a); return (a2)*(a2);}
	inline double POW5(const double& a)							{const double a2 = (a)*(a); return (a2)*(a2)*(a);}
	inline double POW6(const double& a)							{const double a3 = (a)*(a)*(a); return (a3)*(a3);}
	inline double POW7(const double& a)							{const double a3 = (a)*(a)*(a); return (a3)*(a3)*(a);}
	inline double POW8(const double& a)							{const double a2 = (a)*(a); const double a4 = (a2)*(a2); return a4*a4;}
#endif

// Note : The case a == 0 returns negative since phi <= 0 means inside
#define SIGN(a)													(a > 0 ? 1 : -1)

#define DELETE_POINTER(pointer)									if(pointer != 0) {delete pointer; pointer = 0;}
#define DELETE_ARRAY(pointer)									if(pointer != 0) {delete [] pointer; pointer = 0;}

#define INCREASING_SORT3(a, b, c, a1, a2, a3)					if(a <= b){											\
																	if(b <= c) {a1 = a; a2 = b; a3 = c;}			\
																	else if(a <= c) {a1 = a; a2 = c; a3 = b;}		\
																	else {a1 = c; a2 = a; a3 = b;}}					\
																else{												\
																	if(a <= c) {a1 = b; a2 = a; a3 = c;}			\
																	else if(b <= c) {a1 = b; a2 = c; a3 = a;}		\
																	else {a1 = c; a3 = b; a3 = a;}}

// Iteration Macros
#define INIT_GRIDRANGE_2D(grid_input_2d, i_start, j_start, i_end, j_end)					const GRID_STRUCTURE_2D& grid(grid_input_2d);\
																							const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);

#define INIT_GRIDRANGE_3D(grid_input_3d, i_start, j_start, k_start, i_end, j_end, k_end)	const GRID_STRUCTURE_3D& grid(grid_input_3d);\
																							const int i_start(grid.i_start), j_start(grid.j_start), k_start(grid.k_start), i_end(grid.i_end), j_end(grid.j_end), k_end(grid.k_end);

#define LOOPS_2D(i, j, i_start, j_start, i_end, j_end)										for((j) = (j_start) ; (j) <= (j_end); ++(j)) for((i) = (i_start); (i) <= (i_end); ++(i))

#define LOOPS_3D(i, j, k, i_start, j_start, k_start, i_end, j_end, k_end)					for((k) = (k_start) ; (k) <= (k_end); ++(k)) for((j) = (j_start); (j) <= (j_end); ++(j)) for((i) = (i_start); (i) <= (i_end); ++(i))

#define GRID_ITERATION_2D(grid_2d_input)													for(int j_start = (grid_2d_input).j_start, j_end = (grid_2d_input).j_end, i_start = (grid_2d_input).i_start, i_end = (grid_2d_input).i_end, i, j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)

#define GRID_ITERATION_3D(grid_3d_input)													for(int k_start = (grid_3d_input).k_start, k_end = (grid_3d_input).k_end, j_start = (grid_3d_input).j_start, j_end = (grid_3d_input).j_end, i_start = (grid_3d_input).i_start, i_end = (grid_3d_input).i_end, i, j, k = k_start; k <= k_end; ++k) for(j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)

#define BEGIN_GRID_ITERATION_2D(grid_2d_input)												{GRID_STRUCTURE_2D& grid_2d_itr(grid_2d_input);																					\
																							 int i, j;																														\
																							 const int j_start = grid_2d_itr.j_start, j_end = grid_2d_itr.j_end, i_start = grid_2d_itr.i_start, i_end = grid_2d_itr.i_end;  \
																							 for(j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)
#define END_GRID_ITERATION_2D																 multithreading->Sync(thread_id);}

#define BEGIN_GRID_ITERATION_3D(grid_3d_input)												{GRID_STRUCTURE_3D& grid_3d_itr(grid_3d_input);																					\
																							 int i, j, k;																														\
																							 const int k_start = grid_3d_itr.k_start, k_end = grid_3d_itr.k_end, j_start = grid_3d_itr.j_start, j_end = grid_3d_itr.j_end, i_start = grid_3d_itr.i_start, i_end = grid_3d_itr.i_end;  \
																							 for(k = k_start; k <= k_end; ++k) for(j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)
#define END_GRID_ITERATION_3D																 multithreading->Sync(thread_id);}
#define END_GRID_ITERATION_MAX_3D(sync_value)												 multithreading->SyncMax(thread_id, sync_value);}
#define END_GRID_ITERATION_SUM(sync_value)													 multithreading->SyncSum(thread_id, sync_value);}

#define PREPARE_FOR_1D_ITERATION(num)														 multithreading->SplitDomainIndex1D(thread_id, 0, num);

#define BEGIN_1D_ITERATION																	{const int p_start(multithreading->start_ix_1D[thread_id]), p_end(multithreading->end_ix_1D[thread_id]);						\
																							 for(int p = p_start; p <= p_end; p++)
#define END_1D_ITERATION																	 multithreading->Sync(thread_id);}

#define BEGIN_1D_ITERATION_ASYNC(num_elements_input, increase_input)						{const int _num_elements = num_elements_input;									\
																							const int _increase = increase_input;											\
																							multithreading->InitializeSyncIndex(thread_id);								\
																							int p;																			\
																							while(true)																		\
																							{																				\
																								multithreading->GetIterationIndex(p, _increase);							\
																								if(p >= _num_elements) break;												\
																								const int _p_start = p, _p_end = MIN(_p_start + _increase, _num_elements);	\
																								for(p = _p_start; p < _p_end; ++p)
#define END_1D_ITERATION_ASYNC																}																				\
																							multithreading->Sync(thread_id);}

#define BEGIN_PARTICLE_ITERATION															{const int particle_start(multithreading->start_ix_1D[thread_id]), particle_end(multithreading->end_ix_1D[thread_id]); \
																							for(int i = particle_start; i <= particle_end; i++)
#define END_PARTICLE_ITERATION																multithreading->Sync(thread_id);}

#define BEGIN_SCOPED_LOCK																	{boost::mutex::scoped_lock lk(multithreading->sync_mutex);
#define END_SCOPED_LOCK																		}

#define BEGIN_HEAD_THREAD_WORK																if(thread_id == 0)
#define END_HEAD_THREAD_WORK																multithreading->Sync(thread_id);

#define HEAD_THREAD_WORK(expression)														if(thread_id == 0){expression;}; multithreading->Sync(thread_id);

#define PRINT_AND_EXIT(variable)															{cout << #variable << " = " << variable << endl; exit(1);}

