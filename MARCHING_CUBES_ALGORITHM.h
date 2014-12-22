#pragma once

#include "TRIANGULAR_SURFACE.h"
#include "LEVELSET_3D.h"
#include "MARCHING_CUBES_ALGORITHM_TABLE.h"

const static float node_lookup[8][3] = {{-1, 1, -1}, {1, 1, -1}, {1, -1, -1}, {-1, -1, -1}, {-1, 1, 1}, {1, 1, 1}, {1, -1, 1}, {-1, -1, 1},};
static const int Nodes_of_Edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {0, 3}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {4, 0}, {5, 1}, {2, 6}, {3, 7}};
static const int power_of_two[15] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};

class MARCHING_CUBE
{
public: 
	VERTEX* edgevertices[12];

public: // Constructor and Destructor
	MARCHING_CUBE(void)
	{
		for(int i = 0; i < 12; i++)
			edgevertices[i] = NULL;
	}
};

class MARCHING_CUBES_ALGORITHM
{
public: // Essential Data
	MULTITHREADING* multithreading;
	ARRAY<TRIANGULAR_SURFACE*> surfaces;
	ARRAY<VERTEX*> vertices;
	ARRAY<TRIANGLE*> triangles;

	ARRAY<VERTEX**> x_edge_vertices, y_edge_vertices, z_edge_vertices, x_edge_vertices_k0, y_edge_vertices_k0;

	ARRAY<int> vertex_start_indices, vertex_end_indices;
	ARRAY<int> triangle_start_indices, triangle_end_indices;

	T dx, dy, dz;
	T half_dx, half_dy, half_dz;
	T quater_dx;

	int num_threads;
	ARRAY<VT> min_sizes, max_sizes;
	VT global_translation;
	int resolution_x, resolution_y;
	ARRAY<int> resolutions_z;

	VT deviation;

	T isolevel;

public: // Constructors and Destructor
	MARCHING_CUBES_ALGORITHM(void)
	{}

	MARCHING_CUBES_ALGORITHM(MULTITHREADING* multithreading_input, GRID_STRUCTURE_3D& grid, T grid_scale = (T)1, const T x_dev_input = (T)0, const T y_dev_input = (T)0, const T z_dev_input = (T)0)
		: isolevel((T)0)
	{
		Initialize(multithreading_input, grid, grid_scale, x_dev_input, y_dev_input, z_dev_input);
	}

	~MARCHING_CUBES_ALGORITHM(void)
	{
		if(surfaces.num_elements != 0)
		{
			for(int i = 0; i < surfaces.num_elements; i++)
				DELETE_POINTER(surfaces[i]);
		}
	}
	
public: // Initialization Function
	void Initialize(MULTITHREADING* multithreading_input, GRID_STRUCTURE_3D& grid, T grid_scale = (T)1, const T x_dev_input = (T)0, const T y_dev_input = (T)0, const T z_dev_input = (T)0)
	{
		multithreading = multithreading_input;

		int grid_max_res = MIN3(grid.i_res, grid.j_res, grid.k_res);

		if(multithreading->num_threads < grid_max_res)
			num_threads = multithreading->num_threads;
		else
			num_threads = grid_max_res;

		if(surfaces.num_elements != 0)
		{
			for(int i = 0; i < surfaces.num_elements; i++)
				DELETE_POINTER(surfaces[i]);
		}

		surfaces.Initialize(num_threads, 0);

		vertex_start_indices.Initialize(num_threads);
		vertex_end_indices.Initialize(num_threads);
		triangle_start_indices.Initialize(num_threads);
		triangle_end_indices.Initialize(num_threads);

		x_edge_vertices.Initialize(num_threads, 0);
		y_edge_vertices.Initialize(num_threads, 0);
		z_edge_vertices.Initialize(num_threads, 0);
		x_edge_vertices_k0.Initialize(num_threads, 0);
		y_edge_vertices_k0.Initialize(num_threads, 0);

		ARRAY<GRID_STRUCTURE_3D> partial_grids;

		grid.SplitInZDirection(num_threads, partial_grids);

		dx = grid.dx/grid_scale;
		dy = grid.dy/grid_scale;
		dz = grid.dz/grid_scale;

		half_dx = grid.dx*(T)0.5;
		half_dy = grid.dy*(T)0.5;
		half_dz = grid.dz*(T)0.5;

		quater_dx = MAX3(dx, dy, dz)*4;

		resolution_x = grid.i_res*(int)grid_scale;
		resolution_y = grid.j_res*(int)grid_scale;

		min_sizes.Initialize(num_threads);
		max_sizes.Initialize(num_threads);
		resolutions_z.Initialize(num_threads);

		for(int i = 0; i < num_threads; i++)
		{
			resolutions_z[i] = partial_grids[i].k_res;

			min_sizes[i].x = partial_grids[i].min[0];
			max_sizes[i].x = partial_grids[i].max[0];
			min_sizes[i].y = partial_grids[i].min[1];
			max_sizes[i].y = partial_grids[i].max[1];
			min_sizes[i].z = partial_grids[i].min[2];
			max_sizes[i].z = partial_grids[i].max[2];

			surfaces[i] = new TRIANGULAR_SURFACE();
		}

		deviation.x = x_dev_input;
		deviation.y = y_dev_input;
		deviation.z = z_dev_input;

		global_translation = VT();

		isolevel = (T)0;
	}

public: // Member Functions
	int index(const int i, const int j) const
	{
		return i + resolution_x*j;
	}

	VT center(const int thread_id, const int i, const int j, const int k) const
	{
		VT resolution_offset = VT();
		return VT(min_sizes[thread_id].x + half_dx + dx*(T)i, min_sizes[thread_id].y + half_dy + dy*(T)j, min_sizes[thread_id].z + half_dz + dz*(T)k) + global_translation;
	}

	VT NodePosition(const int thread_id, const int n, const int i, const int j, const int k) const
	{
		VT center_position = center(thread_id, i, j, k);
		VT node_position;

		node_position.x = center_position.x + (T)node_lookup[n][0]*half_dx;
		node_position.y = center_position.y + (T)node_lookup[n][1]*half_dy;
		node_position.z = center_position.z + (T)node_lookup[n][2]*half_dz;

		return node_position + deviation;
	}

	int index_edge_x(int i, int j, int k) const
	{
		return i + resolution_x*j + resolution_x*(resolution_y + 1)*k;
	}

	int index_edge_y(int i, int j, int k) const
	{
		return i + (resolution_x + 1)*j + (resolution_x + 1)*resolution_y*k;
	}

	int index_edge_z(int i, int j) const
	{
		return i + (resolution_x + 1)*j;
	}

	VERTEX* GetEdgeVertex(int thread_id, int e, int i, int j) const
	{
		if(e == 0) return x_edge_vertices[thread_id][index_edge_x(i, j+1, 0)];
		if(e == 1) return y_edge_vertices[thread_id][index_edge_y(i+1, j, 0)];
		if(e == 2) return x_edge_vertices[thread_id][index_edge_x(i, j, 0)];
		if(e == 3) return y_edge_vertices[thread_id][index_edge_y(i, j, 0)];
		if(e == 4) return x_edge_vertices[thread_id][index_edge_x(i, j+1, 1)];
		if(e == 5) return y_edge_vertices[thread_id][index_edge_y(i+1, j, 1)];
		if(e == 6) return x_edge_vertices[thread_id][index_edge_x(i, j, 1)];
		if(e == 7) return y_edge_vertices[thread_id][index_edge_y(i, j, 1)];
		if(e == 8) return z_edge_vertices[thread_id][index_edge_z(i, j+1)];
		if(e == 9) return z_edge_vertices[thread_id][index_edge_z(i+1, j+1)];
		if(e == 10) return z_edge_vertices[thread_id][index_edge_z(i+1, j)];
		if(e == 11) return z_edge_vertices[thread_id][index_edge_z(i, j)];
		
		cout << "Null getEdgeVertex " << e << " " << i << " " << j << endl;
		exit(1);
	}

	void SetEdgeVertex(int thread_id, int e, int i, int j, VERTEX* vertex) const
	{
		if(e == 0) {x_edge_vertices[thread_id][index_edge_x(i, j+1, 0)] = vertex; return;}
		if(e == 1) {y_edge_vertices[thread_id][index_edge_y(i+1, j, 0)] = vertex; return;}
		if(e == 2) {x_edge_vertices[thread_id][index_edge_x(i, j, 0)] = vertex; return;}
		if(e == 3) {y_edge_vertices[thread_id][index_edge_y(i, j, 0)] = vertex; return;}
		if(e == 4) {x_edge_vertices[thread_id][index_edge_x(i, j+1, 1)] = vertex; return;}
		if(e == 5) {y_edge_vertices[thread_id][index_edge_y(i+1, j, 1)] = vertex; return;}
		if(e == 6) {x_edge_vertices[thread_id][index_edge_x(i, j, 1)] = vertex; return;}
		if(e == 7) {y_edge_vertices[thread_id][index_edge_y(i, j, 1)] = vertex; return;}
		if(e == 8) {z_edge_vertices[thread_id][index_edge_z(i, j+1)] = vertex; return;}
		if(e == 9) {z_edge_vertices[thread_id][index_edge_z(i+1, j+1)] = vertex; return;}
		if(e == 10) {z_edge_vertices[thread_id][index_edge_z(i+1, j)] = vertex; return;}
		if(e == 11) {z_edge_vertices[thread_id][index_edge_z(i, j)] = vertex; return;}
		
		cout << "Null setEdgeVertex " << e << " " << i << " " << j << endl;
		exit(1);
	}

	void PolygonizeThreadPeriodic(int& thread_id, LEVELSET_OBJECT* levelset, int periodic_num_x = 1, int periodic_num_y = 1)
	{
		VERTEX **x_edge_vertices_old, **x_edge_vertices_temp;
		VERTEX **y_edge_vertices_old, **y_edge_vertices_temp;

		surfaces[thread_id]->Reset();

		DELETE_ARRAY(x_edge_vertices[thread_id]);
		DELETE_ARRAY(y_edge_vertices[thread_id]);
		
		int num_01 = resolution_x*(resolution_y + 1);
		int num_02 = (resolution_x + 1)*resolution_y;
		int num_03 = resolution_x*(resolution_y + 1)*2;
		int num_04 = (resolution_x + 1)*resolution_y*2;
		int num_05 = (resolution_x + 1)*(resolution_y + 1);

		x_edge_vertices[thread_id] = new VERTEX* [num_03];
		y_edge_vertices[thread_id] = new VERTEX* [num_04];
		z_edge_vertices[thread_id] = new VERTEX* [num_05];

		x_edge_vertices_old = new VERTEX* [num_03];
		y_edge_vertices_old = new VERTEX* [num_04];

		x_edge_vertices_k0[thread_id] = new VERTEX* [num_01];
		y_edge_vertices_k0[thread_id] = new VERTEX* [num_02];

		for(int i = 0; i < num_01; i++)		x_edge_vertices_k0[thread_id][i] = 0;
		for(int i = 0; i < num_02; i++)		y_edge_vertices_k0[thread_id][i] = 0;

		for(int i = 0; i < num_03; i++)		x_edge_vertices_old[i] = 0;
		for(int i = 0; i < num_04; i++)		y_edge_vertices_old[i] = 0;

		for(int k = 0; k < resolutions_z[thread_id]; k++)
		{
			for(int i = 0; i < num_03; i++)		x_edge_vertices[thread_id][i] = 0;
			for(int i = 0; i < num_04; i++)		y_edge_vertices[thread_id][i] = 0;
			for(int i = 0; i < num_05; i++)		z_edge_vertices[thread_id][i] = 0;

			for(int i = 0; i < resolution_x; i++)
				for(int j = 0; j < (resolution_y + 1); j++)
					x_edge_vertices[thread_id][index_edge_x(i, j, 0)] = x_edge_vertices_old[index_edge_x(i, j, 1)];

			for(int i = 0; i < (resolution_x + 1); i++)
				for(int j = 0; j < resolution_y; j++)
					y_edge_vertices[thread_id][index_edge_y(i, j, 0)] = y_edge_vertices_old[index_edge_y(i, j, 1)];


			// Set thread boundary edge
			if(resolutions_z[thread_id] - 1 == k)
			{
				multithreading->Sync(thread_id, num_threads);
				if(thread_id != num_threads - 1)
				{
					for(int i = 0; i < resolution_x; i++)
						for(int j = 0; j < resolution_y; j++)
						{
							x_edge_vertices[thread_id][index_edge_x(i, j, 1)] = x_edge_vertices_k0[thread_id + 1][index_edge_x(i, j, 0)];
							y_edge_vertices[thread_id][index_edge_y(i, j, 1)] = y_edge_vertices_k0[thread_id + 1][index_edge_y(i, j, 0)];
						}
				}
			}

			for(int i = 0; i < resolution_x; i++)
				for(int j = 0; j < resolution_y; j++)
				{
					T phi[8];
					bool is_skip = false;

					// Sample 8 node values at corners of a cube
					for(int node = 0; node < 8; node++)
					{
						int i_per = i%(resolution_x/periodic_num_x), j_per = j%(resolution_y/periodic_num_y);
						
						if (i >= resolution_x/periodic_num_x)
						{
							phi[node] = (*levelset)(NodePosition(thread_id, node, i_per + 1, j, k));
						}
						else if (j >= resolution_y/periodic_num_y)
						{
							phi[node] = (*levelset)(NodePosition(thread_id, node, i, j_per + 1, k));
						}
						else
						{
							phi[node] = (*levelset)(NodePosition(thread_id, node, i, j, k));
						}
						
						if(phi[node] == (T)0) phi[node] = -(T)1e-8;			// To avoid making zero area triangles

						if(phi[node] >= quater_dx)
						{
							is_skip = true;
							break;
						}
					}

					if(is_skip == true) continue;

					// Find cube index
					int cube_index(0);
					if(phi[0] <= (T)0) cube_index |= 1;
					if(phi[1] <= (T)0) cube_index |= 2;
					if(phi[2] <= (T)0) cube_index |= 4;
					if(phi[3] <= (T)0) cube_index |= 8;
					if(phi[4] <= (T)0) cube_index |= 16;
					if(phi[5] <= (T)0) cube_index |= 32;
					if(phi[6] <= (T)0) cube_index |= 64;
					if(phi[7] <= (T)0) cube_index |= 128;

					if(MCTABLES::EdgeTable[cube_index] == 0) continue;

					// Make vertices
					for(int edge = 0; edge < 12; edge++)
					{
						if(GetEdgeVertex(thread_id, edge, i, j) != 0) continue;		// If vertex on this edge was generated by other cubes, don't need to generate it again

						if(MCTABLES::EdgeTable[cube_index] & power_of_two[edge])
						{
							// Node indices of this edge
							const int n0 = Nodes_of_Edges[edge][0];
							const int n1 = Nodes_of_Edges[edge][1];

							// Phi values of edge nodes
							const T phi0 = phi[Nodes_of_Edges[edge][0]];
							const T phi1 = phi[Nodes_of_Edges[edge][1]];

							VT vertex_position = (ABS(phi1)*NodePosition(thread_id, n0, i, j, k) + ABS(phi0)*NodePosition(thread_id, n1, i, j, k))/(ABS(phi0) + ABS(phi1));
							VT vertex_normal = levelset->UnitNormal(vertex_position);

							VERTEX* vertex = new VERTEX(vertex_position.values, vertex_normal.values);

							SetEdgeVertex(thread_id, edge, i, j, vertex);

							surfaces[thread_id]->AddVertex(vertex);
						}
					}

					for(int l = 0; l < 5; l++)
						if(MCTABLES::TriTable[cube_index][l*3 + 0] != -1)
							surfaces[thread_id]->AddTriangle(GetEdgeVertex(thread_id, MCTABLES::TriTable[cube_index][l*3 + 0], i, j), GetEdgeVertex(thread_id, MCTABLES::TriTable[cube_index][l*3 + 1], i, j), GetEdgeVertex(thread_id, MCTABLES::TriTable[cube_index][l*3 + 2], i, j));
				} // End of resolution_x, resolution_y loop

				if(k == 0)
				{
					for(int i = 0; i < resolution_x; i++)
						for(int j = 0; j < (resolution_y + 1); j++)
							x_edge_vertices_k0[thread_id][index_edge_x(i, j, 0)] = x_edge_vertices[thread_id][index_edge_x(i, j, 0)];

					for(int i = 0; i < (resolution_x + 1); i++)
						for(int j = 0; j < resolution_y; j++)
							y_edge_vertices_k0[thread_id][index_edge_y(i, j, 0)] = y_edge_vertices[thread_id][index_edge_y(i, j, 0)];
				}

				x_edge_vertices_temp = x_edge_vertices_old;
				x_edge_vertices_old = x_edge_vertices[thread_id];
				x_edge_vertices[thread_id] = x_edge_vertices_temp;
				x_edge_vertices_temp = 0;

				y_edge_vertices_temp = y_edge_vertices_old;
				y_edge_vertices_old = y_edge_vertices[thread_id];
				y_edge_vertices[thread_id] = y_edge_vertices_temp;
				y_edge_vertices_temp = 0;

		} // End of resolution_z loop
		multithreading->Sync(thread_id, num_threads);

		for (int i = 0; i < num_threads; i++)
		{
			vertex_start_indices[thread_id] = 0;
			vertex_end_indices[thread_id] = 0;
			triangle_start_indices[thread_id] = 0;
			triangle_end_indices[thread_id] = 0;
		}

		for (int i = 0; i < thread_id; i++)
		{
			vertex_start_indices[thread_id] += (int)surfaces[i]->vertices.size();
			triangle_start_indices[thread_id] += (int)surfaces[i]->triangles.size();
		}

		for (int i = 0; i <= thread_id; i++)
		{
			vertex_end_indices[thread_id] += (int)surfaces[i]->vertices.size();
			triangle_end_indices[thread_id] += (int)surfaces[i]->triangles.size();
		}
		multithreading->Sync(thread_id, num_threads);

		if(thread_id == 0)
		{
			vertices.Initialize(vertex_end_indices[num_threads - 1]);
			triangles.Initialize(triangle_end_indices[num_threads - 1]);

			// multithreading->SplitDomainIndex1D(0, triangles_end_indices[num_threads-1])
			const int k_start = 0;
			const int k_res = triangle_end_indices[num_threads - 1];
			const int k_end = k_start + k_res - 1;
			const int quotient = k_res/num_threads;
			const int remainder = k_res%num_threads;

			int k_start_p = k_start;

			for(int i = 0; i < num_threads; i++)
			{
				int k_depth = i < remainder ? (quotient + 1) : quotient;

				multithreading->start_ix_1D[i] = k_start_p;
				multithreading->end_ix_1D[i] = k_start_p + k_depth - 1;
				
				k_start_p += k_depth;
			}
		}
		multithreading->Sync(thread_id, num_threads);

		list<TRIANGLE*>::iterator itr_triangle;
		itr_triangle = surfaces[thread_id]->triangles.begin();
		for(int i = triangle_start_indices[thread_id]; i < triangle_end_indices[thread_id]; i++)
		{
			triangles[i] = (*itr_triangle);
			itr_triangle++;
		}
		multithreading->Sync(thread_id, num_threads);

		int counter = 0;
		for(int i = vertex_start_indices[thread_id]; i < vertex_end_indices[thread_id]; i++)
		{
			vertices[i] = surfaces[thread_id]->vertices[counter];
			surfaces[thread_id]->vertices[counter]->index = i + 1;
			counter++;
		}
		multithreading->Sync(thread_id, num_threads);

		BEGIN_1D_ITERATION
		{
			triangles[p]->CorrectCCW();
			triangles[p]->DetermineNormal();
		}
		multithreading->Sync(thread_id, num_threads);}

		if(thread_id == 0)
		{
			// multithreading->SplitDomainIndex1D(0, vertex_end_indices[num_threads - 1])
			const int k_start = 0;
			const int k_res = vertex_end_indices[num_threads - 1];
			const int k_end = k_start + k_res - 1;
			const int quotient = k_res/num_threads;
			const int remainder = k_res%num_threads;

			int k_start_p = k_start;

			for(int i = 0; i < num_threads; i++)
			{
				int k_depth = i < remainder ? (quotient + 1) : quotient;

				multithreading->start_ix_1D[i] = k_start_p;
				multithreading->end_ix_1D[i] = k_start_p + k_depth - 1;

				k_start_p += k_depth;
			}
		}
		multithreading->Sync(thread_id, num_threads);

		BEGIN_1D_ITERATION
		{
			vertices[p]->DetermineNormal();
		}
		multithreading->Sync(thread_id, num_threads);}

		DELETE_ARRAY(x_edge_vertices[thread_id]);
		DELETE_ARRAY(y_edge_vertices[thread_id]);
		DELETE_ARRAY(z_edge_vertices[thread_id]);
		DELETE_ARRAY(x_edge_vertices_old);
		DELETE_ARRAY(y_edge_vertices_old);
		DELETE_ARRAY(x_edge_vertices_k0[thread_id]);
		DELETE_ARRAY(y_edge_vertices_k0[thread_id]);
	}
	
	void PolygonizeThread(int& thread_id, LEVELSET_OBJECT* levelset)
	{
		VERTEX **x_edge_vertices_old, **x_edge_vertices_temp;
		VERTEX **y_edge_vertices_old, **y_edge_vertices_temp;

		surfaces[thread_id]->Reset();

		DELETE_ARRAY(x_edge_vertices[thread_id]);
		DELETE_ARRAY(y_edge_vertices[thread_id]);

		int num_01 = resolution_x*(resolution_y + 1);
		int num_02 = (resolution_x + 1)*resolution_y;
		int num_03 = resolution_x*(resolution_y + 1)*2;
		int num_04 = (resolution_x + 1)*resolution_y*2;
		int num_05 = (resolution_x + 1)*(resolution_y + 1);

		x_edge_vertices[thread_id] = new VERTEX* [num_03];
		y_edge_vertices[thread_id] = new VERTEX* [num_04];
		z_edge_vertices[thread_id] = new VERTEX* [num_05];

		x_edge_vertices_old = new VERTEX* [num_03];
		y_edge_vertices_old = new VERTEX* [num_04];

		x_edge_vertices_k0[thread_id] = new VERTEX* [num_01];
		y_edge_vertices_k0[thread_id] = new VERTEX* [num_02];

		for(int i = 0; i < num_01; i++)		x_edge_vertices_k0[thread_id][i] = 0;
		for(int i = 0; i < num_02; i++)		y_edge_vertices_k0[thread_id][i] = 0;

		for(int i = 0; i < num_03; i++)		x_edge_vertices_old[i] = 0;
		for(int i = 0; i < num_04; i++)		y_edge_vertices_old[i] = 0;

		for(int k = 0; k < resolutions_z[thread_id]; k++)
		{
			for(int i = 0; i < num_03; i++)		x_edge_vertices[thread_id][i] = 0;
			for(int i = 0; i < num_04; i++)		y_edge_vertices[thread_id][i] = 0;
			for(int i = 0; i < num_05; i++)		z_edge_vertices[thread_id][i] = 0;

			for(int i = 0; i < resolution_x; i++)
				for(int j = 0; j < (resolution_y + 1); j++)
					x_edge_vertices[thread_id][index_edge_x(i, j, 0)] = x_edge_vertices_old[index_edge_x(i, j, 1)];

			for(int i = 0; i < (resolution_x + 1); i++)
				for(int j = 0; j < resolution_y; j++)
					y_edge_vertices[thread_id][index_edge_y(i, j, 0)] = y_edge_vertices_old[index_edge_y(i, j, 1)];


			// Set thread boundary edge
			if(resolutions_z[thread_id] - 1 == k)
			{
				multithreading->Sync(thread_id, num_threads);
				if(thread_id != num_threads - 1)
				{
					for(int i = 0; i < resolution_x; i++)
						for(int j = 0; j < resolution_y; j++)
						{
							x_edge_vertices[thread_id][index_edge_x(i, j, 1)] = x_edge_vertices_k0[thread_id + 1][index_edge_x(i, j, 0)];
							y_edge_vertices[thread_id][index_edge_y(i, j, 1)] = y_edge_vertices_k0[thread_id + 1][index_edge_y(i, j, 0)];
						}
				}
			}

			for(int i = 0; i < resolution_x; i++)
				for(int j = 0; j < resolution_y; j++)
				{
					T phi[8];
					bool is_skip = false;

					// Sample 8 node values at corners of a cube
					for(int node = 0; node < 8; node++)
					{
						phi[node] = (*levelset)(NodePosition(thread_id, node, i, j, k));

						if(phi[node] == (T)0) phi[node] = -(T)1e-8;			// To avoid making zero area triangles

						if(phi[node] >= quater_dx)
						{
							is_skip = true;
							break;
						}
					}

					if(is_skip == true) continue;

					// Find cube index
					int cube_index(0);
					if(phi[0] <= (T)0) cube_index |= 1;
					if(phi[1] <= (T)0) cube_index |= 2;
					if(phi[2] <= (T)0) cube_index |= 4;
					if(phi[3] <= (T)0) cube_index |= 8;
					if(phi[4] <= (T)0) cube_index |= 16;
					if(phi[5] <= (T)0) cube_index |= 32;
					if(phi[6] <= (T)0) cube_index |= 64;
					if(phi[7] <= (T)0) cube_index |= 128;

					if(MCTABLES::EdgeTable[cube_index] == 0) continue;

					// Make vertices
					for(int edge = 0; edge < 12; edge++)
					{
						if(GetEdgeVertex(thread_id, edge, i, j) != 0) continue;		// If vertex on this edge was generated by other cubes, don't need to generate it again

						if(MCTABLES::EdgeTable[cube_index] & power_of_two[edge])
						{
							// Node indices of this edge
							const int n0 = Nodes_of_Edges[edge][0];
							const int n1 = Nodes_of_Edges[edge][1];

							// Phi values of edge nodes
							const T phi0 = phi[Nodes_of_Edges[edge][0]];
							const T phi1 = phi[Nodes_of_Edges[edge][1]];

							VT vertex_position = (ABS(phi1)*NodePosition(thread_id, n0, i, j, k) + ABS(phi0)*NodePosition(thread_id, n1, i, j, k))/(ABS(phi0) + ABS(phi1));
							VT vertex_normal = levelset->UnitNormal(vertex_position);

							VERTEX* vertex = new VERTEX(vertex_position.values, vertex_normal.values);

							SetEdgeVertex(thread_id, edge, i, j, vertex);

							surfaces[thread_id]->AddVertex(vertex);
						}
					}

					for(int l = 0; l < 5; l++)
						if(MCTABLES::TriTable[cube_index][l*3 + 0] != -1)
							surfaces[thread_id]->AddTriangle(GetEdgeVertex(thread_id, MCTABLES::TriTable[cube_index][l*3 + 0], i, j), GetEdgeVertex(thread_id, MCTABLES::TriTable[cube_index][l*3 + 1], i, j), GetEdgeVertex(thread_id, MCTABLES::TriTable[cube_index][l*3 + 2], i, j));
				} // End of resolution_x, resolution_y loop

				if(k == 0)
				{
					for(int i = 0; i < resolution_x; i++)
						for(int j = 0; j < (resolution_y + 1); j++)
							x_edge_vertices_k0[thread_id][index_edge_x(i, j, 0)] = x_edge_vertices[thread_id][index_edge_x(i, j, 0)];

					for(int i = 0; i < (resolution_x + 1); i++)
						for(int j = 0; j < resolution_y; j++)
							y_edge_vertices_k0[thread_id][index_edge_y(i, j, 0)] = y_edge_vertices[thread_id][index_edge_y(i, j, 0)];
				}

				x_edge_vertices_temp = x_edge_vertices_old;
				x_edge_vertices_old = x_edge_vertices[thread_id];
				x_edge_vertices[thread_id] = x_edge_vertices_temp;
				x_edge_vertices_temp = 0;

				y_edge_vertices_temp = y_edge_vertices_old;
				y_edge_vertices_old = y_edge_vertices[thread_id];
				y_edge_vertices[thread_id] = y_edge_vertices_temp;
				y_edge_vertices_temp = 0;

		} // End of resolution_z loop
		multithreading->Sync(thread_id, num_threads);

		for (int i = 0; i < num_threads; i++)
		{
			vertex_start_indices[thread_id] = 0;
			vertex_end_indices[thread_id] = 0;
			triangle_start_indices[thread_id] = 0;
			triangle_end_indices[thread_id] = 0;
		}

		for (int i = 0; i < thread_id; i++)
		{
			vertex_start_indices[thread_id] += (int)surfaces[i]->vertices.size();
			triangle_start_indices[thread_id] += (int)surfaces[i]->triangles.size();
		}

		for (int i = 0; i <= thread_id; i++)
		{
			vertex_end_indices[thread_id] += (int)surfaces[i]->vertices.size();
			triangle_end_indices[thread_id] += (int)surfaces[i]->triangles.size();
		}
		multithreading->Sync(thread_id, num_threads);

		if(thread_id == 0)
		{
			vertices.Initialize(vertex_end_indices[num_threads - 1]);
			triangles.Initialize(triangle_end_indices[num_threads - 1]);

			// multithreading->SplitDomainIndex1D(0, triangles_end_indices[num_threads-1])
			const int k_start = 0;
			const int k_res = triangle_end_indices[num_threads - 1];
			const int k_end = k_start + k_res - 1;
			const int quotient = k_res/num_threads;
			const int remainder = k_res%num_threads;

			int k_start_p = k_start;

			for(int i = 0; i < num_threads; i++)
			{
				int k_depth = i < remainder ? (quotient + 1) : quotient;

				multithreading->start_ix_1D[i] = k_start_p;
				multithreading->end_ix_1D[i] = k_start_p + k_depth - 1;
				
				k_start_p += k_depth;
			}
		}
		multithreading->Sync(thread_id, num_threads);

		list<TRIANGLE*>::iterator itr_triangle;
		itr_triangle = surfaces[thread_id]->triangles.begin();
		for(int i = triangle_start_indices[thread_id]; i < triangle_end_indices[thread_id]; i++)
		{
			triangles[i] = (*itr_triangle);
			itr_triangle++;
		}
		multithreading->Sync(thread_id, num_threads);

		int counter = 0;
		for(int i = vertex_start_indices[thread_id]; i < vertex_end_indices[thread_id]; i++)
		{
			vertices[i] = surfaces[thread_id]->vertices[counter];
			surfaces[thread_id]->vertices[counter]->index = i + 1;
			counter++;
		}
		multithreading->Sync(thread_id, num_threads);

		BEGIN_1D_ITERATION
		{
			triangles[p]->CorrectCCW();
			triangles[p]->DetermineNormal();
		}
		multithreading->Sync(thread_id, num_threads);}

		if(thread_id == 0)
		{
			// multithreading->SplitDomainIndex1D(0, vertex_end_indices[num_threads - 1])
			const int k_start = 0;
			const int k_res = vertex_end_indices[num_threads - 1];
			const int k_end = k_start + k_res - 1;
			const int quotient = k_res/num_threads;
			const int remainder = k_res%num_threads;

			int k_start_p = k_start;

			for(int i = 0; i < num_threads; i++)
			{
				int k_depth = i < remainder ? (quotient + 1) : quotient;

				multithreading->start_ix_1D[i] = k_start_p;
				multithreading->end_ix_1D[i] = k_start_p + k_depth - 1;

				k_start_p += k_depth;
			}
		}
		multithreading->Sync(thread_id, num_threads);

		BEGIN_1D_ITERATION
		{
			vertices[p]->DetermineNormal();
		}
		multithreading->Sync(thread_id, num_threads);}

		DELETE_ARRAY(x_edge_vertices[thread_id]);
		DELETE_ARRAY(y_edge_vertices[thread_id]);
		DELETE_ARRAY(z_edge_vertices[thread_id]);
		DELETE_ARRAY(x_edge_vertices_old);
		DELETE_ARRAY(y_edge_vertices_old);
		DELETE_ARRAY(x_edge_vertices_k0[thread_id]);
		DELETE_ARRAY(y_edge_vertices_k0[thread_id]);
	}

	void PolygonizeUseFillGhostCellsFromThread(int& thread_id, LEVELSET_3D* levelset)
	{
		levelset->FillGhostCellsFrom(thread_id, levelset->phi, false);
	}

	void Polygonize(LEVELSET_OBJECT& levelset)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
			multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_CUBES_ALGORITHM::PolygonizeThread, this, thread_id, &levelset);

		// multithreading->JoinAll()
		for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

		// DeleteAllThread()
		for(int i = 0; i < num_threads; i++)
		{
			if(multithreading->thread_list[i] != 0)
			{
				delete multithreading->thread_list[i];
				multithreading->thread_list[i] = 0;
			}
		}
		multithreading->num_of_waiting_threads = 0;
	}

	void Polygonize(LEVELSET_3D& levelset, bool use_fill_ghost_cells_from)
	{
		if(use_fill_ghost_cells_from == true)
		{
			for(int thread_id = 0; thread_id < num_threads; thread_id++)
				multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_CUBES_ALGORITHM::PolygonizeUseFillGhostCellsFromThread, this, thread_id, &levelset);

			// multithreading->JoinAll();
			for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

			// DeleteAllThreads()
			for(int i = 0; i < num_threads; i++)
			{
				if(multithreading->thread_list[i] != 0)
				{
					delete multithreading->thread_list[i];
					multithreading->thread_list[i] = 0;
				}
			}
			multithreading->num_of_waiting_threads = 0;
		}

		for(int thread_id = 0; thread_id < num_threads; thread_id++)
			multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_CUBES_ALGORITHM::PolygonizeThread, this, thread_id, &levelset);

		// multithreading->JoinAll();
		for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

		// DeleteAllThreads()
		for(int i = 0; i < num_threads; i++)
		{
			if(multithreading->thread_list[i] != 0)
			{
				delete multithreading->thread_list[i];
				multithreading->thread_list[i] = 0;
			}
		}
		multithreading->num_of_waiting_threads = 0;
	}

	void PolygonizePeriodic(LEVELSET_3D& levelset, bool use_fill_ghost_cells_from, int periodic_num_x = 1, int periodic_num_y = 1)
	{
		if(use_fill_ghost_cells_from == true)
		{
			for(int thread_id = 0; thread_id < num_threads; thread_id++)
				multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_CUBES_ALGORITHM::PolygonizeUseFillGhostCellsFromThread, this, thread_id, &levelset);

			// multithreading->JoinAll();
			for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

			// DeleteAllThreads()
			for(int i = 0; i < num_threads; i++)
			{
				if(multithreading->thread_list[i] != 0)
				{
					delete multithreading->thread_list[i];
					multithreading->thread_list[i] = 0;
				}
			}
			multithreading->num_of_waiting_threads = 0;
		}

		for(int thread_id = 0; thread_id < num_threads; thread_id++)
			multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_CUBES_ALGORITHM::PolygonizeThreadPeriodic, this, thread_id, &levelset, periodic_num_x, periodic_num_y);

		// multithreading->JoinAll();
		for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

		// DeleteAllThreads()
		for(int i = 0; i < num_threads; i++)
		{
			if(multithreading->thread_list[i] != 0)
			{
				delete multithreading->thread_list[i];
				multithreading->thread_list[i] = 0;
			}
		}
		multithreading->num_of_waiting_threads = 0;
	}

	bool IsInterfacial(const T phi[8])
	{
		// Check if these 8 phi values are interfacial
		if(phi[0] > (T)0)
		{
			for(int i = 1; i < 8; i++)
				if(phi[i] <= (T)0) return true;
		}
		else 
		{
			for(int i = 1; i < 8; i++)
				if(phi[i] > (T)0) return true;
		}

		return false;
	}
};






