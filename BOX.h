#pragma once

#include "PLANE.h"
#include "ARRAY.h"

// Need to check the way of signing index
static const int edgevertex_lookup[12][2] = {{0,4}, {0,2}, {0,3}, {4,7}, {2,7}, {4,5}, {1,7}, {1,5}, {1,6}, {3,5}, {2,6}, {3,6}};

// Rectangular parallelepiped
class BOX : public LEVELSET_OBJECT
{
public: // Essential Data
	VT min;
	VT max;
	
public: // Constructors and Destructor
	BOX(void)
	{}

	BOX(const BOX& box_input)
		: min(box_input.min), max(box_input.max)
	{}

	BOX(const VT& min_input, const VT& max_input)
	{
		Initialize(min_input, max_input);
	}

	BOX(const T& min_x, const T& min_y, const T& min_z, const T& max_x, const T& max_y, const T& max_z)
	{
		Initialize(min_x, min_y, min_z, max_x, max_y, max_z);
	}

	BOX(const T min_input[3], const T max_input[3])
		:min(min_input), max(max_input)
	{}
	
	~BOX(void)
	{}

public: // Initialization Functions
	void Initialize(const VT& min_input, const VT& max_input)
	{
		min.Assign(min_input);
		max.Assign(max_input);
	}

	void Initialize(const T& min_x, const T& min_y, const T& min_z, const T& max_x, const T& max_y, const T& max_z)
	{
		min.Assign(min_x, min_y, min_z);
		max.Assign(max_x, max_y, max_z);
	}

public: // Member Functions
	inline bool Inside(const VT& position, const T& padding_width) const
	{

		if(position.x < min.x + padding_width) return false;
		else if(position.y < min.y + padding_width) return false;
		else if(position.z < min.z + padding_width) return false;
		else if(position.x > max.x - padding_width) return false;
		else if(position.y > max.y - padding_width) return false;
		else if(position.z > max.z - padding_width) return false;
		return true;
	}

	inline bool Inside(const VT& position) const
	{
		if(position.x < min.x) return false;
		else if(position.y < min.y) return false;
		else if(position.z < min.z) return false;
		else if(position.x > max.x) return false;
		else if(position.y > max.y) return false;
		else if(position.z > max.z) return false;
		return true;
	}

	inline void Extend(const VT& position)				
	{
		min.x = MIN(min.x, position.x);
		min.y = MIN(min.y, position.y);
		min.z = MIN(min.z, position.z);

		max.x = MAX(max.x, position.x);
		max.y = MAX(max.y, position.y);
		max.z = MAX(max.z, position.z);
	}

	inline void Enlarge(const T& width) 
	{
		min.x -= width;
		min.y -= width;
		min.z -= width;

		max.x += width;
		max.y += width;
		max.z += width;
	}

	inline void Translate(const VT& deviation)
	{
		min += deviation;
		max += deviation;
	}

	inline void GetCornerVertices(ARRAY<VT>& corner_vertices)
	{
		corner_vertices.Initialize(8);

		corner_vertices[0].Assign(min);
		corner_vertices[1].Assign(max);
		corner_vertices[2].Assign(corner_vertices[0].x, corner_vertices[0].y, corner_vertices[1].z);
		corner_vertices[3].Assign(corner_vertices[0].x, corner_vertices[1].y, corner_vertices[0].z);
		corner_vertices[4].Assign(corner_vertices[1].x, corner_vertices[0].y, corner_vertices[0].z);
		corner_vertices[5].Assign(corner_vertices[1].x, corner_vertices[1].y, corner_vertices[0].z);
		corner_vertices[6].Assign(corner_vertices[0].x, corner_vertices[1].y, corner_vertices[1].z);
		corner_vertices[7].Assign(corner_vertices[1].x, corner_vertices[0].y, corner_vertices[1].z);
	}

	inline bool Intersection(const PLANE& plane)
	{
		ARRAY<VT> corner_vertices;
		GetCornerVertices(corner_vertices);
		ARRAY<T> phi(corner_vertices.num_elements);
		bool negative_phi_exist = false;
		bool if_all_negative = true;
		int first_phi_index = -1;

		for(int i = 0; i < corner_vertices.num_elements; i++)
		{
			phi[i] = plane.SignedDistance(corner_vertices[i]);
			if(phi[i] <= (T)0)
			{
				negative_phi_exist = true;
				if(first_phi_index == -1) first_phi_index = i;
			}
			else if_all_negative = false;
		}

		if(negative_phi_exist == false)						// If there is no intersection, returns zero volume box
		{
			min.AssignZeroVector();
			max.AssignZeroVector();
			return false;
		}

		if(if_all_negative == true) return true;			// If all of this box is inside the plane volume, returns as is.

		min.Assign(corner_vertices[first_phi_index]);
		max.Assign(corner_vertices[first_phi_index]);

		for(int i = 0; i < corner_vertices.num_elements; i++)
		{
			if(phi[i] <= 0) Extend(corner_vertices[i]);
		}

		// Edge Intersection
		for(int edge = 0; edge < 12; edge++)
		{
			const int index0 = edgevertex_lookup[edge][0], index1 = edgevertex_lookup[edge][1];
			const T phi0 = phi[index0], phi1 = phi[index1];
			if((phi0 <= (T)0 && phi1 > (T)0) || (phi1 <= (T)0 && phi0 > (T)0))
			{
				Extend((corner_vertices[index0]*ABS(phi1) + corner_vertices[index1]*ABS(phi0))/(ABS(phi0) + ABS(phi1)));
			}
		}

		return true;
	}

	inline void Scale(const T& s)
	{
		min *= s;
		max *= s;
	}
};
		
		


		

