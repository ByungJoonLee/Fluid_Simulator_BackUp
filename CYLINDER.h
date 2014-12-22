#pragma once

#include "COMMON_DEFINITIONS.h"
#include "LEVELSET_OBJECT.h"
#include "TRIANGULAR_SURFACE.h"
#include "LEVELSET_3D.h"
#include "QUATERNION.h"

class CYLINDER : public LEVELSET_OBJECT
{
public: // Essential Data
	TRIANGULAR_SURFACE* surface;
	LEVELSET_3D*		levelset;

public: // Constructor and Destructor
	CYLINDER(void)
		: surface(0), levelset(0)
	{}

	~CYLINDER(void)
	{
		DELETE_POINTER(surface);
		DELETE_POINTER(levelset);
	}

public: // Initialization Functions
	void Initialize(const VT& position, const T& radius, const T& half_length, const int& slices, MULTITHREADING* multithreading, const GRID_STRUCTURE_3D& grid, const int& ghost_cell_width)
	{
		// Generate cylinder surface structure
		DELETE_POINTER(surface);
		surface = new TRIANGULAR_SURFACE();
		surface->Reset();

		T radius_slice = (T)2*PI/(T)slices;

		for (int i = 0; i < slices; i++)
		{
			T theta0 = (T)i*radius_slice;
			T theta1 = ((T)i + 1)*radius_slice;

			VT p0(0, half_length, 0);
			VT p1(radius*cos(theta0), half_length, radius*sin(theta0));
			VT p2(radius*cos(theta1), half_length, radius*sin(theta1));
			VT p3(radius*cos(theta0), -half_length, radius*sin(theta0));
			VT p4(radius*cos(theta1), -half_length, radius*sin(theta1));
			VT p5(0, -half_length, 0);

			p0 += position;
			p1 += position;
			p2 += position;
			p3 += position;
			p4 += position;
			p5 += position;

			VERTEX* v0 = surface->AddVertex(p0.values);
			VERTEX* v1 = surface->AddVertex(p1.values);
			VERTEX* v2 = surface->AddVertex(p2.values);
			VERTEX* v3 = surface->AddVertex(p3.values);
			VERTEX* v4 = surface->AddVertex(p4.values);
			VERTEX* v5 = surface->AddVertex(p5.values);

			TRIANGLE* t0 = surface->AddTriangle(v2, v1, v0);
			TRIANGLE* t1 = surface->AddTriangle(v2, v3, v1);
			TRIANGLE* t2 = surface->AddTriangle(v4, v3, v2);
			TRIANGLE* t3 = surface->AddTriangle(v4, v5, v3);

			t0->DetermineNormal();
			t1->DetermineNormal();
			t2->DetermineNormal();
			t3->DetermineNormal();
		}

		// Generate levelset surface
		DELETE_POINTER(levelset);
		levelset = new LEVELSET_3D();
		levelset->Initialize(multithreading, grid, ghost_cell_width);

		levelset->AssignMinNegativeValuesLevelsetThreaded();
		levelset->AssignTriangularSurfaceLevelSet(*surface);

		levelset->FloodFillNonRecursive(0, 0, 0);

		levelset->FastSweepingMethodThreaded();
		levelset->UpdateThreaded();
	}

	void Initialize(const VT& start_position, const VT& end_position, const T& radius, const int& slices, MULTITHREADING* multithreading, const GRID_STRUCTURE_3D& grid, const int& ghost_cell_width)
	{
		VT start_to_end_vector = end_position - start_position;
		T half_length = start_to_end_vector.Magnitude()*(T)0.5;
		start_to_end_vector.Normalize();

		VT position = start_position + start_to_end_vector*half_length;

		VT rotation_normal = CrossProduct(VT(0, 1, 0), start_to_end_vector);
		rotation_normal.Normalize();

		T rotation_radiun = acos(DotProduct(start_to_end_vector, VT(0, 1, 0)));

		QUATERNION quat = QUATERNION::FromRotationVector(rotation_normal*rotation_radiun);

		// Generate cylinder surface structure
		surface = new TRIANGULAR_SURFACE();
		surface->Reset();

		T radius_slice = (T)2*PI/(T)slices;

		for (int i = 0; i < slices; i++)
		{
			T theta0 = (T)i*radius_slice;
			T theta1 = ((T)i+1)*radius_slice;

			VT p0(0, half_length, 0);
			VT p1(radius*cos(theta0), half_length, radius*sin(theta0));
			VT p2(radius*cos(theta1), half_length, radius*sin(theta1));
			VT p3(radius*cos(theta0), -half_length, radius*sin(theta0));
			VT p4(radius*cos(theta1), -half_length, radius*sin(theta1));
			VT p5(0, -half_length, 0);

			p0 = quat.Rotate(p0);
			p1 = quat.Rotate(p1);
			p2 = quat.Rotate(p2);
			p3 = quat.Rotate(p3);
			p4 = quat.Rotate(p4);
			p5 = quat.Rotate(p5);

			p0 += position;
			p1 += position;
			p2 += position;
			p3 += position;
			p4 += position;
			p5 += position;

			VERTEX* v0 = surface->AddVertex(p0.values);
			VERTEX* v1 = surface->AddVertex(p1.values);
			VERTEX* v2 = surface->AddVertex(p2.values);
			VERTEX* v3 = surface->AddVertex(p3.values);
			VERTEX* v4 = surface->AddVertex(p4.values);
			VERTEX* v5 = surface->AddVertex(p5.values);

			TRIANGLE* t0 = surface->AddTriangle(v2, v1, v0);
			TRIANGLE* t1 = surface->AddTriangle(v2, v3, v1);
			TRIANGLE* t2 = surface->AddTriangle(v4, v3, v2);
			TRIANGLE* t3 = surface->AddTriangle(v4, v5, v3);

			t0->DetermineNormal();
			t1->DetermineNormal();
			t2->DetermineNormal();
			t3->DetermineNormal();
		}

		// Generate levelset surface
		levelset = new LEVELSET_3D();
		levelset->Initialize(multithreading, grid, ghost_cell_width);

		levelset->AssignMinNegativeValuesLevelsetThreaded();
		levelset->AssignTriangularSurfaceLevelSet(*surface);

		levelset->FloodFillNonRecursive(0, 0, 0);

		levelset->FastSweepingMethodThreaded();
		levelset->UpdateThreaded();
	}

public: // Virtual Functions
	virtual const T SignedDistance(const VT& position) const
	{
		return levelset->SignedDistance(position);
	}

	virtual void Normal(const VT& position, VT& normal) const
	{
		normal = levelset->Normal(position);
	}

	virtual const VT Normal(const VT& position) const
	{
		return levelset->Normal(position);
	}

	virtual void UnitNormal(const VT& position, VT& normal) const
	{
		normal = levelset->UnitNormal(position);
	}

	virtual const VT UnitNormal(const VT& position) const
	{
		return levelset->UnitNormal(position);
	}
};













