#pragma once

#include "MATRIX_3X3.h"
#include "LEVELSET_3D.h"
#include "SIMULATION_OBJECT.h"
#include "SCRIPT_READER.h"
#include "SCAN_LINE_ALGORITHM.h"

#include <boost/filesystem.hpp>
#include <vector>

// Based on David Baraff's lecture notes http://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf

class RIGID_OBJECT : public SIMULATION_OBJECT
{
public: // Essential Data
	// Constant quantities
	T					mass;
	T					dm;
	MATRIX_3X3			inertia_body;
	MATRIX_3X3			inverse_inertia_body;
	
	// State variables
	VT					position;

	QUATERNION			rotation;
	VT					linear_momentum;
	VT					angular_momentum;

	VT					velocity;
	VT					angular_velocity;

	// Computed quantities
	VT					force;
	VT					torque;

	// Geometry
	TRIANGULAR_SURFACE*	triangular_surface;
	LEVELSET_3D*		levelset_3d;

	// Options for rigid object
	bool				static_object;
	bool				no_gravity;

	T					restitution;
	T					friction;

	// For spatial
	VT					spiral_center;
	VT					spiral_linear_velocity;
	VT					spiral_angular_velocity;
	T					spiral_outside_velocity;

	T					spiral_radius;
	
	T					spring_coef;
	T					damping_coef;

public: // Constructor and Destructor
	RIGID_OBJECT(void)
		: triangular_surface(0), levelset_3d(0), static_object(false), no_gravity(false), spring_coef((T)100), damping_coef((T)0.001)
	{
		object_type = RIGID;

		mass = (T)0;
		dm = (T)0;
		position = VT(0, 0, 0);

		restitution = (T)0;
		friction = (T)0;

		rotation = MATRIX_3X3::Identity();
		velocity = VT();
		angular_velocity = VT();
		force = VT();
		torque = VT();

		linear_momentum = VT(0, 0, 0);
		angular_momentum = VT(0, 0, 0);

		inertia_body = MATRIX_3X3::Identity();
		inverse_inertia_body = inertia_body.Inversed();
	}

	~RIGID_OBJECT(void)
	{
		DELETE_POINTER(triangular_surface);
		DELETE_POINTER(levelset_3d);
	}

public: // Initialization Functions
	void InitializeFromSolidFile(const char* solid_filename, MULTITHREADING* multithreading_input, const T& mass_input, const GRID_STRUCTURE_3D& grid_input, const int& ghost_cell_width)
	{
		// TODO : Implements to initialize from solid file
	}

	void InitializeFromObjectFile(const char* obj_filename, MULTITHREADING* multithreading_input, const T& mass_input, const GRID_STRUCTURE_3D& grid_input, const int& ghost_cell_width)
	{
		// Need to be updated
	}

	void InitializeFromScriptBlock(MULTITHREADING* multithreading_input, const SCRIPT_BLOCK& script_block, const GRID_STRUCTURE_3D& grid_input, const int& ghost_cell_width)
	{
		const char* obj_filename = script_block.GetString("file", 0);

		cout << "Loading " << obj_filename << endl;

		DELETE_POINTER(triangular_surface);
		triangular_surface = new TRIANGULAR_SURFACE;
		triangular_surface->ReadOBJ(obj_filename);
		triangular_surface->Scale(script_block.GetFloat("obj_scale", (T)1));

		position = script_block.GetVector3("position", VT());

		mass = script_block.GetFloat("mass", (T)1);
		rotation = MATRIX_3X3::Identity();
		velocity = VT();
		angular_velocity = VT();
		force = VT();
		torque = VT();

		spring_coef = script_block.GetFloat("spring_coef", spring_coef);
		damping_coef = script_block.GetFloat("damping_coef", damping_coef);

		restitution = script_block.GetFloat("restitution", (T)0);
		friction = script_block.GetFloat("friction", (T)0);

		spiral_center = script_block.GetVector3("spiral_center", VT());
		spiral_linear_velocity = script_block.GetVector3("spiral_linear_velocity", VT());
		spiral_angular_velocity = script_block.GetVector3("spiral_angular_velocity", VT());
		spiral_outside_velocity = script_block.GetFloat("spiral_outside_velocity", (T)0);

		spiral_radius = (position - spiral_center).Magnitude();

		BOX domain(triangular_surface->bounding_box);
		domain.Enlarge(grid_input.dx*(T)5);
		GRID_STRUCTURE_3D levelset_grid;
		levelset_grid.Initialize(domain, grid_input.dx*script_block.GetFloat("dx_scale", (T)1));

		// Generate levelset 
		levelset_3d = new LEVELSET_3D;
		levelset_3d->Initialize(multithreading_input, levelset_grid, ghost_cell_width);

		T old_value = levelset_grid.dx*(T)-3;
		T new_value = levelset_grid.dx*(T)3;

		levelset_3d->AssignAllValuesLevelsetThreaded(old_value);
		levelset_3d->AssignTriangularSurfaceLevelSet(*triangular_surface);

		SCAN_LINE_ALGORITHM<T> scan_line_algorithm;
		scan_line_algorithm.Initialize(multithreading_input, levelset_3d->grid.ijk_res);

		// Register first seed cells
		for (int t = 0; t < multithreading_input->num_threads; t++)
		{
			GRID_STRUCTURE_3D& bc_grid = levelset_3d->partial_grids[t];
			int i, j, k;
			LOOPS_3D(i, j, k, bc_grid.i_start, bc_grid.j_start, bc_grid.k_start, bc_grid.i_end, bc_grid.j_end, bc_grid.k_end)
			{
				scan_line_algorithm.stack_from_thread[t].Push(VI(i, j, k));
			}
		}

		scan_line_algorithm.ScanFieldNoBlockThreaded(&(levelset_3d->signed_distance_field), old_value, new_value);

		levelset_3d->FullFastSweepingMethodThreaded();
		levelset_3d->UpdateThreaded();

		levelset_object = levelset_3d;

		const VT center_of_mass(ComputeMassCenter());

		// Note : Translate rigid body geometry so that the center of mass = (0, 0, 0) in object coordinates
		levelset_3d->Translate(-center_of_mass);
		triangular_surface->Translate(-center_of_mass.x, -center_of_mass.y, -center_of_mass.z);
		triangular_surface->bounding_box.Translate(-center_of_mass);

		ComputeInertiaTensor();
	}

public: // Member Functions
	VT ComputeMassCenter()
	{
		assert(levelset_3d != 0);

		int inner_cell_count = 0;
		VT center_of_mass((T)0, (T)0, (T)0);

		GRID_ITERATION_3D(levelset_3d->grid)
		{
			if ((*levelset_3d)(i, j, k) <= (T)0)
			{
				inner_cell_count++;
				center_of_mass += levelset_3d->CellCenter(i, j, k);
			}
		}

		if (inner_cell_count > 0)
		{
			center_of_mass /= (T)inner_cell_count;
			dm = mass / (T)inner_cell_count;
		}

		return center_of_mass;
	}

	void ComputeInertiaTensor()
	{
		if (!levelset_3d)
		{
			return;
		}

		T i_xx(0), i_yy(0), i_zz(0), i_xy(0), i_xz(0), i_yz(0);

		GRID_ITERATION_3D(levelset_3d->grid)
		{
			if ((*levelset_3d)(i, j, k) <= (T)0)
			{
				VT center_of_cell = levelset_3d->CellCenter(i, j, k);
				VT correct_position = center_of_cell;

				i_xx += dm*(POW2(correct_position.y) + POW2(correct_position.z));
				i_yy += dm*(POW2(correct_position.z) + POW2(correct_position.x));
				i_zz += dm*(POW2(correct_position.x) + POW2(correct_position.y));

				i_xy += dm*(correct_position.x*correct_position.y);
				i_xz += dm*(correct_position.x*correct_position.z);
				i_yz += dm*(correct_position.y*correct_position.z);
			}
		}

		inertia_body = MATRIX_3X3(i_xx, -i_xy, -i_xz, -i_xz, i_yy, -i_yz, -i_xz, -i_yz, i_zz);
		inverse_inertia_body = inertia_body.Inversed();
	}

	VT ObjectSpaceVector(const VT& world_space_vector) const
	{
		return rotation.Inverse_Rotate(world_space_vector);
	}

	VT WorldSpaceVector(const VT& object_space_vector) const
	{
		return rotation.Rotate(object_space_vector);
	}

	VT WorldVelocity(const VT& sampling_position) const
	{
		return velocity + CrossProduct(angular_velocity, sampling_position - position);
	}

	VT WorldForce(const VT& sampling_position) const
	{
		return force + CrossProduct(torque, sampling_position - position);
	}

	MATRIX_3X3 Star(const VT& omega) const
	{
		const T* a = omega.values;
		return MATRIX_3X3((T)0, -a[2], a[1], a[2], (T)0, -a[0], -a[1], a[0], (T)0);
	}

	VT Object2World(const VT& position_input) const
	{
		return WorldSpaceVector(position_input) + position;
	}

	VT World2Object(const VT& position_input) const
	{
		return ObjectSpaceVector(position_input - position);
	}

	const T SignedDistance(const VT& position_input) const
	{
		assert(levelset_object != 0);

		return levelset_object->SignedDistance(World2Object(position_input));
	}

	const void Normal(const VT& position_input, VT& normal) const
	{
		assert(levelset_object != 0);

		normal = WorldSpaceVector(levelset_object->UnitNormal(World2Object(position_input)));
	}

	const VT Normal(const VT& position_input) const
	{
		assert(levelset_object != 0);

		return WorldSpaceVector(levelset_object->UnitNormal(World2Object(position_input)));
	}

	const void UnitNormal(const VT& position_input, VT& normal) const
	{
		assert(levelset_object != 0);

		normal = WorldSpaceVector(levelset_object->UnitNormal(World2Object(position_input)));
	}

	const VT UnitNormal(const VT& position_input) const
	{
		assert(levelset_object != 0);

		return WorldSpaceVector(levelset_object->UnitNormal(World2Object(position_input)));
	}

	const VT ClosestPoint(const VT& position_input) const
	{
		assert(levelset_object != 0);

		VT levelset_position = World2Object(position_input);
		T closest_distance = levelset_object->SignedDistance(levelset_position);
		VT unit_normal = levelset_object->UnitNormal(levelset_position);
		unit_normal = WorldSpaceVector(unit_normal);
		return position_input - unit_normal*closest_distance;
	}

	const VT Velocity(const VT& position_input) const
	{
		return WorldVelocity(position_input);
	}

	const bool Inside(const VT& position_input) const
	{
		assert(levelset_object != 0);

		if (levelset_object->SignedDistance(World2Object(position_input)) <= (T)0)
		{
			return true;
		}

		return false;
	}

	inline bool InsideAABB(const VT& position) const
	{
		return aabb.Inside(position);
	}

	inline bool InsideOBB(const VT& position) const
	{
		assert(levelset_object != 0);
		return obb.Inside(World2Object(position));
	}

	void EulerUpdate(const T& dt)
	{
		if (static_object == true)
		{
			return;
		}

		// For Spatial Velocity
		VT direction = position - spiral_center;
		VT spiral_velocity = CrossProduct(spiral_velocity, direction);
		T spiral_velocity_magnitude = spiral_velocity.Magnitude();

		if (spiral_velocity_magnitude != (T)0)
		{
			linear_momentum = spiral_linear_velocity*mass;
			linear_momentum += spiral_velocity*mass;
		}

		// Update momentum from forces
		// NOTE : momentum change by collision was handled in collision response part
		linear_momentum += force*dt;
		angular_momentum += torque*dt;

		// Update linear and angular velocities from momentums
		velocity = linear_momentum/mass;
		angular_velocity = WorldSpaceVector(inverse_inertia_body*ObjectSpaceVector(angular_momentum));

		// Update position and orientation from velocities
		position += velocity*dt;
		rotation = QUATERNION::FromRotationVector(dt*angular_velocity)*rotation;
		if (rotation.s != (T)0)
		{
			rotation.Normalize();
		}

		// Reset Forces
		force = VT();
		torque = VT();

		if (spiral_velocity_magnitude != (T)0)
		{
			position = spiral_center + (position - spiral_center).Normalized()*spiral_radius;
			spiral_center += spiral_linear_velocity*dt;
		}
	}

	void VelocityUpdate(const T& dt)
	{
		if (static_object == true)
		{
			return;
		}

		// Updates momentums from objects
		linear_momentum += force*dt;
		angular_momentum += torque*dt;

		// Calculate velocities from momentums
		velocity = linear_momentum/mass;
		angular_velocity = WorldSpaceVector(inverse_inertia_body*ObjectSpaceVector(angular_momentum));

		// Reset forces
		force = VT();
		torque = VT();
	}

	void UpdateBoundingBoxes()
	{
		aabb.Initialize(VT(levelset_3d->grid.min), VT(levelset_3d->grid.max));
		ARRAY<VT> obb_vertices;
		aabb.GetCornerVertices(obb_vertices);

		for (int i = 0; i < obb_vertices.num_elements; i++)
		{
			obb_vertices[i].Assign(Object2World(obb_vertices[i]));
		}

		aabb.Initialize(obb_vertices[0], obb_vertices[0]);

		for (int i = 1; i < obb_vertices.num_elements; i++)
		{
			aabb.Extend(obb_vertices[i]);
		}

		obb.Initialize(levelset_3d->grid.min, levelset_3d->grid.max);
	}

	ofstream& Write(ofstream& os)
	{
		os.write((char*)&position, sizeof(position));
		os.write((char*)&rotation, sizeof(rotation));

		os.write((char*)&linear_momentum, sizeof(linear_momentum));
		os.write((char*)&angular_momentum, sizeof(angular_momentum));

		os.write((char*)&velocity, sizeof(velocity));
		os.write((char*)&angular_velocity, sizeof(angular_velocity));

		return os;
	}

	ifstream& Read(ifstream& is)
	{
		is.read((char*)&position, sizeof(position));
		is.read((char*)&rotation, sizeof(rotation));

		is.read((char*)&linear_momentum, sizeof(linear_momentum));
		is.read((char*)&angular_momentum, sizeof(angular_momentum));

		is.read((char*)&velocity, sizeof(velocity));
		is.read((char*)&angular_velocity, sizeof(velocity));

		return is;
	}
};

class SIMULATION_OBJECT_CONTACT
{
public: // Essential Data
	SIMULATION_OBJECT* object_a;
	SIMULATION_OBJECT* object_b;

	VT contact_normal;
	VT contact_position;

	T restitution;
	T friction;

public: // Constructor and Destructor
	SIMULATION_OBJECT_CONTACT(void)
		: object_a(0), object_b(0), restitution(0), friction(0)
	{}

	void operator = (const SIMULATION_OBJECT_CONTACT& c)
	{
		object_a = c.object_a;
		object_b = c.object_b;

		contact_normal = c.contact_normal;
		contact_position = c.contact_position;
		
		restitution = c.restitution;
		friction = c.friction;
	}
};

class SEQUENCE_OBJECT : public SIMULATION_OBJECT
{
public: // Essential Data
	LEVELSET_3D*							levelset_3d;
	TRIANGULAR_SURFACE*						triangular_surface;

	FIELD_STRUCTURE_3D<VT>					velocity_field;

	string									sequence_data_directory;
	string									sequence_data_directory_name;

	T										sequence_start_time;
	T										sequence_end_time;
	T										sequence_current_time;
	T										sequence_frame_time;

	int										sequence_start_frame;
	int										sequence_end_frame;
	int										sequence_current_frame;

	T										switch_time;

	T										velocity_magnitude;
	T										sequence_velocity_scale;

	VT										scale;
	VT										translation;

	MULTITHREADING*							multithreading;

	DYNAMIC_ARRAY<boost::filesystem::path>	surface_file_list;
	DYNAMIC_ARRAY<boost::filesystem::path>	levelset_file_list;
	DYNAMIC_ARRAY<boost::filesystem::path>	velocity_file_list;

public: // Constructor and Destructor
	SEQUENCE_OBJECT(void)
		: levelset_3d(0), sequence_start_time((T)0), sequence_end_time((T)0), sequence_current_time((T)0), triangular_surface(0),
		  sequence_start_frame(0), sequence_end_frame(0), sequence_current_frame(0), sequence_velocity_scale(1), multithreading(0), switch_time((T)0),
		  scale(VT((T)1, (T)1, (T)1)), translation(VT())
	{}

	~SEQUENCE_OBJECT(void)
	{}

public: // Initialization Function - Don't know this! After you study more, review this one
	void Initialize(MULTITHREADING* multithreading_input, const char* sequence_data_directory_path, const char* sequence_data_directory_name,
					T frame_time, T start_time, T end_time, int start_frame, int end_frame)
	{
		sequence_data_directory = sequence_data_directory_path;
		sequence_data_directory_name = sequence_data_directory_name;

		sequence_start_time = start_time;
		sequence_end_time = end_time;
		sequence_frame_time = frame_time;

		sequence_start_frame = start_frame;
		sequence_end_frame = end_frame;
		sequence_current_frame = start_frame - 1;

		multithreading = multithreading_input;

		boost::filesystem::path p(sequence_data_directory_name);
		if (exists(p) && is_directory(p))
		{
			vector<boost::filesystem::path> file_list;
			copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(file_list));
			sort(file_list.begin(), file_list.end());
			
			int array_size = (int)file_list.size();
			levelset_file_list.Initialize(array_size);
			velocity_file_list.Initialize(array_size);

			vector<boost::filesystem::path>::iterator it(file_list.begin());
			for (; it != file_list.end(); it++)
			{
				boost::filesystem::path file_path = absolute(*it);

				if (file_path.extension().string() == ".levelset" || file_path.extension().string() == ".LEVELSET")
				{
					levelset_file_list.Push(file_path);
				}
				else if (file_path.extension().string() == ".velocity" || file_path.extension().string() == ".VELOCITY")
				{
					velocity_file_list.Push(file_path);
				}
				else if (file_path.extension().string() == ".obj" || file_path.extension().string() == ".OBJ")
				{
					surface_file_list.Push(file_path);
				}
			}
		}

		if (sequence_end_frame == 0)
		{
			sequence_end_frame = MAX(levelset_file_list.num_of_elements, surface_file_list.num_of_elements) - 1;
		}

		if (levelset_file_list.num_of_elements > 0)
		{
			levelset_3d = new LEVELSET_3D;
			levelset_object = levelset_3d;

			levelset_3d->Read(levelset_file_list.values[0].string().data(), multithreading);

			levelset_3d->FullFastSweepingMethodThreaded();
			levelset_3d->FillGhostCellsFromThreaded(&(levelset_3d->signed_distance_field.array_for_this), false);
			levelset_3d->ComputeNormalsThreaded();

			levelset_3d->Translate(translation);

			aabb.Initialize(levelset_3d->grid.min, levelset_3d->grid.max);
			obb.Initialize(levelset_3d->grid.min, levelset_3d->grid.max);
		}

		if (velocity_file_list.num_of_elements > 0)
		{
			velocity_field.Read(velocity_file_list.values[0].string().data(), multithreading);
			velocity_field.Translate(translation);
		}
		
		if (surface_file_list.num_of_elements > 0)
		{
			DELETE_POINTER(triangular_surface);
			triangular_surface = new TRIANGULAR_SURFACE;
			triangular_surface->ReadOBJ(surface_file_list.values[0].string().data(), (T)0, VT((T)1, (T)1, (T)1), scale, translation);
		}
	}

	void AdvanceSequence(const T& dt)
	{
		if ((sequence_current_time += dt) < sequence_start_time)
		{
			return;
		}

		if ((switch_time += dt) < sequence_frame_time)
		{
			return;
		}

		if ((++sequence_current_frame) > sequence_end_frame)
		{
			return;
		}

		switch_time = (T)0;

		if (levelset_file_list.num_of_elements > 0)
		{
			levelset_3d->Read(levelset_file_list.values[sequence_current_frame].string().data(), multithreading);
			levelset_3d->Translate(translation);

			aabb.Initialize(levelset_3d->grid.min, levelset_3d->grid.max);
			obb.Initialize(levelset_3d->grid.min, levelset_3d->grid.max);
		}

		if (velocity_file_list.num_of_elements > 0)
		{
			velocity_field.Read(velocity_file_list.values[sequence_current_frame].string().data(), multithreading);
			velocity_field.Translate(translation);
		}

		if (surface_file_list.num_of_elements > 0)
		{
			DELETE_POINTER(triangular_surface);
			triangular_surface = new TRIANGULAR_SURFACE;
			triangular_surface->ReadOBJ(surface_file_list.values[sequence_current_frame].string().data(), (T)0, VT((T)1, (T)1, (T)1), scale, translation);
		}
	}

	void UpdateLevelset(const int& thread_id)
	{
		if (levelset_3d)
		{
			levelset_3d->FillGhostCellsFromPointer(thread_id, &(levelset_3d->signed_distance_field.array_for_this), false);
			levelset_3d->ComputeNormals(thread_id);
		}
	}

	const T SignedDistance(const VT& position) const
	{
		if (!levelset_3d)
		{
			return (T)0;
		}
		
		return levelset_3d->SignedDistance(position);
	}

	const bool Inside(const VT& position) const
	{
		if (!levelset_3d)
		{
			return false;
		}

		if (levelset_3d->SignedDistance(position) <= 0)
		{
			return true;
		}

		else
		{
			return true;
		}
	}

	const VT ClosestPoint(const VT& position) const
	{
		if (!levelset_3d)
		{
			return VT();
		}
		
		return position - levelset_3d->UnitNormal(position)*levelset_3d->SignedDistance(position);
	}

	const VT Velocity(const VT& position) const
	{
		if (!levelset_3d)
		{
			return VT();
		}

		return velocity_field(position)*sequence_velocity_scale;
	}

	const void Normal(const VT& position, VT& normal) const
	{
		if (!levelset_3d)
		{
			return;
		}
		
		normal = levelset_3d->Normal(position);
	}

	const VT Normal(const VT& position) const
	{
		if (!levelset_3d)
		{
			return VT();
		}
	
		return levelset_3d->Normal(position);
	}

	const void UnitNormal(const VT& position, VT& normal) const
	{
		if (!levelset_3d)
		{
			return;
		}

		levelset_3d->UnitNormal(position, normal);
	}

	const VT UnitNormal(const VT& position) const
	{
		if (!levelset_3d)
		{
			return VT();
		}

		return levelset_3d->UnitNormal(position);
	}
};
