#pragma once

#include "OPENGL_OBJECT_BASE.h"

class OPENGL_VECTORFIELD : public OPENGL_OBJECT_BASE
{
public: // Enumerates
	enum FIELD
	{
		FIELD_X = 0,
		FIELD_Y,
		FIELD_Z,

		FIELD_NUM,
	};

	enum VECTORFIELD_DRAW_TYPE
	{
		VECTORFIELD_HIDE	= 0,
		VECTORFIELD_DRAW_SHOW_X,
		VECTORFIELD_DRAW_SHOW_Y,
		VECTORFIELD_DRAW_SHOW_Z,
	};

public: // Essential Data
	FIELD_STRUCTURE_3D<VT>*			vector_field;
	
	FIELD_STRUCTURE_3D<T>*			vector_field_mac;
	
	GRID_STRUCTURE_3D&				grid;

	int								min[FIELD_NUM];
	int								max[FIELD_NUM];
	int								index[FIELD_NUM];

	int								name_base;

	float							length_scale;

public: // Constructor and Destructor
	OPENGL_VECTORFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_3D<VT>* vector_field_input)
		: OPENGL_OBJECT_BASE(display_name, driver), vector_field(vector_field_input), vector_field_mac(0), grid(vector_field_input->grid), length_scale(1.0)
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;

		RegisterDrawType((int) VECTORFIELD_HIDE, "HIDE");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_X, "DRAW_X (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Y, "DRAW_Y (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Z, "DRAW_Z (-,+,/,*)");
		SetDrawType((int) VECTORFIELD_HIDE);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		min[FIELD_X] = vector_field->i_start_g;
		max[FIELD_X] = vector_field->i_end_g;
		index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

		min[FIELD_Y] = vector_field->j_start_g;
		max[FIELD_Y] = vector_field->j_end_g;
		index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

		min[FIELD_Z] = vector_field->k_start_g;
		max[FIELD_Z] = vector_field->k_end_g;
		index[FIELD_Z] = (min[FIELD_Z] + max[FIELD_Z])/2;
		is_velocity = true;
	}

	OPENGL_VECTORFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_3D<VT>* vector_field_input, FIELD_STRUCTURE_3D<T>* vector_field_mac_input)
		: OPENGL_OBJECT_BASE(display_name, driver), vector_field(vector_field_input), vector_field_mac(vector_field_mac_input), length_scale(1.0), grid(vector_field->grid)
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;

		RegisterDrawType((int) VECTORFIELD_HIDE, "HIDE");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_X, "DRAW_X (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Y, "DRAW_Y (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Z, "DRAW_Z (-,+,/,*)");
		SetDrawType((int) VECTORFIELD_HIDE);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		// Speed up variable
		int i_res(vector_field_mac_input->grid.i_res), j_res(vector_field_mac_input->grid.j_res), k_res(vector_field_mac_input->grid.k_res);

		// MAC grid
		// x-component
		if (i_res > j_res && i_res > k_res)
		{
			min[FIELD_X] = vector_field_mac->i_start_g;
			max[FIELD_X] = vector_field_mac->i_end_g;
			index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

			min[FIELD_Y] = vector_field_mac->j_start_g;
			max[FIELD_Y] = vector_field_mac->j_end_g;
			index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

			min[FIELD_Z] = vector_field_mac->k_start_g;
			max[FIELD_Z] = vector_field_mac->k_end_g;
			index[FIELD_Z] = (min[FIELD_Z] + max[FIELD_Z])/2;
			
			is_velocity_x = true;
		}		
			
		// y-component
		if (j_res > i_res && j_res > k_res)
		{
			min[FIELD_X] = vector_field_mac->i_start_g;
			max[FIELD_X] = vector_field_mac->i_end_g;
			index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

			min[FIELD_Y] = vector_field_mac->j_start_g;
			max[FIELD_Y] = vector_field_mac->j_end_g;
			index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

			min[FIELD_Z] = vector_field_mac->k_start_g;
			max[FIELD_Z] = vector_field_mac->k_end_g;
			index[FIELD_Z] = (min[FIELD_Z] + max[FIELD_Z])/2;
			
			is_velocity_y = true;
		}
		
		// z-component
		if (k_res > i_res && k_res > j_res)
		{
			min[FIELD_X] = vector_field_mac->i_start_g;
			max[FIELD_X] = vector_field_mac->i_end_g;
			index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

			min[FIELD_Y] = vector_field_mac->j_start_g;
			max[FIELD_Y] = vector_field_mac->j_end_g;
			index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

			min[FIELD_Z] = vector_field_mac->k_start_g;
			max[FIELD_Z] = vector_field_mac->k_end_g;
			index[FIELD_Z] = (min[FIELD_Z] + max[FIELD_Z])/2;
			
			is_velocity_z = true;
		}
	}

	~OPENGL_VECTORFIELD(void)
	{}

public: // Initialization Function
	void Initialize()
	{}

	void DrawLine(VI& ijk, VT& ces, VT& vec, int name_cell, bool draw_with_name)
	{
		if (draw_with_name)
		{
			glPushName(name_cell);
		}

		if (!draw_with_name && (driver->GetSelectedNameId() == name_cell))
		{
			glColor3f(0.0f, 1.0f, 0.0f);

			glLineWidth(2.0);
			glBegin(GL_LINES);
				glVertex3f(ces.x, ces.y, ces.z);
				glVertex3f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale), ces.z + (vec.z*length_scale));
			glEnd();
			glLineWidth(1.0);
		}
		else
		{
			glColor3f(0.65f, 0.0f, 0.56f);

			glBegin(GL_LINES);
				glVertex3f(ces.x, ces.y, ces.z);
				glVertex3f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale), ces.z + (vec.z*length_scale));
			glEnd();
		}

		if (draw_with_name)
		{
			glPopName();
		}
	}

	void DrawLine(VI& ijk, VT& ces, T& vec_x, T& vec_y, T& vec_z, int name_cell, bool draw_with_name)
	{
		if (draw_with_name)
		{
			glPushName(name_cell);
		}

		if (!draw_with_name && (driver->GetSelectedNameId() == name_cell))
		{
			glColor3f(0.0f, 1.0f, 0.0f);

			glLineWidth(2.0);
			glBegin(GL_LINES);
				glVertex3f(ces.x, ces.y, ces.z);
				glVertex3f(ces.x + (vec_x*length_scale), ces.y + (vec_y*length_scale), ces.z + (vec_z*length_scale));
			glEnd();
			glLineWidth(1.0);
		}
		else
		{
			glColor3f(0.65f, 0.0f, 0.56f);

			glBegin(GL_LINES);
				glVertex3f(ces.x, ces.y, ces.z);
				glVertex3f(ces.x + (vec_x*length_scale), ces.y + (vec_y*length_scale), ces.z + (vec_z*length_scale));
			glEnd();
		}

		if (draw_with_name)
		{
			glPopName();
		}
	}

	/*void RenderXField(bool draw_with_name = false)
	{
		int i, j, k;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_3D(i, j, k, index[FIELD_X], min[FIELD_Y], min[FIELD_Z], index[FIELD_X], max[FIELD_Y], max[FIELD_Z])
		{
			VT& vec = vector_field->array_for_this(i, j, k);
			VT& ces = (VT)grid.CellCenter(i, j, k);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j, k), ces, vec, name_cell, draw_with_name);

			count_cell++;
		}
	}

	void RenderYField(bool draw_with_name = false)
	{
		int i, j, k;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_3D(i, j, k, min[FIELD_X], index[FIELD_Y], min[FIELD_Z], max[FIELD_X], index[FIELD_Y], max[FIELD_Z])
		{
			VT& vec = vector_field->array_for_this(i, j, k);
			VT& ces = (VT)grid.CellCenter(i, j, k);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j, k), ces, vec, name_cell, draw_with_name);

			count_cell++;
		}
	}

	void RenderZField(bool draw_with_name = false)
	{
		int i, j, k;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_3D(i, j, k, min[FIELD_X], min[FIELD_Y], index[FIELD_Z], max[FIELD_X], max[FIELD_Y], index[FIELD_Z])
		{
			VT& vec = vector_field->array_for_this(i, j, k);
			VT& ces = (VT)grid.CellCenter(i, j, k);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j, k), ces, vec, name_cell, draw_with_name);

			count_cell++;
		}
	}*/

	void RenderXField(bool draw_with_name = false)
	{
		int i, j, k;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_3D(i, j, k, index[FIELD_X], min[FIELD_Y], min[FIELD_Z], index[FIELD_X], max[FIELD_Y], max[FIELD_Z])
		{
			T& x_component = vector_field_mac->array_for_this(i, j, k);
			VT& vec_x = VT(x_component, 0, 0);
			VT& ces = (VT)grid.CellCenter(i, j, k);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j, k), ces, vec_x, name_cell, draw_with_name);

			count_cell++;
		}
	}

	void RenderYField(bool draw_with_name = false)
	{
		int i, j, k;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_3D(i, j, k, min[FIELD_X], index[FIELD_Y], min[FIELD_Z], max[FIELD_X], index[FIELD_Y], max[FIELD_Z])
		{
			T& y_component = vector_field_mac->array_for_this(i, j, k);
			VT& vec_y = VT(0, y_component, 0);
			VT& ces = (VT)grid.CellCenter(i, j, k);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j, k), ces, vec_y, name_cell, draw_with_name);

			count_cell++;
		}
	}
	
	void RenderZField(bool draw_with_name = false)
	{
		int i, j, k;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_3D(i, j, k, min[FIELD_X], min[FIELD_Y], index[FIELD_Z], max[FIELD_X], max[FIELD_Y], index[FIELD_Z])
		{
			T& z_component = vector_field_mac->array_for_this(i, j, k);
			VT& vec_z = VT(0, 0, z_component);
			VT& ces = (VT)grid.CellCenter(i, j, k);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j, k), ces, vec_z, name_cell, draw_with_name);

			count_cell++;
		}
	}

	void IncrementValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		length_scale += 0.1f;

		cout << "Vector field length scale = " << length_scale << endl;
	}

	void DecrementValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		length_scale -= 0.1f;
		if (length_scale < 0.1f)
		{
			length_scale = 0.1f;
		}

		cout << "Vector field length scale = " << length_scale << endl;
	}

	void LeftValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		FIELD field = FIELD_X;
		switch (GetDrawType())
		{
		case VECTORFIELD_DRAW_SHOW_X:
			field = FIELD_X;
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			field = FIELD_Y;
			break;
		case VECTORFIELD_DRAW_SHOW_Z:
			field = FIELD_Z;
			break;
		}

		index[field] -= 1;
		index[field] = CLAMP(index[field], min[field], max[field]);
	}

	void RightValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		FIELD field = FIELD_X;
		switch (GetDrawType())
		{
		case VECTORFIELD_DRAW_SHOW_X:
			field = FIELD_X;
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			field = FIELD_Y;
			break;
		case VECTORFIELD_DRAW_SHOW_Z:
			field = FIELD_Z;
			break;
		}

		index[field] += 1;
		index[field] = CLAMP(index[field], min[field], max[field]);
	}

public: // Virtual Functions
	virtual int NextDrawType()
	{
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_X);
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_Y);
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_Z);
			break;
		case VECTORFIELD_DRAW_SHOW_Z:
			SetDrawType((int) VECTORFIELD_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_Z);
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			SetDrawType((int) VECTORFIELD_HIDE);
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_X);
			break;
		case VECTORFIELD_DRAW_SHOW_Z:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_Y);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{}

	virtual void Render()
	{
		GetDriver()->SetDefaultRenderStatesLineMode();
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			RenderXField();
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			RenderYField();
			break;
		case VECTORFIELD_DRAW_SHOW_Z:
			RenderZField();
			break;
		}
	}

	virtual void RenderWithName()
	{
		GetDriver()->SetDefaultRenderStatesLineMode();
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			RenderXField(true);
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			RenderYField(true);
			break;
		case VECTORFIELD_DRAW_SHOW_Z:
			RenderZField(true);
			break;
		}
	}

	virtual void UserAction(USER_ACTION_TYPE user_action)
	{
		switch (user_action	)
		{
		case OPENGL_OBJECT_BASE::ACTION_INCREMENT:
			IncrementValue();
			break;
		case OPENGL_OBJECT_BASE::ACTION_DECREMENT:
			DecrementValue();
			break;
		case OPENGL_OBJECT_BASE::ACTION_LEFT:
			LeftValue();
			break;
		case OPENGL_OBJECT_BASE::ACTION_RIGHT:
			RightValue();
			break;
		}
	}
};

