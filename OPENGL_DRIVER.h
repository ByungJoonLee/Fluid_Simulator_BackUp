#pragma once

#include "OPENGL_COMMON.h"
#include "COMMON_DEFINITIONS.h"
#include "MARCHING_CUBES_ALGORITHM.h"
#include "MATERIAL_POINT_3D.h"
#include <Cg/cg.h>
#include <math.h>

class OPENGL_DRIVER
{
public: // Essential Data
	CGcontext				cg_context;
	int						selected_name_id;

	static OPENGL_MATERIAL	material_2Dmode;
	static OPENGL_MATERIAL	material_3Dmode;
	static OPENGL_MATERIAL	material_line;
	
public: // Constructor and Destructor
	OPENGL_DRIVER(void)
	{
		
		// Default material
		material_3Dmode.color_material = COLOR_MATERIAL_DIFFUSE_AND_AMBIENT;

		material_2Dmode.is_lighting = false;
		material_2Dmode.z_buffer = COMPARISON_NEVER;
		material_2Dmode.color_material = COLOR_MATERIAL_DIFFUSE_AND_AMBIENT;

		material_line.is_lighting = false;
		material_line.color_material = COLOR_MATERIAL_DIFFUSE_AND_AMBIENT;

		cg_context = cgCreateContext();

		SetSelectedNameId(-1);
	}

	~OPENGL_DRIVER(void)
	{}

public: // Initialization Function
	void InitGenericGL()
	{
		// Lighting
		GLfloat data[4] = {0.0f, 0.0f, 0.0f, 0.0f,};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, data);
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);

		// Anti-alising
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);

		// ETC
		glClearDepth(1.0);

		glDepthFunc(GL_LEQUAL);

		glFrontFace(GL_CW);

		glAlphaFunc(GL_GREATER, 0.0f);

		SetAntialising(ANTIALISING_OFF);
	}
	
public: // Member Functions
	void SetRenderStatesByMaterial(const OPENGL_MATERIAL& material)
	{
		// Color material
		switch(material.color_material)
		{
		case COLOR_MATERIAL_NONE:
			glDisable(GL_COLOR_MATERIAL);
			break;
		case COLOR_MATERIAL_DIFFUSE:
			glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
			break;
		case COLOR_MATERIAL_AMBIENT:
			glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
			break;
		case COLOR_MATERIAL_EMISSIVE:
			glColorMaterial(GL_FRONT_AND_BACK, GL_EMISSION);
			break;
		case COLOR_MATERIAL_SPECULAR:
			glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR);
			break;
		case COLOR_MATERIAL_DIFFUSE_AND_AMBIENT:
			glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
			break;
		};

		if (material.color_material != COLOR_MATERIAL_NONE)
		{
			glEnable(GL_COLOR_MATERIAL);
		}

		GLfloat color[4];
		if ((material.color_material != COLOR_MATERIAL_AMBIENT) && (material.color_material != COLOR_MATERIAL_DIFFUSE_AND_AMBIENT))
		{
			color[0] = material.ambient_color.GetRed();
			color[1] = material.ambient_color.GetGreen();
			color[2] = material.ambient_color.GetBlue();
			color[3] = material.ambient_color.GetAlpha();
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
		}

		if((material.color_material != COLOR_MATERIAL_DIFFUSE) && (material.color_material != COLOR_MATERIAL_DIFFUSE_AND_AMBIENT))
		{
			color[0] = material.diffuse_color.GetRed();
			color[1] = material.diffuse_color.GetGreen();
			color[2] = material.diffuse_color.GetBlue();
			color[3] = material.diffuse_color.GetAlpha();
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
		}

		if(material.color_material != COLOR_MATERIAL_EMISSIVE)
		{
			color[0] = material.emissive_color.GetRed();
			color[1] = material.emissive_color.GetGreen();
			color[2] = material.emissive_color.GetBlue();
			color[3] = material.emissive_color.GetAlpha();
			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, color);
		}

		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);

		GLfloat color2[4] = {0.f, 0.f, 0.f, 1.f};
		if ((material.shininess != 0.0f) && (material.color_material != COLOR_MATERIAL_SPECULAR))
		{
			color2[0] = material.specular_color.GetRed();
			color2[1] = material.specular_color.GetGreen();
			color2[2] = material.specular_color.GetBlue();
			color2[3] = material.specular_color.GetAlpha();
		}

		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color2);

		// Polygon Mode
		switch(material.polygon_mode)
		{
			case POLYGON_SOLID:
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				break;
			case POLYGON_WIREFRAME:
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				break;
			case POLYGON_POINT:
				glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
				break;
		}

		// Shade Mode
		if (material.is_gouraud_shading)
		{
			glShadeModel(GL_SMOOTH);
		}
		else
		{
			glShadeModel(GL_FLAT);
		}

		// Lighting
		if (material.is_lighting)
		{
			glEnable(GL_LIGHTING);
		}
		else
		{
			glDisable(GL_LIGHTING);
		}

		// z buffer
		switch(material.z_buffer)
		{
		case COMPARISON_NEVER:
			glDisable(GL_DEPTH_TEST);
			break;
		case COMPARISON_LESSEQUAL:
			glDepthFunc(GL_LEQUAL);
			break;
		case COMPARISON_EQUAL:
			glDepthFunc(GL_EQUAL);
			break;
		case COMPARISON_LESS:
			glDepthFunc(GL_LESS);
			break;
		case COMPARISON_NOTEQUAL:
			glDepthFunc(GL_NOTEQUAL);
			break;
		case COMPARISON_GREATEREQUAL:
			glDepthFunc(GL_GREATER);
			break;
		case COMPARISON_ALWAYS:
			glDepthFunc(GL_ALWAYS);
			break;
		}

		if (material.z_buffer != COMPARISON_NEVER)
		{
			glEnable(GL_DEPTH_TEST);
		}

		// z write
		if(material.is_z_writable)	
		{
			glDepthMask(GL_TRUE);
		}
		else
		{
			glDepthMask(GL_FALSE);
		}

		// Face culling
		if ((material.is_frontface_culling) && (material.is_backface_culling))
		{
			glCullFace(GL_FRONT_AND_BACK);
			glEnable(GL_CULL_FACE);
		}
		else
		{
			if (material.is_backface_culling)
			{
				glCullFace(GL_BACK);
				glEnable(GL_CULL_FACE);
			}
			else
			{
				if(material.is_frontface_culling)
				{
					glCullFace(GL_FRONT);
					glEnable(GL_CULL_FACE);
				}
				else
				{
					glDisable(GL_CULL_FACE);
				}
			}
		}

		// Normalization 
		if (material.is_normalize_normals)
		{
			glEnable(GL_NORMALIZE);
		}
		else
		{
			glDisable(GL_NORMALIZE);
		}

		// Thickness
		glPointSize(material.thickness);
		glLineWidth(material.thickness);
	}

	void SetDefaultRenderStatesLineMode(OPENGL_COLOR line_color = OPENGL_COLOR(0.0f, 0.0f, 0.0f), GLfloat thickness_input = 1.0f)
	{
		GLfloat temp = material_line.thickness;
		material_line.thickness = thickness_input;
		SetRenderStatesByMaterial(material_line);
		material_line.thickness = temp;
		glColor4f(line_color.GetRed(), line_color.GetGreen(), line_color.GetBlue(), line_color.GetAlpha());
	}

	void SetDefaultRenderStates2DMode()
	{
		SetRenderStatesByMaterial(material_2Dmode);
		glClear(GL_DEPTH_BUFFER_BIT);
	}

	void SetDefaultRenderStates3DMode()
	{
		SetRenderStatesByMaterial(material_3Dmode);
	}

	// Need to review after you study OPEN_GL
	void SetAntialising(ANTIALISING_MODE anti_aliasing)
	{
		// Anti aliasing
		if (GLEW_ARB_multisample && (anti_aliasing && ANTIALISING_MULTISAMPLE))
		{
			if (anti_aliasing & ANTIALISING_ALPHA_TO_COVERAGE)
			{
				glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE_ARB);
			}
			else
			{
				glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE_ARB);
			}

			glEnable(GL_MULTISAMPLE_ARB);
			if (GLEW_NV_multisample_filter_hint)
			{
				glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
			}
		}
		else
		{
			glDisable(GL_MULTISAMPLE_ARB);
		}

		if (anti_aliasing & ANTIALISING_POINT_SMOOTH)
		{
			glEnable(GL_POINT_SMOOTH);
		}
		else
		{
			glDisable(GL_POINT_SMOOTH);
		}

		if (anti_aliasing & ANTIALISING_LINE_SMOOTH)
		{
			glEnable(GL_LINE_SMOOTH);
		}
		else
		{
			glDisable(GL_LINE_SMOOTH);
		}
	}

	// Draw mesh
	void DrawSolidTriangles(list<TRIANGLE*>& triangles)
	{
		list<TRIANGLE*>::iterator it;
		for (it = triangles.begin(); it !=triangles.end(); it++)
		{
			glBegin(GL_TRIANGLES);
			{
				#ifdef USE_FLOAT_T
					glNormal3fv((float*) (*it)->vertices[0]->GetNormal());
					glTexCoord2f((*it)->uv[0].x, (*it)->uv[0].y);
					glVertex3fv((float*) (*it)->vertices[0]->GetPosition());
					glNormal3fv((float*) (*it)->vertices[1]->GetNormal());
					glTexCoord2f((*it)->uv[1].x, (*it)->uv[1].y);
					glVertex3fv((float*) (*it)->vertices[1]->GetPosition());
					glNormal3fv((float*) (*it)->vertices[2]->GetNormal());
					glTexCoord2f((*it)->uv[2].x, (*it)->uv[2].y);
					glVertex3fv((float*) (*it)->vertices[2]->GetPosition());
				#else
					glNormal3dv((double*) (*it)->vertices[0]->GetNormal());
					glTexCoord2f((*it)->uv[0].x, (*it)->uv[0].y);
					glVertex3dv((double*) (*it)->vertices[0]->GetNormal());
					glNormal3dv((double*) (*it)->vertices[1]->GetNormal());
					glTexCoord2f((*it)->uv[1].x, (*it)->uv[1].y);
					glVertex3dv((double*) (*it)->vertices[1]->GetNormal());
					glNormal3dv((double*) (*it)->vertices[2]->GetNormal());
					glTexCoord2f((*it)->uv[2].x, (*it)->uv[2].y);
					glVertex3dv((double*) (*it)->vertices[2]->GetNormal());
				#endif
			}
			glEnd();
		}
	}

	void DrawWireTriangles(list<TRIANGLE*>& triangles)
	{
		list<TRIANGLE*>::iterator it;
		for (it = triangles.begin(); it != triangles.end(); it++)
		{
			glBegin(GL_LINES);
			{
				#ifdef USE_FLOAT_T
					glVertex3fv((*it)->vertices[1]->GetPosition());
					glVertex3fv((*it)->vertices[2]->GetPosition());

					glVertex3fv((*it)->vertices[2]->GetPosition());
					glVertex3fv((*it)->vertices[0]->GetPosition());

					glVertex3fv((*it)->vertices[0]->GetPosition());
					glVertex3fv((*it)->vertices[1]->GetPosition());
				#else
					glVertex3dv((*it)->vertices[1]->GetPosition());
					glVertex3dv((*it)->vertices[2]->GetPosition());

					glVertex3dv((*it)->vertices[2]->GetPosition());
					glVertex3dv((*it)->vertices[0]->GetPosition());

					glVertex3dv((*it)->vertices[0]->GetPosition());
					glVertex3dv((*it)->vertices[1]->GetPosition());
				#endif
			}
			glEnd();
		}
	}
		
	// Draw basic shape
	void DrawSolidSphere(GLdouble radius, GLint slices = 40, GLint stacks = 40)
	{
		glutSolidSphere(radius, slices, stacks);
	}

	void DrawWireSphere(GLdouble radius, GLint slices = 40, GLint stacks = 40)
	{
		glutWireSphere(radius, slices, stacks);
	}

	void DrawSolidCube(GLdouble size)
	{
		glutSolidCube(size);
	}

	void DrawWireCube(GLdouble size)
	{
		glutWireCube(size);
	}

	void DrawSolidBox(GLfloat dx, GLfloat dy, GLfloat dz)
	{
		glBegin(GL_QUADS);
			// -ve y plane
			glNormal3f(0.0f, -1.0f, 0.0f);
			glVertex3f(-dx, -dy, dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(dx, -dy, dz);

			// +ve y plane
			glNormal3f(0.0f, 1.0f, 0.0f);
			glVertex3f(-dx, dy, dz);
			glVertex3f(-dx, dy, -dz);
			glVertex3f(dx, dy, -dz);
			glVertex3f(dx, dy, dz);

			// +ve x plane
			glNormal3f(1.0f, 0.0f, 0.0f);
			glVertex3f(dx, dy, dz);
			glVertex3f(dx, dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(dx, -dy, dz);

			// -ve x plane
			glNormal3f(-1.0f, 0.0f, 0.0f);
			glVertex3f(-dx, dy, dz);
			glVertex3f(-dx, dy, -dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(-dx, -dy, dz);
			
			// Top
			glNormal3f(0.0f, 0.0f, 1.0f);
			glVertex3f(-dx, dy, dz);
			glVertex3f(-dx, -dy, dz);
			glVertex3f(dx, -dy, dz);
			glVertex3f(dx, dy, dz);

			// Bottom 
			glNormal3f(0.0f, 0.0f, -1.0f);
			glVertex3f(dx, dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(-dx, dy, -dz);
		glEnd();
	}

	void DrawWireBox(GLfloat dx, GLfloat dy, GLfloat dz)
	{
		// -ve y plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(-dx, -dy, dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(dx, -dy, dz);
		glEnd();

		// +ve y plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(dx, dy, dz);
			glVertex3f(dx, dy, -dz);
			glVertex3f(-dx, dy, -dz);
			glVertex3f(-dx, dy, dz);
		glEnd();

		// +ve x plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(dx, -dy, dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(dx, dy, -dz);
			glVertex3f(dx, dy, dz);
		glEnd();

		// -ve x plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(-dx, dy, dz);
			glVertex3f(-dx, dy, -dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(-dx, -dy, dz);
		glEnd();

		
		// Top
		glBegin(GL_LINE_LOOP);
			glVertex3f(-dx, dy, dz);
			glVertex3f(-dx, -dy, dz);
			glVertex3f(dx, -dy, dz);
			glVertex3f(dx, dy, dz);
		glEnd();

		// Bottom
		glBegin(GL_LINE_LOOP);
			glVertex3f(dx, dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(-dx, dy, -dz);
		glEnd();
	}
	
	void DrawWireBox(GLfloat dx, GLfloat dy, GLfloat dz, GLint y_period)
	{
		// -ve y plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(-dx, -dy, dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(dx, -dy, dz);
		glEnd();

		// +ve y plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(dx, dy + (y_period - 1)*2*dy, dz);
			glVertex3f(dx, dy + (y_period - 1)*2*dy, -dz);
			glVertex3f(-dx, dy + (y_period - 1)*2*dy, -dz);
			glVertex3f(-dx, dy + (y_period - 1)*2*dy, dz);
		glEnd();

		// +ve x plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(dx, -dy, dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(dx, dy + (y_period - 1)*2*dy, -dz);
			glVertex3f(dx, dy + (y_period - 1)*2*dy, dz);
		glEnd();

		// -ve x plane
		glBegin(GL_LINE_LOOP);
			glVertex3f(-dx, dy + (y_period - 1)*2*dy, dz);
			glVertex3f(-dx, dy + (y_period - 1)*2*dy, -dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(-dx, -dy, dz);
		glEnd();

		
		// Top
		glBegin(GL_LINE_LOOP);
			glVertex3f(-dx, dy + (y_period - 1)*2*dy, dz);
			glVertex3f(-dx, -dy, dz);
			glVertex3f(dx, -dy, dz);
			glVertex3f(dx, dy + (y_period - 1)*2*dy, dz);
		glEnd();

		// Bottom
		glBegin(GL_LINE_LOOP);
			glVertex3f(dx, dy + (y_period - 1)*2*dy, -dz);
			glVertex3f(dx, -dy, -dz);
			glVertex3f(-dx, -dy, -dz);
			glVertex3f(-dx, dy + (y_period - 1)*2*dy, -dz);
		glEnd();
	}

	void CircleTable(double** sint, double** cost, const int n)
	{
		int i;

		// Table Size
		const int size = abs(n);

		// Determine the angle between samples
		const double angle = 2*PI / (double)((n == 0) ? 1:n);

		// Allocate memory for n samples
		*sint = (double*) new double[size + 1];
		*cost = (double*) new double[size + 1];

		// Bail out if memory allocation fails, fgError never returns
		if(!(*sint) || !(*cost))
		{
			delete (*sint);
			delete (*cost);
		}

		// Compute cos and sin around the circle
		(*sint)[0] = 0.0;
		(*cost)[0] = 1.0;

		for (i = 1; i < size; i++)
		{
			(*sint)[i] = sin(angle*i);
			(*cost)[i] = cos(angle*i);
		}

		// Last sample is duplicate of the first
		(*sint)[size] = (*sint)[0];
		(*cost)[size] = (*cost)[0];
	}

	void DrawSolidCylinder(GLdouble radius, GLdouble height, GLint slices = 20, GLint stacks = 10)
	{
		int i, j;

		// Step in z and radius as stacks are drawn
		double z0, z1;
		const double zStep = height / ((stacks > 0) ? stacks : 1);

		// Pre-computed circle
		double *sint, *cost;
		CircleTable(&sint, &cost, slices);

		// Cover the base and top
		// Base
		glBegin(GL_TRIANGLE_FAN);
			glNormal3d(0.0, 0.0, -1.0);
			glVertex3d(0.0, 0.0, 0.0);
			for (j = 0; j <= slices; j++)
			{
				glVertex3d(cost[j]*radius, sint[j]*radius, 0.0);
			}
		glEnd();

		glBegin(GL_TRIANGLE_FAN);
			glNormal3d(0.0, 0.0, 1.0);
			glVertex3d(0.0, 0.0, height);
			for (j = 0; j <= slices; j++)
			{
				glVertex3d(cost[j]*radius, sint[j]*radius, height);
			}
		glEnd();

		// Do the stacks
		z0 = 0.0;
		z1 = zStep;

		for (i = 0; i <= stacks; i++)
		{
			if (i == stacks)
			{
				z1 = height;
			}
			
			glBegin(GL_QUAD_STRIP);
			for (j = 0; j <= slices; j++)
			{
				glNormal3d(cost[j], sint[j], 0.0);
				glVertex3d(cost[j]*radius, sint[j]*radius, z0);
				glVertex3d(cost[j]*radius, sint[j]*radius, z1);
			}
			glEnd();

			z0 = z1;
			z1 += zStep;
		}

		// Release sin and cos tables
		delete sint;
		delete cost;
	}

	void DrawWireCylinder(GLdouble radius, GLdouble height, GLint slices = 20, GLint stacks = 10)
	{
		int i, j;

		// Step in z and radius as stacks are drawn
		double z = 0.0;
		const double zStep = height / ((stacks > 0) ? stacks : 1);

		// Pre-computed circle
		double *sint, *cost;
		CircleTable(&sint, &cost, slices);

		// Draw the stacks;
		for (i = 0; i <= stacks; i++)
		{
			if (i == stacks)
			{
				z = height;
			}

			glBegin(GL_LINE_LOOP);
			for (j = 0; j < slices; j++)
			{
				glNormal3d(cost[j], sint[j], 0.0);
				glVertex3d(cost[j]*radius, sint[j]*radius, z);
			}
			glEnd();

			z += zStep;
		}

		// Draw the slices
		glBegin(GL_LINES);
		for (j = 0; j < slices; j++)
		{
			glNormal3d(cost[j], sint[j], 0.0);
			glVertex3d(cost[j]*radius, sint[j]*radius, 0.0);
			glVertex3d(cost[j]*radius, sint[j]*radius, height);
		}

		// Release sin and cos tables
		delete sint;
		delete cost;
	}	

	void DrawSolidCone(GLdouble base, GLdouble height, GLint slices, GLint stacks)
	{
		// glutSolidCone(base, height, slices, stacks);

		const double zStep = height / ((stacks > 0) ? stacks : 1);
		const double rStep = base / ((stacks > 0) ? stacks : 1);

		// Scaling factors for vertex normals
		const double cosn = (height / sqrt(height*height + base*base));
		const double sinn = (base / sqrt(height*height + base*base));

		// Pre-computed circle
		double *sint, *cost;
		CircleTable(&sint, &cost, slices);

		// Cover the circular base with a triangle fan
		double z0 = 0.0;
		double z1 = zStep;

		double r0 = base;
		double r1 = r0 - rStep;

		glBegin(GL_TRIANGLE_FAN);
			glNormal3d(0.0, 0.0, -1.0);
			glVertex3d(0.0, 0.0, z0);
			for (int j = 0; j <= slices; j++)
			{
				glVertex3d(cost[j]*r0, sint[j]*r0, z0);
			}
		glEnd();

		// Cover each stack with a quad strip, except the top stack
		for (int i = 0; i < stacks - 1; i++)
		{
			glBegin(GL_QUAD_STRIP);
			for (int j = 0; j <= slices; j++)
			{
				glNormal3d(cost[j]*cosn, sint[j]*cosn, sinn);
				glVertex3d(cost[j]*r0  , sint[j]*r0  , z0  );
				glVertex3d(cost[j]*r1  , sint[j]*r1  , z1  );
			}
			z0 = z1;
			z1 += zStep;
			r0 = r1;
			r1 -= rStep;
			glEnd();
		}

		// The top stack is covered with individual triangles
		glBegin(GL_TRIANGLES);
			glNormal3d(cost[0]*sinn, sint[0]*sinn, cosn);
			for (int j = 0; j < slices; j++)
			{
				glVertex3d(cost[j + 0]*r0, sint[j + 0]*r0, z0);
				glVertex3d(0, 0, height);

				glNormal3d(cost[j + 1]*sinn, sint[j + 1]*sinn, cosn);
				glVertex3d(cost[j + 1]*r0  , sint[j + 1]*r0  , z0  );
			}
		glEnd();		

		// Release sin and cos tables
		delete sint;
		delete cost;
	}

	void DrawWireCone(GLdouble base, GLdouble height, GLint slices, GLint stacks)
	{
		// glutWireCone(base, height, slices, stacks);

		double z = 0.0;
		double r = base;
		
		const double zStep = height /((stacks > 0) ? stacks : 1);
		const double rStep = base / ((stacks > 0) ? stacks : 1);
		const double cosn = (height / sqrt(height*height + base*base));
		const double sinn = (base / sqrt(height*height + base*base));

		// Pre-computed circle
		double *sint, *cost;
		CircleTable(&sint, &cost, slices);

		// Draw the stacks
		for (int i = 0; i < stacks; i++)
		{
			glBegin(GL_LINE_LOOP);
			for (int j = 0; j < slices; j++)
			{
				glNormal3d(cost[j]*sinn, sint[j]*sinn, cosn);
				glVertex3d(cost[j]*r   , sint[j]*r   , z   );
			}
			glEnd();

			z += zStep;
			r -= rStep;
		}

		// Draw the slices
		r = base;
		glBegin(GL_LINES);
		for (int j = 0; j < slices; j++)
		{
			glNormal3d(cost[j]*sinn, sint[j]*sinn, cosn);
			glVertex3d(cost[j]*r   , sint[j]*r   , 0.0 );
			glVertex3d(0.0		   , 0.0         , height);
		}
		glEnd();

		// Release sin and cos tables
		delete sint;
		delete cost;
	}

	void DrawSolidSqaurePlane(GLfloat width, GLint slices)
	{
		GLfloat fExtent = width / 2.0f;
		GLfloat fStep = width / (float)slices;

		const GLfloat z = -0.4f;
		GLfloat iStrip, iRun;

		for (iStrip = -fExtent; iStrip <= fExtent; iStrip += fStep)
		{
			glBegin(GL_TRIANGLE_STRIP);
				glNormal3f(0.0f, 1.0f, 0.0f);
				for (iRun = fExtent; iRun >= -fExtent; iRun -= fStep)
				{
					glVertex3f(iStrip, iRun, z);
					glVertex3f(iStrip + fStep, iRun, z);
				}
			glEnd();
		}
	}

	void DrawWireSquarePlane(GLfloat width, GLint slices)
	{
		if((slices % 2) != 0)
		{
			slices++;
		}

		GLboolean lighting;
		glGetBooleanv(GL_LIGHTING, &lighting);
		glDisable(GL_LIGHTING);

		GLfloat fExtent = width / 2.0f;
		GLfloat fStep = width / (float)slices;
		GLfloat y = 0.0f;

		glColor3f(0.5, 0.5, 0.6);
		glLineWidth(1.0);
		
		glBegin(GL_LINES);
		for (GLfloat iLine = fStep - fExtent; iLine <= (fExtent - fStep); iLine += fStep)
		{
			glVertex3f(iLine, y, fExtent);
			glVertex3f(iLine, y,-fExtent);
			glVertex3f(fExtent, y, iLine);
			glVertex3f(-fExtent,y, iLine);
		}
		glEnd();

		glColor3f(0.4, 0.4, 0.5);
		glLineWidth(2.0);
		
		glBegin(GL_LINES);
			glVertex3f(fExtent, y, fExtent);
			glVertex3f(fExtent, y,-fExtent);
			glVertex3f(-fExtent,y, fExtent);
			glVertex3f(-fExtent,y,-fExtent);

			glVertex3f(fExtent, y, fExtent);
			glVertex3f(-fExtent,y, fExtent);
			glVertex3f(fExtent,y, -fExtent);
			glVertex3f(-fExtent,y,-fExtent);

			glVertex3f(0, y, fExtent);
			glVertex3f(0, y, -fExtent);
			glVertex3f(fExtent, y, 0);
			glVertex3f(-fExtent, y, 0);
		glEnd();

		glLineWidth(1.0);
		if (lighting)
		{
			glEnable(GL_LIGHTING);
		}
	}

	void DrawSolid2DCircle(GLfloat radius, GLint slices)
	{
		double *sint, *cost;
		CircleTable(&sint, &cost, slices);

		glBegin(GL_TRIANGLE_FAN);
			glNormal3d(0.0, 0.0, -1.0);
			glVertex3d(0.0, 0.0, 0.0);
			for (int i = 0; i <= slices; i++)
			{
				glVertex3d(cost[i]*radius, sint[i]*radius, 0.0);
			}
		glEnd();

		// Release sin and cos tables
		delete sint;
		delete cost;
	}

	void DrawWire2DCircle(GLfloat radius, GLint slices)
	{
		double *sint, *cost;
		CircleTable(&sint, &cost, slices);

		glBegin(GL_LINE_LOOP);
		for (int i = 0; i <= slices; i++)
		{
			glVertex3d(cost[i]*radius, sint[i]*radius, 0.0);
		}
		glEnd();

		// Release sin and cos tables
		delete sint;
		delete cost;
	}

	void DrawMaterialPoint(const ARRAY_3D<MATERIAL_POINT_3D*>& first_particle_array)
	{
		glBegin(GL_POINTS);
		for (int i = 0; i < first_particle_array.ijk_res; i++)
		{
			MATERIAL_POINT_3D* particle = first_particle_array.values[i];

			while (particle != 0)
			{
				#ifdef USE_FLOAT_T
					glVertex3fv(particle->position.values);
				#else
					glVertex3dv(particle->position.values);
				#endif

				particle = particle->next;
			}
		}
		glEnd();
	}

	void DrawCenterAxis(GLfloat size = 1.2f)
	{
		GLfloat color_x[3] = {1.0f, 0.2f, 0.2f};
		GLfloat color_y[3] = {0.2f, 1.0f, 0.2f};
		GLfloat color_z[3] = {0.2f, 0.2f, 1.0f};
		GLfloat color_char[3] = {0.7f, 0.7f, 0.3f};

		// Draw Cylinder
		GLfloat cyn_width = size/100.0f;
		
		glColor3fv(color_z);
		glPushMatrix();
			DrawSolidCylinder(cyn_width, size);
		glPopMatrix();

		glColor3fv(color_x);
		glPushMatrix();
			glRotatef(90.0, 0.0, 1.0, 0.0);
			DrawSolidCylinder(cyn_width, size);
		glPopMatrix();

		glColor3fv(color_y);
		glPushMatrix();
			glRotatef(-90.0, 1.0, 0.0, 0.0);
			DrawSolidCylinder(cyn_width, size);
		glPopMatrix();

		// Draw Cone
		GLfloat con_radius = size/50.0f;
		GLfloat con_height = size/30.0f;

		glColor3fv(color_z);
		glPushMatrix();
			glTranslatef(0.0, 0.0, size);
			DrawSolidCone(con_radius, con_height, 10, 10);
		glPopMatrix();

		glColor3fv(color_x);
		glPushMatrix();
			glTranslatef(size, 0.0, 0.0);
			glRotatef(90.0, 0.0, 1.0, 0.0);
			DrawSolidCone(con_radius, con_height, 10, 10);
		glPopMatrix();

		glColor3fv(color_y);
		glPushMatrix();
			glTranslatef(0.0, size, 0.0);
			glRotatef(-90.0, 1.0, 0.0, 0.0);
			DrawSolidCone(con_radius, con_height, 10, 10);
		glPopMatrix();

		// Draw Character
		glDisable(GL_LIGHTING);
		glColor3fv(color_x);
		glRasterPos3f(size + 0.07f, 0.0f, 0.0f);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'x');

		glColor3fv(color_y);
		glRasterPos3f(0.0f, size + 0.07f, 0.0f);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'y');

		glColor3fv(color_z);
		glRasterPos3f(0.0f, 0.0f, size + 0.07f);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'z');
	}

	void DrawCornerAxis(int size = 100)
	{
		int viewport[4];
		int scissor[4];

		GLfloat color_x[3] = {1.0f, 0.2f, 0.2f};
		GLfloat color_y[3] = {0.2f, 1.0f, 0.2f};
		GLfloat color_z[3] = {0.2f, 0.2f, 1.0f};
		GLfloat color_char[3] = {0.7f, 0.7f, 0.3f};
		GLfloat modelview_mat[16] = {0.0, };

		GLboolean isLighting, isColorMaterial;
		glGetBooleanv(GL_LIGHTING, &isLighting);
		glGetBooleanv(GL_COLOR_MATERIAL, &isColorMaterial);
		glGetFloatv(GL_MODELVIEW_MATRIX, modelview_mat);
		modelview_mat[12] = 0.0f;
		modelview_mat[13] = 0.0f;
		modelview_mat[14] = 0.0f;

		// The viewport and the scissor are changed to fit the lower left corner. Original values are saved.
		glGetIntegerv(GL_VIEWPORT, viewport);
		glGetIntegerv(GL_SCISSOR_BOX, scissor);

		// Axis viewport size, in pixels
		glViewport(0, 0, size, size);
		glScissor(0, 0, size, size);

		// The Z-buffer is cleared to make the axis appear over the original image.
		glLineWidth(2.0);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
			glLoadIdentity();
			glOrtho(-1.1, 1.1, -1.1, 1.1, -1.1, 1.1);

			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
				glLoadIdentity();
				glMultMatrixf(modelview_mat);

				glBegin(GL_LINES);
					glColor3fv(color_x);
					glVertex3f(0.0, 0.0, 0.0);
					glVertex3f(1.0, 0.0, 0.0);

					glColor3fv(color_y);
					glVertex3f(0.0, 0.0, 0.0);
					glVertex3f(0.0, 1.0, 0.0);

					glColor3fv(color_z);
					glVertex3f(0.0, 0.0, 0.0);
					glVertex3f(0.0, 0.0, 1.0);
				glEnd();

				glColor3fv(color_x);
				glRasterPos3f(1.02f, 0.0f, 0.0f);
				glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'x');

				glColor3fv(color_y);
				glRasterPos3f(0.0f, 1.02f, 0.0f);
				glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'y');

				glColor3fv(color_z);
				glRasterPos3f(0.0f, 0.0f, 1.02f);
				glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'z');

				// The viewport and the scissor are restored.
				glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
				glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
	}

	void DrawGradientBackground(const GLfloat* upperColor, const GLfloat* lowerColor)
	{
		GLboolean isDepthTest, isLighting;
		glGetBooleanv(GL_DEPTH_TEST, &isDepthTest);
		glGetBooleanv(GL_LIGHTING, &isLighting);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);

		// glClear(GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glLoadIdentity();

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
				glLoadIdentity();

				glBegin(GL_QUADS);
					// Lower color
					glColor3fv(lowerColor);
					glVertex2f(-1.0, -1.0);
					glVertex2f(1.0, -1.0);
					
					// Upper color
					glColor3fv(upperColor);
					glVertex2f(1.0, 1.0);
					glVertex2f(-1.0, 1.0);
				glEnd();
	
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		if (isDepthTest)
		{
			glEnable(GL_DEPTH_TEST);
		}				
		if (isLighting)
		{
			glEnable(GL_LIGHTING);
		}	
	}

	void DrawGradientBackground(GLfloat upperR, GLfloat upperG, GLfloat upperB, GLfloat lowerR, GLfloat lowerG, GLfloat lowerB)
	{
		GLfloat upper[3] = {upperR, upperG, upperB};
		GLfloat lower[3] = {lowerR, lowerG, lowerB};

		DrawGradientBackground(&upper[0], &lower[0]);
	}


	void DrawGradientBackground(const GLubyte* upperColor, const GLubyte* lowerColor)
	{
		GLboolean isDepthTest, isLighting;
		glGetBooleanv(GL_DEPTH_TEST, &isDepthTest);
		glGetBooleanv(GL_LIGHTING, &isLighting);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);

		// glClear(GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glLoadIdentity();

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
				glLoadIdentity();

				glBegin(GL_QUADS);
					// Lower color
					glColor3ubv(lowerColor);
					glVertex2f(-1.0, -1.0);
					glVertex2f(1.0, -1.0);
					
					// Upper color
					glColor3ubv(upperColor);
					glVertex2f(1.0, 1.0);
					glVertex2f(-1.0, 1.0);
				glEnd();
	
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		if (isDepthTest)
		{
			glEnable(GL_DEPTH_TEST);
		}				
		if (isLighting)
		{
			glEnable(GL_LIGHTING);
		}	
	}

	void DrawGradientBackground(GLubyte upperR, GLubyte upperG, GLubyte upperB, GLubyte lowerR, GLubyte lowerG, GLubyte lowerB)
	{
		GLubyte upper[3] = {upperR, upperG, upperB};
		GLubyte lower[3] = {lowerR, lowerG, lowerB};

		DrawGradientBackground(&upper[0], &lower[0]);
	}

	void DrawSolidTeapot(GLdouble size)
	{
		glutSolidTeapot(size);
	}

	void DrawWireTeapot(GLdouble size)
	{
		glutWireTeapot(size);
	}
	
	CGcontext GetCGContext()
	{
		return cg_context;
	}

	void SetSelectedNameId(int name_id)
	{
		selected_name_id = name_id;
	}

	int GetSelectedNameId()
	{
		return selected_name_id;
	}
};


