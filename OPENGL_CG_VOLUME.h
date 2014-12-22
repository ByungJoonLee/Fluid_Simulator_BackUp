#pragma once

#include "OPENGL_OBJECT_BASE.h"
#include <Cg/cgGL.h>

enum CG_RGB_MODE
{
	CG_R = 0x0001,
	CG_G = 0x0002,
	CG_B = 0x0004,
};

class OPENGL_CG_VOLUME : public OPENGL_OBJECT_BASE
{
public: // Enumerates
	enum CG_VOLUME_DRAW_TYPE
	{
		CG_VOLUME_DRAW_HIDE = 0,
		CG_VOLUME_DRAW_TRANSPARENT,
		CG_VOLUME_DRAW_OPAQUE,
	};

public: // Essential Data
	CGcontext				cg_context;
	CGprofile				cg_vprofile;
	CGprofile				cg_fprofile;
	CGprogram				raymarch_vprog;
	CGprogram				raymarch_fprog;

	CGparameter				density_param;
	CGparameter				brightness_param;

	GLuint					texture_name;
	float					density;
	float					brightness;

	FIELD_STRUCTURE_3D<T>*	volume_object;
	float					fscale;
	int						rgb_mode;
	bool					is_gen_3d_tex;

	float*					values;

public: // Constructor and Destructor
	OPENGL_CG_VOLUME(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_3D<T>* volume_object_input, float fscale_input = 1.0f, int rgb_mode_input = (CG_R|CG_G|CG_B))
		: OPENGL_OBJECT_BASE(display_name, driver)
		, cg_context(driver->GetCGContext())
		, raymarch_vprog(0)
		, raymarch_fprog(0)
		, density(0.1)
		, brightness(3.0)
		, volume_object(volume_object_input)
		, fscale(fscale_input)
		, rgb_mode(rgb_mode_input)
		, is_gen_3d_tex(false)
		, values(0)
	{
		RegisterDrawType((int) CG_VOLUME_DRAW_HIDE, "HIDE");
		RegisterDrawType((int) CG_VOLUME_DRAW_TRANSPARENT, "DRAW_TRANSPARENT");
		RegisterDrawType((int) CG_VOLUME_DRAW_OPAQUE, "DRAW_OPAQUE");
		SetDrawType((int) CG_VOLUME_DRAW_TRANSPARENT);

		LoadCGScript();
	}

	~OPENGL_CG_VOLUME(void)
	{
		if (raymarch_vprog != 0)
		{
			cgDestroyProgram(raymarch_vprog);
		}
		if (raymarch_fprog != 0)
		{
			cgDestroyProgram(raymarch_fprog);
		}

		if (values)
		{
			delete [] values;
		}
	}

public: // Member Functions
	float GetDensity()
	{
		return density;
	}

	float GetBrightness()
	{
		return brightness;
	}

	void SetDensity(float density_input)
	{
		density = density_input;
	}

	void SetBrightness(float bright_input)
	{
		brightness = bright_input;
	}

	void RenderVolume()
	{
		cgGLBindProgram(raymarch_vprog);
		cgGLEnableProfile(cg_vprofile);

		cgGLBindProgram(raymarch_fprog);
		cgGLEnableProfile(cg_fprofile);

		cgGLSetParameter1f(density_param, density);
		cgGLSetParameter1f(brightness_param, brightness);

		glActiveTextureARB(GL_TEXTURE0_ARB);

		glBindTexture(GL_TEXTURE_3D, texture_name);

		glPushMatrix();
			glTranslatef(GetCenter().x, GetCenter().y, GetCenter().z);
			GetDriver()->DrawSolidBox(GetLength().x, GetLength().y, GetLength().z);
		glPopMatrix();

		cgGLDisableProfile(cg_vprofile);
		cgGLDisableProfile(cg_fprofile);
	}

	void RenderTransparent()
	{
		if (!is_gen_3d_tex)
		{
			Update();
		}

		GetDriver()->SetRenderStatesByMaterial(material);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_COLOR);

		RenderVolume();

		glDisable(GL_BLEND);
	}

	void RenderOpaque()
	{
		if (!is_gen_3d_tex)
		{
			Update();
		}

		GetDriver()->SetRenderStatesByMaterial(material);
		
		RenderVolume();
	}

	void LoadCGScript()
	{
		cg_vprofile = cgGLGetLatestProfile(CG_GL_VERTEX);
		cg_fprofile = cgGLGetLatestProfile(CG_GL_FRAGMENT);

		string resolved_path = "./cg/raymarch.cg";
		raymarch_vprog = cgCreateProgramFromFile(cg_context, CG_SOURCE, resolved_path.c_str(), cg_vprofile, "RayMarchVP", 0);
		cgGLLoadProgram(raymarch_vprog);

		raymarch_fprog = cgCreateProgramFromFile(cg_context, CG_SOURCE, resolved_path.c_str(), cg_fprofile, "RayMarchFP", 0);
		cgGLLoadProgram(raymarch_fprog);

		density_param = cgGetNamedParameter(raymarch_fprog, "density");
		brightness_param = cgGetNamedParameter(raymarch_fprog, "brightness");
	}

	void Generate3DTextureFrom1DArray(const int width, const int height, const int depth, const T* density_array)
	{
		// Convert 1D simulation data to RGB data
		const int num = width*height*depth*4;
		if (!values)
		{
			values = new float[num];
			glGenTextures(1, &texture_name);
		}

		int count_rgb(0), count_density(0);
		while (count_rgb < num)
		{
			float density = (float)density_array[count_density++] / fscale;
			if (rgb_mode & CG_R)	// Red
			{
				values[count_rgb++] = density;
			}
			else
			{
				values[count_rgb++] = 0.0f;
			}

			if (rgb_mode && CG_G)	// Green
			{
				values[count_rgb++] = density;
			}
			else
			{
				values[count_rgb++] = 0.0f;
			}

			if (rgb_mode && CG_B)	// Blue
			{
				values[count_rgb++] = density;
			}
			else
			{
				values[count_rgb++] = 0.0f;
			}
			
			values[count_rgb++] = 0.1f;
		}

		glBindTexture(GL_TEXTURE_3D, texture_name);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLint mode = GL_CLAMP_TO_BORDER;
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, mode);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, mode);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, mode);

		glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA16F_ARB, width, height, depth, 0, GL_RGBA, GL_FLOAT, values);
	}

public: // Virtual Functions
	virtual int NextDrawType()
	{
		switch (GetDrawType())
		{
		case CG_VOLUME_DRAW_HIDE:
			SetDrawType((int) CG_VOLUME_DRAW_TRANSPARENT);
			break;
		case CG_VOLUME_DRAW_TRANSPARENT:
			SetDrawType((int) CG_VOLUME_DRAW_OPAQUE);
			break;
		case CG_VOLUME_DRAW_OPAQUE:
			SetDrawType((int) CG_VOLUME_DRAW_HIDE);
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch (GetDrawType())
		{
		case CG_VOLUME_DRAW_HIDE:
			SetDrawType((int) CG_VOLUME_DRAW_OPAQUE);
			break;
		case CG_VOLUME_DRAW_TRANSPARENT:
			SetDrawType((int) CG_VOLUME_DRAW_HIDE);
			break;
		case CG_VOLUME_DRAW_OPAQUE:
			SetDrawType((int) CG_VOLUME_DRAW_TRANSPARENT);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{
		if (!volume_object)
		{
			return;
		}

		if (GetDrawType() == CG_VOLUME_DRAW_HIDE)
		{
			is_gen_3d_tex = false;
		}
		else
		{
			Generate3DTextureFrom1DArray(volume_object->grid.i_res, volume_object->grid.j_res, volume_object->grid.k_res, volume_object->array_for_this.values);
			is_gen_3d_tex = true;
		}
	}

	virtual void Render()
	{
		switch (GetDrawType())
		{
		case CG_VOLUME_DRAW_HIDE:
			break;
		case CG_VOLUME_DRAW_TRANSPARENT:
			RenderTransparent();
			break;
		case CG_VOLUME_DRAW_OPAQUE:
			RenderOpaque();
			break;
		}
	}

	virtual void RenderWithName()
	{
		Render();
	}
};

