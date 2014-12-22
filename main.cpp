#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <ctype.h>

#include "PROJECT_INFO.h"
#include "SIMULATION_MANAGER.h"
#include "OPENGL_GLUT_CAMERA.h"
#include "CAPTURE_MANAGER.h"

SIMULATION_MANAGER simulation;
CAPTURE_MANAGER capture;
OPENGL_GLUT_CAMERA trackball;

int width = 900;
int height = 300;

// Control Position
int last_frame = -1;
bool is_playing = false;
bool is_video_falg = false;
bool is_captur_image = false;
bool is_capture_falg = false;

void set_2d_display_info()
{
	// Set 2D display info
	DISPLAY_INFO info;
	info.camera_mode = "Mouse";
	if (trackball.GetCameraMode() == OPENGL_GLUT_CAMERA::FILE_MODE)
	{
		ostringstream string_stream;
		string_stream << "File " << trackball.GetCameraNum();
		info.camera_mode = string_stream.str();
	}

	info.is_capture_image = is_captur_image;
	info.is_capture_move_at_last_frame = capture.IsAutoCaptureMovieAtLastFrame();
	simulation.Set2DDisplayInfo(info);
}

void reset_config()
{
	last_frame = simulation.GetSimulationWorld()->last_frame;
	is_playing = simulation.GetSimulationWorld()->auto_run;

	// Capture
	capture.Initialize(PROJECT_INFO::GetScriptAbsPath());
	capture.AutoCopyScript();
	is_captur_image = capture.IsAutoCaptureImage();

	trackball.AutoLoadCamera();

	set_2d_display_info();

	is_capture_falg = false;
	is_video_falg = true;
}

void reset_simulation()
{
	simulation.ResetSimulation();
	reset_config();
}

void one_step_simulation()
{
	simulation.OneStepSimulation();
	is_capture_falg = true;
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	trackball.GlutDisplay();
	simulation.Render();

	if (is_captur_image && is_capture_falg)
	{
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		capture.CaptureImage(simulation.GetCurrentFrame(), viewport[2] - viewport[0], viewport[3] - viewport[1]);
		is_capture_falg = false;
	}

	glutSwapBuffers();

	if (is_video_falg && (last_frame > 0) && (last_frame <= simulation.GetCurrentFrame()))
	{
		simulation.FinalizeSimulation();
		is_video_falg = false;
		is_playing = false;						// Stop Simulation
		capture.MakeVideoAtLastFrame();
	}
}

void select_object(int x, int y)
{
	if (0 == simulation.GetOpenGLWorld())
	{
		return;
	}

	trackball.BeginSelect(x, y);
	trackball.GlutDisplay();
	simulation.GetOpenGLWorld()->DrawWithName();
	trackball.EndSelect();
	simulation.GetOpenGLWorld()->SetPickObjectID(trackball.GetSelectedId());
}

void mouse(int button, int state, int x, int y)
{
	trackball.GlutMouse(button, state, x, y);

	// Selection Mode
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && glutGetModifiers() == GLUT_ACTIVE_ALT)
	{
		select_object(x, y);
		glutPostRedisplay();
	}
}

void motion(int x, int y)
{
	trackball.GlutMotion(x, y);
	set_2d_display_info();
}

void idle()
{
	if (is_playing)
	{
		one_step_simulation();
	}

	glutPostRedisplay();
}

void specialKey(int key, int x, int y)
{
	simulation.SpecialKey(key);
	glutPostRedisplay();
}

void key(unsigned char c, int x, int y)
{
	c = tolower(c);

	switch (c)
	{
	case 'q':
		exit(0);
		break;
	// Simulation control
	case 'r':
		reset_simulation();
		break;
	case 's':
		one_step_simulation();
		break;
	case 'p':
		is_playing = !is_playing;
		break;

	// Capture
	case 'c':
		is_captur_image = !is_captur_image;
		break;

	case 'v':
		capture.MakeVideo();
		break;

	// Save / Load camera
	case 'z':
		trackball.SaveCamera();
		set_2d_display_info();
		break;
	case 'x':
		trackball.LoadCamera();
		set_2d_display_info();
		break;

	// For pbrt camera
	case 'f':
		trackball.UpAngle();
		break;
	case 'g':
		trackball.DownAngle();
		break;
	case 'h':
		trackball.SavePBRTCamera();
		break;
	
	// Save / Load simulation data
	case 't':
		simulation.GetSimulationWorld()->Save();
		break;
	case '/':
		simulation.GetSimulationWorld()->Load();
		simulation.GetOpenGLWorld()->Update();
		break;
	}

	int special_key = glutGetModifiers();
	if (special_key == GLUT_ACTIVE_ALT)
	{
		simulation.KeyWithAlt(c);
	}
	else
	{
		simulation.Key(c);
	}

	glutPostRedisplay();
}

void reshape(int w, int h)
{
	width = w;
	height = h;
	trackball.GlutReshape(w, h);
}

void initGL()
{
	glewInit();

	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
}

void exitFunction()
{}

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "Usage: application <script file path>" << endl;
		return 1;
	}

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(width, height);
	glutCreateWindow("Fluid Simulator by BJ");

	glutDisplayFunc(display);
	glutKeyboardFunc(key);
	glutSpecialFunc(specialKey);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutIdleFunc(idle);
	initGL();

	// Key binding info
	cout << "[ShortCut Key]" << endl;
	cout << "    p: Start/Stop Simulation"	<< endl;
	cout << "	 s: One Step Simulation"	<< endl;
	cout << "    c: On/Off Capture Image"   << endl;
	cout << "	 r: Reset Simulation"		<< endl;
	cout << "	 z: Save Camera"			<< endl;
	cout <<	"	 x: Load Camera"			<< endl;
	cout << "	 t: Save Simulation Data"   << endl;
	cout << "	 b: Load Simulation Data"   << endl;
	cout << "	 i: Show/Hide Information"  << endl;
	
	// Simulation
	cout << "Script File: " << argv[1] << endl;
	string script_path = argv[1];

	if (false == simulation.Initialize(script_path))
	{
		return 2;
	}

	capture.Initialize(PROJECT_INFO::GetScriptAbsPath());
	trackball.Initialize(PROJECT_INFO::GetScriptAbsPath(), simulation.GetSimulationWorld()->world_discretization.world_grid);
	reset_config();

	glutMainLoop();

	return 0;
}

