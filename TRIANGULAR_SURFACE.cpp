#include "stdafx.h"
#include "TRIANGULAR_SURFACE.h"

#include "GL/glut.h"
#include "GL/gl.h"

#include <math.h>
#include <algorithm>

#include <string.h>

#include "io.h"
#include <boost/format.hpp>
#include <direct.h>

#define TRAVERSE_VERTICES	vector<VERTEX*>::iterator itr_vertex;	for(itr_vertex = vertices.begin(); itr_vertex != vertices.end(); itr_vertex++)
#define TRAVERSE_EDGES		list<EDGE*>::iterator itr_edge;			for(itr_edge = edges.begin(); itr_edge != edges.end(); itr_edge++)
#define TRAVERSE_TRIANGLES  list<TRIANGLE*>::iterator itr_triangle; for(itr_triangle = triangles.begin(); itr_triangle != triangles.end(); itr_triangle++)

//////////////////////////////////////////////////////////////////////// 
//				 VERTEX Constructors and Destructor					 //
///////////////////////////////////////////////////////////////////////
VERTEX::VERTEX(void)
{
	feature = false;
	is_mc_vertex = false;
	is_boundary = false;
	x[0] = x[1] = x[2] = 0;
	ARRAY_VECTOR3::set<T>(curvature_normal, (T)0);
	ARRAY_VECTOR3::set<T>(deviation, (T)0);
	ARRAY_VECTOR3::set<T>(velocity, (T)0);
	normalizer = (T)0;
	normal_deviation = (T)0;
}

VERTEX::VERTEX(T* x)
{
	feature = false;
	is_mc_vertex = false;
	is_boundary = false;
	SetPosition(x);
	ARRAY_VECTOR3::set<T>(curvature_normal, (T)0);
	ARRAY_VECTOR3::set<T>(deviation, (T)0);
	ARRAY_VECTOR3::set<T>(velocity, (T)0);
	normalizer = (T)0;
	normal_deviation = (T)0;
}

VERTEX::VERTEX(T* x, T* n)
{
	feature = false;
	is_mc_vertex = false;
	is_boundary = false;
	SetPosition(x);
	SetNormal(n);
	ARRAY_VECTOR3::set<T>(curvature_normal, (T)0);
	ARRAY_VECTOR3::set<T>(deviation, (T)0);
	ARRAY_VECTOR3::set<T>(velocity, (T)0);
	normalizer = (T)0;
	normal_deviation = (T)0;
}

VERTEX::~VERTEX(void)
{}

//////////////////////////////////////////////////////////////////////// 
//					VERTEX Member Functions						     //
///////////////////////////////////////////////////////////////////////
void VERTEX::SetPosition(T* x_input)
{
	ARRAY_VECTOR3::set<T>(x, x_input);
}

void VERTEX::SetNormal(T* n_input)
{
	ARRAY_VECTOR3::set<T>(n, n_input);
}

T* VERTEX::GetPosition()
{
	return x;
}

T* VERTEX::GetNormal()
{
	return n;
}

void VERTEX::AddTriangle(TRIANGLE* triangle)
{
	triangles.push_back(triangle);
}

void VERTEX::DelTriangle(TRIANGLE* triangle)
{
	triangles.remove(triangle);
}

void VERTEX::DetermineNormal()
{
	if((int)triangles.size() == 0)
	{
		n[0] = (T)0;
		n[1] = (T)0;
		n[2] = (T)0;
		return;
	}

	ARRAY_VECTOR3::set<T>(n, (T)0);

	TRAVERSE_TRIANGLES
	{
		ARRAY_VECTOR3::add<T>(n, (*itr_triangle)->GetNormal(), n);
	}

	ARRAY_VECTOR3::normalize<T>(VERTEX::n);

	ARRAY_VECTOR3::set<T>(curvature_normal, (T)0);			// Initialize curvature

	normalizer = (T)0;
}

void VERTEX::DrawNormal()
{
	glPushMatrix();
	glTranslatef(x[0], x[1], x[2]);
	glBegin(GL_LINES);
		#ifdef USE_FLOAT_T
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3fv((float*)GetNormal());
		#else
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3dv(GetNormal());
		#endif
	glEnd();
	glPopMatrix();
}

void VERTEX::DrawCurvatureNormal()
{
	glPushMatrix();
	glTranslatef(x[0], x[1], x[2]);
	T n[3];
	ARRAY_VECTOR3::set<T>(n, curvature_normal);
	ARRAY_VECTOR3::mul<T>(n, (T)1);
	glBegin(GL_LINES);
		#ifdef USE_FLOAT_T
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3fv((float*)n);
		#else
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3dv(n);
		#endif
	glEnd();
	glPopMatrix();
}

void VERTEX::DrawNeighborConnectivity()
{
	TRAVERSE_TRIANGLES
	{
		T v_center[3];
		ARRAY_VECTOR3::set<T>(v_center, (*itr_triangle)->vertices[0]->GetPosition());
		ARRAY_VECTOR3::add<T>(v_center, (*itr_triangle)->vertices[1]->GetPosition(), v_center);
		ARRAY_VECTOR3::add<T>(v_center, (*itr_triangle)->vertices[2]->GetPosition(), v_center);
		ARRAY_VECTOR3::div<T>(v_center, (T)3);

		glBegin(GL_LINES);
			#ifdef USE_FLOAT_T
				glVertex3fv((float*)GetPosition());
				glVertex3fv((float*)v_center);
			#else
				glVertex3dv(GetPosition());
				glVertex3dv(v_center);
			#endif
		glEnd();
	}
}

void VERTEX::DrawDeviation()
{
	if((int)triangles.size() == 0)
	{
		return;
	}
	glPushMatrix();
	glTranslatef(x[0], x[1], x[2]);
	glBegin(GL_LINES);
		#ifdef USE_FLOAT_T
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3fv((float*)deviation);
		#else 
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3dv(deviation);
		#endif
	glEnd();
	glPopMatrix();
}

void VERTEX::DrawVelocity()
{
	glPushMatrix();
	glTranslatef(x[0], x[1], x[2]);
	glBegin(GL_LINE);
		#ifdef USE_FLOAT_T
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3fv((float*)velocity);
		#else 
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3dv(velocity);
		#endif
	glEnd();
	glPopMatrix();
}

// After defining the Triangle structure, check this one more
void VERTEX::Replace(VERTEX* v_after)
{
	if(v_after == this)
	{
		cout << "Want to replace to the same vertex?" << endl;
		return;
	}

	cout << "Vertex replacing";
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		for(int i = 0; i < 3; i++)
		{
			if(triangle->vertices[i] == this)
			{
				triangle->vertices[i] = v_after;
				v_after->AddTriangle(triangle);
			}
		}
	}
	cout << "End" << endl;
	v_after->triangles.unique();
}

void VERTEX::DetCurvatureNormal()
{
	T A(0);
	TRAVERSE_TRIANGLES
	{
		A += (*itr_triangle)->area;
	}
	ARRAY_VECTOR3::div<T>(curvature_normal, (T)4*A);
}

//////////////////////////////////////////////////////////////////////// 
//					EDGE Constructor and Destructor  				 //
///////////////////////////////////////////////////////////////////////
EDGE::EDGE(void)
{}

EDGE::EDGE(VERTEX* v0, VERTEX* v1)
{
	AddVertex(v0);
	AddVertex(v1);
}

EDGE::~EDGE(void)
{}

//////////////////////////////////////////////////////////////////////// 
//					   EDGE Member Functions						 //
///////////////////////////////////////////////////////////////////////
void EDGE::AddVertex(VERTEX* vertex)
{
	vertices.push_back(vertex);
}

void EDGE::AddTriangle(TRIANGLE* triangle)
{
	cout << "EDGE::AddTriangle" << endl;
	exit(1);
}

bool EDGE::IsSame(EDGE* edge)
{
	VERTEX *v00, *v01, *v10, *v11;
	vector<VERTEX*>::iterator itr;

	itr = EDGE::vertices.begin();
	v00 = *itr;
	itr++;
	v01 = *itr;

	itr = edge->vertices.begin();
	v10 = *itr;
	itr++;
	v11 = *itr;

	if((v00 == v10) && (v01 == v11)) return true;
	else if((v00 == v11) && (v01 == v10)) return true;
	else return false;
}

void EDGE::Draw()
{
	glBegin(GL_LINES);
		TRAVERSE_VERTICES
		{
			#ifdef USE_FLOAT_T
				glVertex3fv((float*)(*itr_vertex)->GetPosition());
			#else
				glVertex3dv((*itr_vertex)->GetPosition());
			#endif
		}
	glEnd();
}

//////////////////////////////////////////////////////////////////////// 
//	              TRIANGLE Constructor and Destructor  			     //
///////////////////////////////////////////////////////////////////////
TRIANGLE::TRIANGLE(void)
{}

TRIANGLE::TRIANGLE(VERTEX* v0, VERTEX* v1, VERTEX* v2)
{
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;

	edge_vertex[0] = NULL;
	edge_vertex[1] = NULL;
	edge_vertex[2] = NULL;

	is_old = false;
	wrong = false;
}

TRIANGLE::~TRIANGLE(void)
{}

//////////////////////////////////////////////////////////////////////// 
//				      TRIANGLE Member Functions						 //
///////////////////////////////////////////////////////////////////////
void TRIANGLE::DelTriangle(TRIANGLE* triangle)
{
	for(int i = 0; i < 3; i++)
	{
		if(triangles[i] == triangle)
		{
			triangles[i] = NULL;
			return;
		}
	}
}

void TRIANGLE::DetermineNormal()
{
	T l0[3];
	T l1[3];
	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices[2]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l0);
	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices[1]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l1);
	ARRAY_VECTOR3::cross<T>(l1, l0, TRIANGLE::n);
	TRIANGLE::area = ARRAY_VECTOR3::det<T>(TRIANGLE::n)/2.0f;

	ARRAY_VECTOR3::normalize<T>(n);
}

T* TRIANGLE::GetNormal()
{
	return n;
}

void TRIANGLE::SetNormal(T* n_input)
{
	ARRAY_VECTOR3::set<T>(TRIANGLE::n, n_input);
}

void TRIANGLE::DrawNormal()
{
	T p[3];
	ARRAY_VECTOR3::set<T>(p, (T)0);
	ARRAY_VECTOR3::add<T>(p, TRIANGLE::vertices[0]->GetPosition(), p);
	ARRAY_VECTOR3::add<T>(p, TRIANGLE::vertices[1]->GetPosition(), p);
	ARRAY_VECTOR3::add<T>(p, TRIANGLE::vertices[2]->GetPosition(), p);
	ARRAY_VECTOR3::div<T>(p, 3.0f);

	glPushMatrix();
	glTranslatef(p[0], p[1], p[2]);
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3fv((float*)TRIANGLE::GetNormal());
	glEnd();
	glPopMatrix();
}

void TRIANGLE::Draw()
{
	glBegin(GL_TRIANGLES);
	{
		#ifdef USE_FLOAT_T
			glNormal3fv(vertices[0]->GetNormal());
			glTexCoord2f(uv[0].x, uv[0].y);
			glVertex3fv(vertices[0]->GetPosition());
			glNormal3fv(vertices[1]->GetNormal());
			glTexCoord2f(uv[1].x, uv[1].y);
			glVertex3fv(vertices[1]->GetPosition());
			glNormal3fv(vertices[2]->GetNormal());
			glTexCoord2f(uv[2].x, uv[2].y);
			glVertex3fv(vertices[2]->GetPosition());
		#else
			glNormal3dv(vertices[0]->GetNormal());
			glTexCoord2d(uv[0].x, uv[0].y);
			glVertex3dv(vertices[0]->GetPosition());
			glNormal3dv(vertices[1]->GetNormal());
			glTexCoord2d(uv[1].x, uv[1].y);
			glVertex3dv(vertices[1]->GetPosition());
			glNormal3dv(vertices[2]->GetNormal());
			glTexCoord2d(uv[2].x, uv[2].y);
			glVertex3dv(vertices[2]->GetPosition());
		#endif
	}
	glEnd();
}

void TRIANGLE::DrawBack()
{
	glBegin(GL_TRIANGLES);
	{
		#ifdef USE_FLOAT_T
			glNormal3f(-vertices[0]->GetNormal()[0], -vertices[0]->GetNormal()[1], -vertices[0]->GetNormal()[2]);
			glTexCoord2f(uv[0].x, uv[0].y);
			glVertex3fv((float*)vertices[0]->GetPosition());
			glNormal3f(-vertices[2]->GetNormal()[0], -vertices[2]->GetNormal()[1], -vertices[2]->GetNormal()[2]);
			glTexCoord2f(uv[2].x, uv[2].y);
			glVertex3fv((float*)vertices[2]->GetPosition());
			glNormal3f(-vertices[1]->GetNormal()[0], -vertices[1]->GetNormal()[1], -vertices[1]->GetNormal()[2]);
			glTexCoord2f(uv[1].x, uv[1].y);
			glVertex3fv((float*)vertices[1]->GetPosition());
		#else
			glNormal3d(-vertices[0]->GetNormal()[0], -vertices[0]->GetNormal()[1], -vertices[0]->GetNormal()[2]);
			glTexCoord2d(uv[0].x, uv[0].y);
			glVertex3dv((double*)vertices[0]->GetPosition());
			glNormal3d(-vertices[2]->GetNormal()[0], -vertices[2]->GetNormal()[1], -vertices[2]->GetNormal()[2]);
			glTexCoord2d(uv[2].x, uv[2].y);
			glVertex3dv((double*)vertices[2]->GetPosition());
			glNormal3d(-vertices[1]->GetNormal()[0], -vertices[1]->GetNormal()[1], -vertices[1]->GetNormal()[2]);
			glTexCoord2d(uv[1].x, uv[1].y);
			glVertex3dv((double*)vertices[1]->GetPosition());
		#endif
	}
	glEnd();
}

void TRIANGLE::DrawEdges()
{
	glBegin(GL_LINES);
	{
		#ifdef USE_FLOAT_T
			glVertex3fv(TRIANGLE::vertices[1]->GetPosition());
			glVertex3fv(TRIANGLE::vertices[2]->GetPosition());
		#else
			glVertex3dv(TRIANGLE::vertices[1]->GetPosition());
			glVertex3dv(TRIANGLE::vertices[2]->GetPosition());
		#endif
	}
	glEnd();

	glBegin(GL_LINES);
	{
		#ifdef USE_FLOAT_T
			glVertex3fv(TRIANGLE::vertices[2]->GetPosition());
			glVertex3fv(TRIANGLE::vertices[0]->GetPosition());
		#else
			glVertex3dv(TRIANGLE::vertices[2]->GetPosition());
			glVertex3dv(TRIANGLE::vertices[0]->GetPosition());
		#endif
	}
	glEnd();
	
	glBegin(GL_LINES);
	{
		#ifdef USE_FLOAT_T
			glVertex3fv(TRIANGLE::vertices[0]->GetPosition());
			glVertex3fv(TRIANGLE::vertices[1]->GetPosition());
		#else
			glVertex3dv(TRIANGLE::vertices[0]->GetPosition());
			glVertex3dv(TRIANGLE::vertices[1]->GetPosition());
		#endif
	}
	glEnd();
}

void TRIANGLE::DrawCenter()
{
	T center[3];
	ARRAY_VECTOR3::set<T>(center, (T)0);
	ARRAY_VECTOR3::add<T>(center, vertices[0]->GetPosition(), center);
	ARRAY_VECTOR3::add<T>(center, vertices[1]->GetPosition(), center);
	ARRAY_VECTOR3::add<T>(center, vertices[2]->GetPosition(), center);
	ARRAY_VECTOR3::div<T>(center, 3.0f);

	glBegin(GL_POINTS);
		glVertex3fv((float*)center);
	glEnd();
}

void TRIANGLE::DrawNeighborConnectivity()
{
	T center[3];
	ARRAY_VECTOR3::set<T>(center, (T)0);
	ARRAY_VECTOR3::add<T>(center, vertices[0]->GetPosition(), center);
	ARRAY_VECTOR3::add<T>(center, vertices[1]->GetPosition(), center);
	ARRAY_VECTOR3::add<T>(center, vertices[2]->GetPosition(), center);
	ARRAY_VECTOR3::div<T>(center, 3.0f);

	for(int i = 0; i < 3; i++)
	{
		if(TRIANGLE::triangles[i] == NULL)
		{
			cout << "No triangles are connected" << endl;
			continue;
		}
		T neighbor[3];
		
		ARRAY_VECTOR3::set<T>(neighbor, (T)0);
		ARRAY_VECTOR3::add<T>(neighbor, TRIANGLE::triangles[i]->vertices[0]->GetPosition(), neighbor);
		ARRAY_VECTOR3::add<T>(neighbor, TRIANGLE::triangles[i]->vertices[1]->GetPosition(), neighbor);
		ARRAY_VECTOR3::add<T>(neighbor, TRIANGLE::triangles[i]->vertices[2]->GetPosition(), neighbor);
		ARRAY_VECTOR3::div<T>(neighbor, 3.0f);

		ARRAY_VECTOR3::add<T>(center, neighbor, neighbor);
		ARRAY_VECTOR3::div<T>(neighbor, 2.0f);

		glBegin(GL_LINES);
			glVertex3fv((float*)center);
			glVertex3fv((float*)neighbor);
		glEnd();
	}
}

// Need to understand this when you have a time
void TRIANGLE::CorrectCCW()
{
	T l0[3];
	T l1[3];

	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices[2]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l0);
	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices[1]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l1);

	T c[3];
	ARRAY_VECTOR3::cross<T>(l1, l0, c);				// Face normal vector

	TRIANGLE::area = ARRAY_VECTOR3::det(c)/2.0f;	// Calculate triangle area

	T n[3];
	ARRAY_VECTOR3::add<T>(TRIANGLE::vertices[0]->GetNormal(), TRIANGLE::vertices[1]->GetNormal(), n);
	ARRAY_VECTOR3::add<T>(TRIANGLE::vertices[2]->GetNormal(), n, n);

	if((c[0]*n[0] + c[1]*n[1] + c[2]*n[2]) < (T)0)
	{
		VERTEX* temp = TRIANGLE::vertices[1];
		TRIANGLE::vertices[1] = TRIANGLE::vertices[2];
		TRIANGLE::vertices[2] = temp;
	}
}

int TRIANGLE::CountFeatureVertex()
{
	int number = 0;
	for(int i = 0; i < 3; i++)
	{
		if(vertices[i]->feature == true)
		{
			number++;
		}
	}
	return number;
}

void TRIANGLE::FindIntersection(int direction, T* p_input, T* p_intersection, T* uv)
{
	T l0[3];
	T l1[3];

	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices[1]->x, TRIANGLE::vertices[0]->x, l0);
	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices[2]->x, TRIANGLE::vertices[0]->x, l1);

	T l0x = l0[(direction+1)%3];
	T l0y = l0[(direction+2)%3];
	T l1x = l0[(direction+1)%3];
	T l1y = l0[(direction+2)%3];

	T p[2];
	p[0] = p_input[(direction+1)%3] - vertices[0]->x[(direction+1)%3];
	p[0] = p_input[(direction+2)%3] - vertices[0]->x[(direction+2)%3];

	T det = l0x*l1y - l1x*l0y;
	det = 1.0f/det;
	uv[0] = det*(l1y*p[0] - l1x*p[1]);
	uv[1] = det*(-l0y*p[0] + l0x*p[1]);

	ARRAY_VECTOR3::set<T>(p_intersection, vertices[0]->x);

	p_intersection[0] += l0[0]*uv[0] + l1[0]*uv[1];
	p_intersection[1] += l0[1]*uv[0] + l1[1]*uv[1];
	p_intersection[2] += l0[2]*uv[0] + l1[2]*uv[1];
	
	return;
}

bool TRIANGLE::IsInside(VERTEX* v)
{
	if(vertices[0] == v)
	{
		return true;
	}
	if(vertices[1] == v)
	{
		return true;
	}
	if(vertices[2] == v)
	{
		return true;
	}
	return false;
}

void TRIANGLE::Flip()
{}

void TRIANGLE::Flip(VERTEX* v0, VERTEX* v1)
{}

VERTEX* TRIANGLE::FindAnotherVertex(VERTEX* v0, VERTEX* v1)
{
	for(int i = 0; i < 3; i++)
	{
		if(vertices[i] != v0 && vertices[i] != v1)
		{
			return vertices[i];
		}
	}
	
	return NULL;
}

int TRIANGLE::GetNeighborIndex(TRIANGLE* triangle)
{
	for(int i = 0; i < 3; i++)
	{
		if(triangles[i] == triangle)
		{
			return i;
		}
	}
	cout << "TRIANGLE::GetNeighborIndex error!" << endl;
	return 0;
}

int TRIANGLE::GetVertexIndex(VERTEX* v)
{
	for(int i = 0; i < 3; i++)
	{
		if(vertices[i] == v)
		{
			return i;
		}
	}
	cout << "TRIANGLE::GetVertexIndex error!" << endl;
	return -1;
}

void TRIANGLE::ChkNeighborConnectivity()
{
	for(int i = 0; i < 3; i++)
	{
		triangles[i] = NULL;
		{
			list<TRIANGLE*>::iterator itr;
			vertices[(i+1)%3]->DelTriangle(this);
			vertices[(i+2)%3]->DelTriangle(this);
			for(itr = vertices[(i+1)%3]->triangles.begin(); itr != vertices[(i+1)%3]->triangles.end(); itr++)
			{
				if(find(vertices[(i+2)%3]->triangles.begin(), vertices[(i+2)%3]->triangles.end(), *itr) != vertices[(i+2)%3]->triangles.end())
				{
					triangles[i] = *itr;
					break;
				}
			}
			vertices[(i+1)%3]->AddTriangle(this);
			vertices[(i+2)%3]->AddTriangle(this);
		}
	}

	for(int i = 0; i < 3; i++)
	{
		if(triangles[i] == NULL)
		{
			continue;
		}
		if(triangles[i] == triangles[(i+1)%3])
		{
			wrong = true;
			triangles[i]->wrong = true;
			triangles[(i+1)%3]->wrong = true;
		}
		if(triangles[i] == triangles[(i+2)%3])
		{
			wrong = true;
			triangles[i]->wrong = true;
			triangles[(i+2)%3]->wrong = true;
		}
		if(vertices[i] == vertices[(i+1)%3] || vertices[i] == vertices[(i+2)%3])
		{
			cout << "Triangle vertices duplication" << endl;
		}
	}
}

T TRIANGLE::GetOppositeEdgeLength(VERTEX* v)
{
	int i;
	for(i = 0;i < 3; i++)
	{
		if(TRIANGLE::vertices[i] == v)
		{
			break;
		}
	}

	T length[3];
	ARRAY_VECTOR3::set<T>(length, vertices[(i+1)%3]->GetPosition());
	ARRAY_VECTOR3::sub<T>(length, vertices[(i+2)%3]->GetPosition(), length);

	return ARRAY_VECTOR3::det(length);
}

void TRIANGLE::AddLocalCurvatureNormal()
{
	for(int i = 0; i < 3; i++)
	{
		T l0[3], l1[3];
		ARRAY_VECTOR3::set<T>(l0, vertices[(i+1)%3]->GetPosition());
		ARRAY_VECTOR3::sub<T>(l0, vertices[i]->GetPosition(), l0);
		ARRAY_VECTOR3::set<T>(l1, vertices[(i+2)%3]->GetPosition());
		ARRAY_VECTOR3::sub<T>(l1, vertices[i]->GetPosition(), l1);

		T length[3];
		ARRAY_VECTOR3::set<T>(length, vertices[(i+1)%3]->GetPosition());
		ARRAY_VECTOR3::sub<T>(length, vertices[(i+2)%3]->GetPosition(), length);

		T weight = ARRAY_VECTOR3::det(length);

		ARRAY_VECTOR3::mul<T>(l0, weight);
		ARRAY_VECTOR3::mul<T>(l1, weight);
		ARRAY_VECTOR3::add<T>(vertices[i]->curvature_normal, l0, vertices[i]->curvature_normal);
		ARRAY_VECTOR3::add<T>(vertices[i]->curvature_normal, l1, vertices[i]->curvature_normal);
		vertices[i]->normalizer += weight;
	}
}

//////////////////////////////////////////////////////////////////////// 
//	                HOLE Constructor and Destructor  			     //
///////////////////////////////////////////////////////////////////////
HOLE::HOLE(void)
{}

HOLE::~HOLE(void)
{}

//////////////////////////////////////////////////////////////////////// 
//	                       HOLE Member Functions  			         //
///////////////////////////////////////////////////////////////////////
EDGE* HOLE::AddEdge(VERTEX* v0, VERTEX* v1)
{
	EDGE* edge = new EDGE(v0, v1);
	HOLE::edges.push_back(edge);
	return edge;
}

EDGE* HOLE::AddEdge(EDGE* edge)
{
	HOLE::edges.push_back(edge);
	return edge;
}

void HOLE::Draw()
{
	TRAVERSE_EDGES
	{
		(*itr_edge)->Draw();
	}
}

//////////////////////////////////////////////////////////////////////// 
//	          Mesh Manager Constructor and Destructor  			     //
///////////////////////////////////////////////////////////////////////	
TRIANGULAR_SURFACE::TRIANGULAR_SURFACE(void)
{}

TRIANGULAR_SURFACE::~TRIANGULAR_SURFACE(void)
{
	Reset();
}

//////////////////////////////////////////////////////////////////////// 
//				     Mesh Manager Member Functions  			     //
///////////////////////////////////////////////////////////////////////	
EDGE* TRIANGULAR_SURFACE::AddEdge(VERTEX* v0, VERTEX* v1)
{
	EDGE* edge = new EDGE(v0, v1);
	edges.push_back(edge);
	
	return edge;
}

TRIANGLE* TRIANGULAR_SURFACE::AddTriangle(int v0, int v1, int v2)
{
	return TRIANGULAR_SURFACE::AddTriangle(TRIANGULAR_SURFACE::vertices[v0], TRIANGULAR_SURFACE::vertices[v1], TRIANGULAR_SURFACE::vertices[v2]);
}

TRIANGLE* TRIANGULAR_SURFACE::AddTriangle(VERTEX* v0, VERTEX* v1, VERTEX* v2)
{
	VERTEX* v[3] = {v0, v1, v2};
	TRIANGLE* triangle = new TRIANGLE(v[0], v[1], v[2]);
	triangles.push_back(triangle);
	
	v[0]->AddTriangle(triangle);
	v[1]->AddTriangle(triangle);
	v[2]->AddTriangle(triangle);

	return triangle;
}

VERTEX* TRIANGULAR_SURFACE::AddVertex(T* x)
{
	VERTEX* vertex = new VERTEX(x);
	vertices.push_back(vertex);
	
	return vertex;
}

VERTEX* TRIANGULAR_SURFACE::AddVertex(T* x, T* n)
{
	VERTEX* vertex = new VERTEX(x, n);
	vertices.push_back(vertex);

	return vertex;
}

VERTEX* TRIANGULAR_SURFACE::AddVertex(VERTEX* vertex)
{
	vertices.push_back(vertex);

	return vertex;
}

void TRIANGULAR_SURFACE::DelAllTriangles()
{
	list<TRIANGLE*>::iterator itr_triangle = TRIANGULAR_SURFACE::triangles.begin();
	
	while(itr_triangle != TRIANGULAR_SURFACE::triangles.end())
	{
		TRIANGLE* triangle = *itr_triangle;
		for(int i = 0; i < 3; i++)
		{
			triangle->vertices[i]->DelTriangle(triangle);
			if(triangle->triangles[i] != NULL)
			{
				triangle->triangles[i]->DelTriangle(triangle);
				triangle->triangles[i] = NULL;
			}
		}

		delete triangle;

		triangles.erase(itr_triangle);

		itr_triangle = triangles.begin();
	}
}

void TRIANGULAR_SURFACE::DelTriangle(TRIANGLE* triangle)
{
	triangles.remove(triangle);

	for(int i = 0; i < 3; i++)
	{
		triangle->vertices[i]->DelTriangle(triangle);
		if(triangle->triangles[i] != NULL)
		{
			triangle->triangles[i]->DelTriangle(triangle);
		}
	}

	delete triangle;
}

void TRIANGULAR_SURFACE::DetCurvatureNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->AddLocalCurvatureNormal();
	}
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DetCurvatureNormal();
	}
}

void TRIANGULAR_SURFACE::DetVertexNormals()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DetermineNormal();
	}
}

void TRIANGULAR_SURFACE::DetFaceNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DetermineNormal();
	}
}

void TRIANGULAR_SURFACE::AverageDuplexPositionNormal(TRIANGULAR_SURFACE* neighbor, float dx)
{
	T critical_value = dx/4.0;
	for(unsigned int i = 0; i < vertices.size(); i++)
	{
		VERTEX* i_vertex = vertices[i];

		for(unsigned int j = 0; j < neighbor->vertices.size(); j++)
		{
			VERTEX* j_vertex = neighbor->vertices[j];
			float dis[3] = {i_vertex->x[0] - j_vertex->x[0], i_vertex->x[1] - j_vertex->x[1], i_vertex->x[2] - j_vertex->x[2]};
			if(abs(dis[0]) < critical_value && abs(dis[1]) < critical_value && abs(dis[2]) < critical_value)
			{
				float avrNormal[3] = {(i_vertex->n[0] + j_vertex->n[0])/2, (i_vertex->n[1] + j_vertex->n[1])/2, (i_vertex->n[2] + j_vertex->n[2])/2 };
				i_vertex->n[0] = avrNormal[0]; i_vertex->n[1] = avrNormal[1]; i_vertex->n[2] = avrNormal[2];
				j_vertex->n[0] = avrNormal[0]; j_vertex->n[1] = avrNormal[1]; j_vertex->n[2] = avrNormal[2];
			}
		}
	}
}

void TRIANGULAR_SURFACE::DetermineNormalDeviation()
{
	TRAVERSE_VERTICES
	{
		T detdeviation = ARRAY_VECTOR3::det<T>((*itr_vertex)->deviation);
		if(ARRAY_VECTOR3::dot<T>((*itr_vertex)->n, (*itr_vertex)->deviation) < (T)0)
		{
			detdeviation *= -(T)1;
		}

		(*itr_vertex)->normal_deviation = detdeviation;
		(*itr_vertex)->deviation[0] = (*itr_vertex)->normal_deviation*(*itr_vertex)->n[0];
		(*itr_vertex)->deviation[1] = (*itr_vertex)->normal_deviation*(*itr_vertex)->n[1];
		(*itr_vertex)->deviation[2] = (*itr_vertex)->normal_deviation*(*itr_vertex)->n[2];
	}
}

void TRIANGULAR_SURFACE::DetTextureCoordinates(T* xy, T* uv, T s)
{
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		for(int i = 0; i < 3; i++)
		{
			VT& uv = triangle->uv[i];
			VERTEX* vertex = triangle->vertices[i];
			uv.x = (vertex->GetPosition()[0] - xy[0])/s*uv.x;
			uv.y = (vertex->GetPosition()[1] - xy[1])/s*uv.y;
			if(uv.x > (T)1)
			{
				uv.x = (T)1;
			}
			if(uv.y > (T)1)
			{
				uv.y = (T)1;
			}
			if(uv.x < (T)0)
			{
				uv.x = (T)0;
			}
			if(uv.y < (T)0)
			{
				uv.y = (T)0;
			}
		}
	}
}

void TRIANGULAR_SURFACE::DetTextureCoordinates(T x, T y, T u, T v, T s)
{
	T xy[2] = {x, y};
	T uv[2] = {u, v};

	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		for(int i = 0; i < 3; i++)
		{
			VT& uv = triangle->uv[i];
			VERTEX* vertex = triangle->vertices[i];
			uv.x = (vertex->GetPosition()[0] - xy[0])/s*uv.x;
			uv.y = (vertex->GetPosition()[1] - xy[1])/s*uv.y;
			if(uv.x > (T)1)
			{
				uv.x = (T)1;
			}
			if(uv.y > (T)1)
			{
				uv.y = (T)1;
			}
			if(uv.x < (T)0)
			{
				uv.x = (T)0;
			}
			if(uv.y < (T)0)
			{
				uv.y = (T)0;
			}
		}
	}
}

void TRIANGULAR_SURFACE::ChkBoundary()
{
	cout << "# Boundary Checking Started" << endl;

	TRAVERSE_TRIANGLES
	{
		if((*itr_triangle)->triangles[0] == NULL)
		{
			AddEdge((*itr_triangle)->vertices[0], (*itr_triangle)->vertices[1]);
		}
		if((*itr_triangle)->triangles[1] == NULL)
		{
			AddEdge((*itr_triangle)->vertices[1], (*itr_triangle)->vertices[2]);
		}
		if((*itr_triangle)->triangles[2] == NULL)
		{
			AddEdge((*itr_triangle)->vertices[2], (*itr_triangle)->vertices[0]);
		}
	}

	cout << "# Finished" << endl;
	cout << "# Number of Edges = " << edges.size() << endl;
	
	list<EDGE*>::iterator itr_edge = edges.begin();

	while(itr_edge != edges.end())
	{
		cout << "# Number of Edges = " << edges.size() << endl;
		HOLE* hole = new HOLE;
		holes.push_back(hole);
		hole->AddEdge((*itr_edge));
		VERTEX* v0 = (*itr_edge)->vertices[0];
		edges.erase(itr_edge);
		cout << "# Number of Edges = " << edges.size() << endl;
		list<EDGE*>::iterator itr_edge = edges.begin();
		while(itr_edge != edges.end())
		{
			if((*itr_edge)->vertices[0] == v0 || (*itr_edge)->vertices[1] == v0)
			{
				hole->AddEdge((*itr_edge));
				if(v0 == (*itr_edge)->vertices[0])
				{
					v0 = (*itr_edge)->vertices[1];
				}
				else
				{
					v0 = (*itr_edge)->vertices[0];
				}
				edges.erase(itr_edge);
				itr_edge = edges.begin();
			}
			else
			{
				itr_edge++;
			}
		}
		itr_edge = edges.begin();
	}
	cout << "# Finished" << endl;
}

void TRIANGULAR_SURFACE::ChkBoundaryVertices(T* size)
{
	TRAVERSE_VERTICES
	{
		for(int i = 0; i < 3; i++)
		{
			if((*itr_vertex)->GetPosition()[i] <= (T)0 || (*itr_vertex)->GetPosition()[i] >= (T)size[i])
			{
				(*itr_vertex)->is_boundary = true;
			}
		}
	}
}

void TRIANGULAR_SURFACE::ChkTrianglesNeighborConnectivity()
{
	{
		TRAVERSE_TRIANGLES
		{
			TRIANGLE* triangle = *itr_triangle;
			triangle->triangles[0] = NULL;
			triangle->triangles[1] = NULL;
			triangle->triangles[2] = NULL;
		}
	}
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		triangle->ChkNeighborConnectivity();
	}
}

VERTEX* TRIANGULAR_SURFACE::ChkNearestTextureCoordinate(T u, T v)
{
	VERTEX* nearest = NULL;
	T distance = (T)100000000;
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		VERTEX** tri_vertices = triangle->vertices;
		for(int i = 0; i < 3; i++)
		{
			T uv[2] = {triangle->uv[i].x, triangle->uv[i].y};
			T dis = (uv[0] - u)*(uv[0] - u) + (uv[1] - v)*(uv[1] - v);
			if(dis < distance)
			{
				distance = dis;
				nearest = tri_vertices[i];
			}
		}
	}
	return nearest;
}

void TRIANGULAR_SURFACE::DrawEdges()
{
	for_each(triangles.begin(), triangles.end(), mem_fun(&TRIANGLE::DrawEdges));
}

void TRIANGULAR_SURFACE::DrawFaceNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawNormal();
	}
}

void TRIANGULAR_SURFACE::DrawHoles()
{
	list<HOLE*>::iterator itr_hole;
	for(itr_hole = holes.begin(); itr_hole != holes.end(); itr_hole++)
	{
		(*itr_hole)->Draw();
	}
}

void TRIANGULAR_SURFACE::DrawHoles(int index)
{
	list<HOLE*>::iterator itr_hole = holes.begin();
	for(int i = 0; i < index; i++)
	{
		itr_hole++;
	}
	(*itr_hole)->Draw();
}

void TRIANGULAR_SURFACE::DrawTriangles(const bool& draw_front)
{
	if(draw_front)
	{
		TRAVERSE_TRIANGLES
		{
			(*itr_triangle)->Draw();
		}
	}
	else
	{
		TRAVERSE_TRIANGLES
		{
			(*itr_triangle)->DrawBack();
		}
	}
}

void TRIANGULAR_SURFACE::DrawVertices()
{
	glBegin(GL_POINTS);
		for(int i = 0; i < (int)vertices.size(); i++)
		{
			glNormal3fv((float*)vertices[i]->n);
			glVertex3fv((float*)vertices[i]->x);
		}
	glEnd();
}

void TRIANGULAR_SURFACE::DrawVertexNormals()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawNormal();
	}
}

void TRIANGULAR_SURFACE::DrawVertexDeviation()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawDeviation();
	}
}

void TRIANGULAR_SURFACE::DrawTrianglesNeighborConnectivity()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawNeighborConnectivity();
	}
}

void TRIANGULAR_SURFACE::DrawTrianglesCenter()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawCenter();
	}
}

void TRIANGULAR_SURFACE::DrawVerticesNeighborConnectivity()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawNeighborConnectivity();
	}
}

void TRIANGULAR_SURFACE::DrawCurvatureNormal()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawCurvatureNormal();
	}
}

void TRIANGULAR_SURFACE::DrawVertexVelocity()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawVelocity();
	}
}

// Need to study about the filtering
void TRIANGULAR_SURFACE::Filtering(T lambda)
{
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles.size() == 0)
		{
			continue;
		}
		if((*itr_vertex)->is_boundary == true)
		{
			continue;
		}

		list<VERTEX*> v;
		list<VERTEX*>::iterator itr_v;
		list<TRIANGLE*>::iterator itr_triangle;

		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));

		T Li[3] = {(T)0, (T)0, (T)0};
		T E = (T)0;
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			T dist[3];
			T c = (T)1;
			ARRAY_VECTOR3::sub<T>((*itr_v)->GetPosition(), (*itr_vertex)->GetPosition(), dist);
			ARRAY_VECTOR3::det<T>(dist, &c);
			E += c;
			Li[0] += ((*itr_v)->GetPosition()[0] - (*itr_vertex)->GetPosition()[0])/c;
			Li[1] += ((*itr_v)->GetPosition()[1] - (*itr_vertex)->GetPosition()[1])/c;
			Li[2] += ((*itr_v)->GetPosition()[2] - (*itr_vertex)->GetPosition()[2])/c;
		}

		Li[0] = (T)2/E*Li[0];
		Li[1] = (T)2/E*Li[1];
		Li[2] = (T)2/E*Li[2];

		(*itr_vertex)->deviation[0] -= lambda*Li[0];
		(*itr_vertex)->deviation[1] -= lambda*Li[1];
		(*itr_vertex)->deviation[2] -= lambda*Li[2];

		(*itr_vertex)->GetPosition()[0] += lambda*Li[0];
		(*itr_vertex)->GetPosition()[1] += lambda*Li[1];
		(*itr_vertex)->GetPosition()[2] += lambda*Li[2];
	}
}

void TRIANGULAR_SURFACE::Filtering2(T lambda)
{
	vector<VERTEX*>::iterator itr_vertex;
	for(itr_vertex = vertices.begin(); itr_vertex != vertices.end(); itr_vertex++)
	{
		list<VERTEX*> v;
		list<VERTEX*>::iterator itr_v;
		list<TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));
		
		ARRAY_VECTOR3::set<T>((*itr_vertex)->s, (T)0);
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			T dist[3];
			ARRAY_VECTOR3::sub<T>((*itr_v)->GetPosition(), (*itr_vertex)->GetPosition(), dist);
			ARRAY_VECTOR3::add<T>((*itr_vertex)->s, dist, (*itr_vertex)->s);
		}
		ARRAY_VECTOR3::div<T>((*itr_vertex)->s, (T)v.size());
		ARRAY_VECTOR3::det<T>((*itr_vertex)->s, &(*itr_vertex)->stress);
		if(((*itr_vertex)->s[0]*(*itr_vertex)->n[0] + (*itr_vertex)->s[1]*(*itr_vertex)->n[1] + (*itr_vertex)->s[2]*(*itr_vertex)->n[2]) < (T)0)
		{
			(*itr_vertex)->stress *= -(T)1;
		}
	}

	for(itr_vertex = vertices.begin(); itr_vertex != vertices.end(); itr_vertex++)
	{
		if((*itr_vertex)->triangles.size() == 0)
		{
			continue;
		}

		list<VERTEX*> v;
		list<VERTEX*>::iterator itr_v;
		list<TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));

		T s_mean[3] = {(T)0, (T)0, (T)0};
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			s_mean[0] += (*itr_v)->s[0];
			s_mean[1] += (*itr_v)->s[1];
			s_mean[2] += (*itr_v)->s[2];
		}
		s_mean[0] /= (T)v.size();
		s_mean[1] /= (T)v.size();
		s_mean[2] /= (T)v.size();

		T dx[3];
		ARRAY_VECTOR3::sub<T>((*itr_vertex)->s, s_mean, dx);
		(*itr_vertex)->GetPosition()[0] += lambda*dx[0];
		(*itr_vertex)->GetPosition()[1] += lambda*dx[1];
		(*itr_vertex)->GetPosition()[2] += lambda*dx[2];
	}
}

void TRIANGULAR_SURFACE::FilteringTangent(T lambda)
{
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles.size() == 0)
		{
			continue;
		}
		if((*itr_vertex)->is_mc_vertex == true)
		{
			continue;
		}

		list <VERTEX*> v;
		list <VERTEX*>::iterator itr_v;
		list <TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));

		T Li[3] = {(T)0, (T)0, (T)0};
		T E = (T)0;
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			T dist[3];
			T c = (T)1;
			ARRAY_VECTOR3::sub<T>((*itr_v)->GetPosition(), (*itr_vertex)->GetPosition(),dist);
			ARRAY_VECTOR3::det<T>(dist,&c);
			E += c;
			Li[0] += ((*itr_v)->GetPosition()[0] - (*itr_vertex)->GetPosition()[0])/c;
			Li[1] += ((*itr_v)->GetPosition()[1] - (*itr_vertex)->GetPosition()[1])/c;
			Li[2] += ((*itr_v)->GetPosition()[2] - (*itr_vertex)->GetPosition()[2])/c;
		}

		Li[0] = (T)2/E*Li[0];
		Li[1] = (T)2/E*Li[1];
		Li[2] = (T)2/E*Li[2];

		T temp = ARRAY_VECTOR3::dot<T>(Li, (*itr_vertex)->GetNormal());
		T n[3];
		ARRAY_VECTOR3::set<T>(n, (*itr_vertex)->GetNormal());
		ARRAY_VECTOR3::mul<T>(n, temp);
		ARRAY_VECTOR3::sub<T>(Li, n, Li);
		(*itr_vertex)->GetPosition()[0] += lambda*Li[0];
		(*itr_vertex)->GetPosition()[1] += lambda*Li[1];
		(*itr_vertex)->GetPosition()[2] += lambda*Li[2];
	}	
}

void TRIANGULAR_SURFACE::FilteringMinimumVariation(T lambda)
{
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles.size() == 0)
		{
			continue;
		}
		if((*itr_vertex)->is_mc_vertex == true)
		{
			continue;
		}

		list <VERTEX*> v;
		list <VERTEX*>::iterator itr_v;
		list <TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));

		T stress = (T)0;
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			if(ARRAY_VECTOR3::dot<T>((*itr_v)->curvature_normal, (*itr_v)->n) < (T)0)
			{
				stress -= ARRAY_VECTOR3::det<T>((*itr_v)->curvature_normal);
			}
			else
			{
				stress += ARRAY_VECTOR3::det<T>((*itr_v)->curvature_normal);
			}
		}

		stress /= (T)v.size();

		T thisstress = (T)0;
		if(ARRAY_VECTOR3::dot<T>((*itr_vertex)->curvature_normal, (*itr_vertex)->n) < (T)0)
		{
			thisstress -= ARRAY_VECTOR3::det<T>((*itr_vertex)->curvature_normal);
		}
		else
		{
			thisstress += ARRAY_VECTOR3::det<T>((*itr_vertex)->curvature_normal);
		}

		T ds = thisstress - stress;

		T dx[3];
		ARRAY_VECTOR3::set<T>(dx, (*itr_vertex)->curvature_normal);
		ARRAY_VECTOR3::normalize<T>(dx);
		ARRAY_VECTOR3::mul<T>(dx, -ds);
		ARRAY_VECTOR3::mul<T>(dx, lambda);
		
		(*itr_vertex)->GetPosition()[0] += dx[0];
		(*itr_vertex)->GetPosition()[1] += dx[1];
		(*itr_vertex)->GetPosition()[2] += dx[2];
	}
}

void TRIANGULAR_SURFACE::Collapse(TRIANGLE* triangle)
{
	VERTEX* v[3];
	v[0] = triangle->vertices[0];
	v[1] = triangle->vertices[1];
	v[2] = triangle->vertices[2];

	T center[3];
	ARRAY_VECTOR3::set<T>(center, v[0]->x);
	ARRAY_VECTOR3::add<T>(center, v[1]->x, center);
	ARRAY_VECTOR3::add<T>(center, v[2]->x, center);
	ARRAY_VECTOR3::div<T>(center, 3.0f);

	T normal[3];
	ARRAY_VECTOR3::set<T>(normal, v[0]->n);
	ARRAY_VECTOR3::add<T>(normal, v[1]->n, normal);
	ARRAY_VECTOR3::add<T>(normal, v[2]->n, normal);
	ARRAY_VECTOR3::div<T>(normal, 3.0f);

	T deviation[3];
	ARRAY_VECTOR3::set<T>(deviation, v[0]->deviation);
	ARRAY_VECTOR3::add<T>(deviation, v[1]->deviation, deviation);
	ARRAY_VECTOR3::add<T>(deviation, v[2]->deviation, deviation);
	ARRAY_VECTOR3::div<T>(deviation, 3.0f);

	VERTEX* v_after = TRIANGULAR_SURFACE::AddVertex(center, normal);
	ARRAY_VECTOR3::set<T>(v_after->deviation, deviation);

	if(v[0]->is_mc_vertex == true || v[1]->is_mc_vertex == true || v[2]->is_mc_vertex == true)
	{
		v_after->is_mc_vertex = true;
	}

	list<TRIANGLE*>::iterator itr_triangle;
	for(int i = 0; i < 3; i++)
	{
		for(itr_triangle = v[i]->triangles.begin(); itr_triangle != v[i]->triangles.end(); itr_triangle++)
		{
			v_after->AddTriangle((*itr_triangle));
		}
		v[i]->triangles.erase(v[i]->triangles.begin(), v[i]->triangles.end());
	}

	for(itr_triangle = v_after->triangles.begin(); itr_triangle != v_after->triangles.end(); itr_triangle++)
	{
		TRIANGLE* triangle = *itr_triangle;
		for(int i = 0; i < 3; i++)
		{
			if((triangle->vertices[i] == v[0]) || (triangle->vertices[i] == v[1]) || (triangle->vertices[i] == v[2]))
			{
				triangle->vertices[i] = v_after;
			}
		}
	}

	for(itr_triangle = v_after->triangles.begin(); itr_triangle != v_after->triangles.end(); itr_triangle++)
	{
		TRIANGLE* triangle = *itr_triangle;
		int num = 0;
		for(int i = 0; i < 3; i++)
		{
			if(triangle->vertices[i] == v_after)
			{
				num++;
			}
		}

		if(num >= 2)
		{
			v_after->triangles.remove(triangle);
			TRIANGULAR_SURFACE::DelTriangle(triangle);
			itr_triangle = v_after->triangles.begin();
		}
	}
}

void TRIANGULAR_SURFACE::Collapse(VERTEX* v0, VERTEX* v1)
{
	T center[3];
	ARRAY_VECTOR3::set<T>(center, v0->x);
	ARRAY_VECTOR3::add<T>(center, v1->x, center);
	ARRAY_VECTOR3::div<T>(center, 2.0f);

	T normal[3];
	ARRAY_VECTOR3::set<T>(normal, v0->n);
	ARRAY_VECTOR3::add<T>(normal, v1->n, normal);
	ARRAY_VECTOR3::div<T>(normal, 2.0f);

	T deviation[3];
	ARRAY_VECTOR3::set<T>(deviation, v0->deviation);
	ARRAY_VECTOR3::add<T>(deviation, v1->deviation, deviation);
	ARRAY_VECTOR3::div<T>(deviation, 2.0f);

	VERTEX* v_after = TRIANGULAR_SURFACE::AddVertex(center, normal);
	ARRAY_VECTOR3::set<T>(v_after->deviation, deviation);

	if(v0->is_mc_vertex == true || v1->is_mc_vertex == true)
	{
		v_after->is_mc_vertex = true;
	}
	
	list<TRIANGLE*>::iterator itr_triangle;
	for(itr_triangle = v0->triangles.begin(); itr_triangle != v0->triangles.end(); itr_triangle++)
	{
		v_after->AddTriangle((*itr_triangle));
	}
	v0->triangles.erase(v0->triangles.begin(), v0->triangles.end());

	for(itr_triangle = v1->triangles.begin(); itr_triangle != v1->triangles.end(); itr_triangle++)
	{
		v_after->AddTriangle((*itr_triangle));
	}
	v1->triangles.erase(v1->triangles.begin(), v1->triangles.end());
	
	for(itr_triangle = v_after->triangles.begin(); itr_triangle != v_after->triangles.end(); itr_triangle++)
	{
		TRIANGLE* triangle = *itr_triangle;

		for(int i = 0; i < 3; i++)
		{
			if((triangle->vertices[i] == v0) || (triangle->vertices[i] == v1))
			{
				triangle->vertices[i] = v_after;
			}
		}
	}

	for(itr_triangle = v_after->triangles.begin(); itr_triangle != v_after->triangles.end(); itr_triangle++)
	{
		TRIANGLE* triangle = *itr_triangle;
		int num = 0;
		for(int i = 0; i < 3; i++)
		{
			if(triangle->vertices[i] == v_after)
			{
				num++;
			}
		}

		if(num >= 2)
		{
			v_after->triangles.remove(triangle);
			TRIANGULAR_SURFACE::DelTriangle(triangle);
			itr_triangle = v_after->triangles.begin();
		}
	}

	TRIANGULAR_SURFACE::ChkTrianglesNeighborConnectivity();
}

void TRIANGULAR_SURFACE::FilpEdges()
{
	TRAVERSE_TRIANGLES
	{
		if((*itr_triangle)->CountFeatureVertex() != 1)
		{
			continue;
		}
		VERTEX* v_feature = NULL;

		for(int i = 0; i < 3; i++)
		{
			if((*itr_triangle)->vertices[i]->feature == true)
			{
				v_feature = (*itr_triangle)->vertices[i];
				break;
			}
		}
		TRIANGLE* t_neighbor = NULL;
		for(int i = 0; i < 3; i++)
		{
			if((*itr_triangle)->triangles[i] == NULL)
			{
				if((*itr_triangle)->triangles[i]->CountFeatureVertex() == 1)
				{
					if((*itr_triangle)->triangles[i]->vertices[0] != v_feature)
					{
						if((*itr_triangle)->triangles[i]->vertices[1] != v_feature)
						{
							if((*itr_triangle)->triangles[i]->vertices[2] != v_feature)
							{
								t_neighbor = (*itr_triangle)->triangles[i];
								break;
							}
						}
					}
				}
			}
		}
		if(t_neighbor == NULL)
		{
			continue;
		}
		VERTEX* v[4], *temp;
		v[0] = (*itr_triangle)->vertices[0];
		v[1] = (*itr_triangle)->vertices[1];
		v[2] = (*itr_triangle)->vertices[2];

		if(v[1]->feature == true)
		{
			temp = v[1];
			v[1] = v[0];
			v[0] = temp;
		}
		if(v[2]->feature == true)
		{
			temp = v[2];
			v[2] = v[0];
			v[0] = temp;
		}
		for(int j = 0; j < 3; j++)
		{
			if(t_neighbor->vertices[j]->feature == true)
			{
				v[3] = t_neighbor->vertices[j];
				break;
			}
		}
		TRIANGLE* t = (*itr_triangle);
		DelTriangle(t);
		DelTriangle(t_neighbor);
		AddTriangle(v[0], v[1], v[2]);
		AddTriangle(v[0], v[1], v[2]);

		itr_triangle = triangles.begin();
	}
}

void TRIANGULAR_SURFACE::CorrectCCW()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->CorrectCCW();
	}
}

void TRIANGULAR_SURFACE::RemoveSmallTriangles(T threshold)
{
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		if(triangle->triangles[0] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles[0]->triangles[0] == NULL) continue;
			if(triangle->triangles[0]->triangles[1] == NULL) continue;
			if(triangle->triangles[0]->triangles[2] == NULL) continue;
		}
		if(triangle->triangles[1] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles[1]->triangles[0] == NULL) continue;
			if(triangle->triangles[1]->triangles[1] == NULL) continue;
			if(triangle->triangles[1]->triangles[2] == NULL) continue;
		}
		if(triangle->triangles[2] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles[2]->triangles[0] == NULL) continue;
			if(triangle->triangles[2]->triangles[1] == NULL) continue;
			if(triangle->triangles[2]->triangles[2] == NULL) continue;
		}

		T l0[3], l1[3], l2[3];
		ARRAY_VECTOR3::sub<T>(triangle->vertices[0]->x, triangle->vertices[1]->x, l0);
		ARRAY_VECTOR3::sub<T>(triangle->vertices[1]->x, triangle->vertices[2]->x, l1);
		ARRAY_VECTOR3::sub<T>(triangle->vertices[2]->x, triangle->vertices[0]->x, l2);

		T det0, det1, det2;
		ARRAY_VECTOR3::det<T>(l0, &det0);
		ARRAY_VECTOR3::det<T>(l1, &det1);
		ARRAY_VECTOR3::det<T>(l2, &det2);

		if(det0 < threshold && det1 < threshold && det2 < threshold)
		{
			TRIANGULAR_SURFACE::Collapse(triangle);
			itr_triangle = triangles.begin();
		}
	}
}

void TRIANGULAR_SURFACE::RemoveSmallEdges(T threshold)
{
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		bool skip = false;
		for(int i = 0; i < 3; i++)
		{
			if(triangle->triangles[i] == NULL)
			{
				skip = true;
				break;
			}
			else
			{
				for(int j = 0; j < 3; j++)
				{
					if(triangle->triangles[i]->triangles[j] == NULL)
					{
						skip = true;
						break;
					}
				}
			}
			if(skip == true) break;
		}
		if(skip == true) continue;

		for(int i = 0; i < 3; i++)
		{
			T l[3];
			ARRAY_VECTOR3::sub<T>(triangle->vertices[i]->x, triangle->vertices[(i+1)%3]->x, l);
			T det;
			ARRAY_VECTOR3::det<T>(l, &det);
			if(det < threshold)
			{
				if(itr_triangle != triangles.begin())
				{
					itr_triangle--;
					Collapse(triangle->vertices[i], triangle->vertices[(i+1)%3]);
				}
				else
				{
					Collapse(triangle->vertices[i], triangle->vertices[(i+1)%3]);
					itr_triangle = triangles.begin();
				}
			}
			continue;
		}
	}
}

void TRIANGULAR_SURFACE::RemoveLongTriangles(T ratio)
{
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		if(triangle->triangles[0] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles[0]->triangles[0] == NULL) continue;
			if(triangle->triangles[0]->triangles[1] == NULL) continue;
			if(triangle->triangles[0]->triangles[2] == NULL) continue;
		}
		if(triangle->triangles[1] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles[1]->triangles[0] == NULL) continue;
			if(triangle->triangles[1]->triangles[1] == NULL) continue;
			if(triangle->triangles[1]->triangles[2] == NULL) continue;
		}
		if(triangle->triangles[2] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles[2]->triangles[0] == NULL) continue;
			if(triangle->triangles[2]->triangles[1] == NULL) continue;
			if(triangle->triangles[2]->triangles[2] == NULL) continue;
		}

		T l0[3], l1[3], l2[3];
		ARRAY_VECTOR3::sub<T>(triangle->vertices[0]->x, triangle->vertices[1]->x, l0);
		ARRAY_VECTOR3::sub<T>(triangle->vertices[1]->x, triangle->vertices[2]->x, l1);
		ARRAY_VECTOR3::sub<T>(triangle->vertices[2]->x, triangle->vertices[0]->x, l2);

		T det0, det1, det2;
		ARRAY_VECTOR3::det<T>(l0, &det0);
		ARRAY_VECTOR3::det<T>(l1, &det1);
		ARRAY_VECTOR3::det<T>(l2, &det2);

		T max = det0;
		if(det1 > max)
		{
			max = det1;
		}
		if(det2 > max)
		{
			max = det2;
		}
		
		T temp = det0 + det1 + det2 - max;

		T ratio = max/temp;

		if(ratio > (T)0.99)
		{
			Collapse(triangle);
			itr_triangle = triangles.begin();
		}
	}
}

void TRIANGULAR_SURFACE::RemoveVertexDeviation()
{
	TRAVERSE_VERTICES
	{
		for(int i = 0; i < 3; i++)
		{
			(*itr_vertex)->GetPosition()[i] += (*itr_vertex)->deviation[i];
			(*itr_vertex)->deviation[i] = (T)0;
		}
	}
}

void TRIANGULAR_SURFACE::UpdateBoundingBox()
{
	if((int)vertices.size() < 1)
	{
		return;
	}

	bounding_box.Initialize(VT(vertices[0]->x), VT(vertices[0]->x));

	for(int i = 1; i < (int)vertices.size(); i++)
	{
		bounding_box.Extend(VT(vertices[i]->x));
	}
}

void TRIANGULAR_SURFACE::Subdivision()
{
	cout << "Subdivision started" << endl;
	cout << "Triangles are " << (int)triangles.size() << endl;
	cout << "Making new edge vertices" << endl;

	TRAVERSE_TRIANGLES
	{
		static int going0 = 0;
		if(going0 > 100)
		{
			cout << ".";
			going0 = 100;
		}
		else
		{
			going0++;
		}
		(*itr_triangle)->is_old = true;
		for(int d = 0; d < 3; d++)
		{
			if((*itr_triangle)->triangles[0] != NULL && (*itr_triangle)->triangles[1] != NULL && (*itr_triangle)->triangles[2])
			{
				// Butterfly Subdivision
				VERTEX* v[8];
				T x[8][3];
				T x_weight[3];
				T n[3];
				v[1]=(*itr_triangle)->vertices[d];
				v[3]=(*itr_triangle)->vertices[(d+1)%3];
				v[4]=(*itr_triangle)->vertices[(d+2)%3];
				v[0]=(*itr_triangle)->triangles[(d+1)%3]->FindAnotherVertex(v[d],v[(d+2)%3]);
				v[2]=(*itr_triangle)->triangles[(d+2)%3]->FindAnotherVertex(v[d],v[(d+1)%3]);
				v[0]=(*itr_triangle)->triangles[(d+1)%3]->FindAnotherVertex(v[1],v[4]);
				v[2]=(*itr_triangle)->triangles[(d+2)%3]->FindAnotherVertex(v[1],v[3]);
				ARRAY_VECTOR3::set<T>(x[1],v[1]->GetPosition());
				ARRAY_VECTOR3::set<T>(x[3],v[3]->GetPosition());
				ARRAY_VECTOR3::set<T>(x[4],v[4]->GetPosition());
				ARRAY_VECTOR3::set<T>(x[0],v[0]->GetPosition());
				ARRAY_VECTOR3::set<T>(x[2],v[2]->GetPosition());
				ARRAY_VECTOR3::div<T>(x[1],(T)8);
				ARRAY_VECTOR3::div<T>(x[3],(T)2*(T)2);
				ARRAY_VECTOR3::div<T>(x[4],(T)2*(T)2);
				ARRAY_VECTOR3::div<T>(x[0],-(T)16);
				ARRAY_VECTOR3::div<T>(x[2],-(T)16);
				ARRAY_VECTOR3::set<T>(x_weight,(T)0);
				ARRAY_VECTOR3::add<T>(x_weight,x[1],x_weight);
				ARRAY_VECTOR3::add<T>(x_weight,x[3],x_weight);
				ARRAY_VECTOR3::add<T>(x_weight,x[4],x_weight);
				ARRAY_VECTOR3::add<T>(x_weight,x[0],x_weight);
				ARRAY_VECTOR3::add<T>(x_weight,x[2],x_weight);
				ARRAY_VECTOR3::set<T>(n,v[3]->GetNormal());
				ARRAY_VECTOR3::add<T>(n,v[4]->GetNormal(),n);
				ARRAY_VECTOR3::normalize<T>(n);
				if((*itr_triangle)->edge_vertex[d]==NULL)
				{
					(*itr_triangle)->edge_vertex[d]=AddVertex(x_weight,n);
					(*itr_triangle)->triangles[d]->edge_vertex[(*itr_triangle)->triangles[d]->GetNeighborIndex((*itr_triangle))]=(*itr_triangle)->edge_vertex[d];
				}
				else
				{
					ARRAY_VECTOR3::add<T>((*itr_triangle)->edge_vertex[d]->GetPosition(),x_weight,(*itr_triangle)->edge_vertex[d]->GetPosition());
				}
			}
			else
			{
				VERTEX *v[8];
				T x[8][3];
				T x_weight[3];
				T n[3];
//				TRIANGLE *thistriangle=*itr_triangle;
				v[3]=(*itr_triangle)->vertices[(d+1)%3];
				v[4]=(*itr_triangle)->vertices[(d+2)%3];
				ARRAY_VECTOR3::set<T>(x[3],v[3]->GetPosition());
				ARRAY_VECTOR3::set<T>(x[4],v[4]->GetPosition());
				ARRAY_VECTOR3::div<T>(x[3],(T)2);
				ARRAY_VECTOR3::div<T>(x[4],(T)2);
				ARRAY_VECTOR3::set<T>(x_weight,(T)0);
				ARRAY_VECTOR3::add<T>(x_weight,x[3],x_weight);
				ARRAY_VECTOR3::add<T>(x_weight,x[4],x_weight);
				ARRAY_VECTOR3::set<T>(n,v[3]->GetNormal());
				ARRAY_VECTOR3::add<T>(n,v[4]->GetNormal(),n);
				ARRAY_VECTOR3::normalize<T>(n);
				if((*itr_triangle)->edge_vertex[d]==NULL)
				{
					(*itr_triangle)->edge_vertex[d]=AddVertex(x_weight,n);
				}
				else
				{
//					ARRAY_VECTOR3::add<T>((*itr_triangle)->edge_vertex_[d]->GetPosition(), x_weight, (*itr_triangle)->edge_vertex_[d]->GetPosition());
					ARRAY_VECTOR3::set<T>((*itr_triangle)->edge_vertex[d]->GetPosition(), x_weight);
				}
			}
		}
	}
	cout<<"splitting triangles"<<endl;
	list<TRIANGLE*>::iterator itr_triangle2=triangles.begin();
	while(itr_triangle2!=triangles.end())
	{	
		static int going1=0;
		if(going1>10000)
		{
			cout<<".";
			going1=0;
		}
		else
		{
			going1++;
		}

		if((*itr_triangle2)->is_old==false)
		{
			break;
		}
//		if((*itr_triangle2)->edge_vertex_[0]==NULL){T temp=(T)3;}
//		if((*itr_triangle2)->edge_vertex_[1] == NULL){T temp=(T)3;}
//		if((*itr_triangle2)->edge_vertex_[2] == NULL){T temp=(T)3;}

		AddTriangle((*itr_triangle2)->vertices[0],(*itr_triangle2)->edge_vertex[2],(*itr_triangle2)->edge_vertex[1])->is_old=false;
		AddTriangle((*itr_triangle2)->edge_vertex[2],(*itr_triangle2)->edge_vertex[1],(*itr_triangle2)->edge_vertex[0])->is_old=false;
		AddTriangle((*itr_triangle2)->vertices[1],(*itr_triangle2)->edge_vertex[0],(*itr_triangle2)->edge_vertex[2])->is_old=false;
		AddTriangle((*itr_triangle2)->vertices[2],(*itr_triangle2)->edge_vertex[0],(*itr_triangle2)->edge_vertex[1])->is_old=false;		

		TRIANGLE *triangle=*itr_triangle2;
		DelTriangle(triangle);
		itr_triangle2=triangles.begin();
	}
	cout<<"Triangles are "<<triangles.size()<<endl;
	cout<<"Subdivision ended"<<endl;
}

void TRIANGULAR_SURFACE::Scale(const T& scale)
{
	for(int i = 0; i < (int)vertices.size(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			vertices[i]->x[j] *= scale;
		}
	}
	bounding_box.Scale(scale);
}

void TRIANGULAR_SURFACE::Translate(T x, T y, T z)
{
	for(int i = 0; i < (int)vertices.size(); i++)
	{
		vertices[i]->x[0] += x;
		vertices[i]->x[1] += y;
		vertices[i]->x[2] += z;
	}
}

// Need to consider this one more time
void TRIANGULAR_SURFACE::FillHoles()
{
	int number_of_holes_filled = 0;

	TRAVERSE_TRIANGLES
	{
		VERTEX* v[4] = {NULL, NULL, NULL, NULL};
		for(int i = 0; i < 3; i++)
		{
			if((*itr_triangle)->triangles[i] == NULL)
			{
				if((*itr_triangle)->vertices[(i+1)%3]->is_boundary == true && (*itr_triangle)->vertices[(i+2)%3]->is_boundary == true) 
				{
					continue;
				}
				v[0] = (*itr_triangle)->vertices[(i+1)%3];
				v[1] = (*itr_triangle)->vertices[(i+2)%3];
				break;
			}
			if(v[0] != NULL && v[1] != NULL)
			{
				break;
			}
		}
		if(v[0] == NULL && v[1] == NULL)
		{
			continue;
		}

		list<TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = v[0]->triangles.begin(); itr_triangle != v[0]->triangles.end(); itr_triangle++)
		{
			for(int i = 0; i < 3; i++)
			{
				if((*itr_triangle)->triangles[i] == NULL)
				{
					if((*itr_triangle)->vertices[(i+1)%3]->is_boundary == true && (*itr_triangle)->vertices[(i+2)%3]->is_boundary == true)
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[0] && (*itr_triangle)->vertices[(i+2)%3] == v[1])
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[1] && (*itr_triangle)->vertices[(i+2)%3] == v[0])
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[0])
					{
						v[2] = (*itr_triangle)->vertices[(i+2)%3];
						break;
					}
					if((*itr_triangle)->vertices[(i+2)%3] == v[0])
					{
						v[2] = (*itr_triangle)->vertices[(i+1)%3];
						break;
					}
					if(v[2] != NULL)
					{
						break;
					}
				}
			}
			if(v[2] != NULL)
			{
				break;
			}
		}
		for(itr_triangle = v[1]->triangles.begin(); itr_triangle != v[1]->triangles.end(); itr_triangle++)
		{
			for(int i = 0; i < 3; i++)
			{
				if((*itr_triangle)->triangles[i] == NULL)
				{
					if((*itr_triangle)->vertices[(i+1)%3]->is_boundary == true && (*itr_triangle)->vertices[(i+2)%3]->is_boundary == true)
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[1] || (*itr_triangle)->vertices[(i+2)%3] == v[0])
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[0] || (*itr_triangle)->vertices[(i+2)%3] == v[1])
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[1] || (*itr_triangle)->vertices[(i+2)%3] == v[2])
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[2] || (*itr_triangle)->vertices[(i+2)%3] == v[1])
					{
						continue;
					}
					if((*itr_triangle)->vertices[(i+1)%3] == v[1])
					{
						v[3] = (*itr_triangle)->vertices[(i+2)%3];
						break;
					}
					else
					{
						v[3] = (*itr_triangle)->vertices[(i+1)%3];
						break;
					}
					if(v[3] != NULL)
					{
						break;
					}
				}
			}
		}
		if(v[2] != NULL)
		{
			AddTriangle(v[0], v[1], v[2])->CorrectCCW();
		}
		if(v[3] != NULL)
		{
			AddTriangle(v[0], v[2], v[3])->CorrectCCW();
		}
		if(v[2] != NULL || v[3] != NULL)
		{
			number_of_holes_filled++;
			ChkTrianglesNeighborConnectivity();
			cout << "." << flush;
		}
	}
	cout << "Total holes filled " << number_of_holes_filled << endl;
}

void TRIANGULAR_SURFACE::RemoveNonManifold()
{
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		if(triangle->vertices[0]->is_boundary == true)
		{
			continue;
		}
		if(triangle->vertices[1]->is_boundary == true)
		{
			continue;
		}
		if(triangle->vertices[2]->is_boundary == true)
		{
			continue;
		}
		if(triangle->triangles[0] == NULL)
		{
			DelTriangle(triangle);
			itr_triangle = triangles.begin();
			continue;
		}
		if(triangle->triangles[1] == NULL)
		{
			DelTriangle(triangle);
			itr_triangle = triangles.begin();
			continue;
		}
		if(triangle->triangles[2] == NULL)
		{
			DelTriangle(triangle);
			itr_triangle = triangles.begin();
			continue;
		}
	}
}

void TRIANGULAR_SURFACE::Reset()
{
	TRAVERSE_VERTICES
	{
		DELETE_POINTER(*itr_vertex);
	}

	TRAVERSE_EDGES
	{
		DELETE_POINTER(*itr_edge);
	}

	TRAVERSE_TRIANGLES
	{
		DELETE_POINTER(*itr_triangle);
	}

	vertices.clear();
	edges.clear();
	triangles.clear();
}

void TRIANGULAR_SURFACE::MinMaxPosition(VT& min, VT& max)
{
	T min_compare_x = vertices[0]->x[0];
	T min_compare_y = vertices[0]->x[1];
	T min_compare_z = vertices[0]->x[2];

	T max_compare_x = vertices[0]->x[0];
	T max_compare_y = vertices[0]->x[1];
	T max_compare_z = vertices[0]->x[2];

	for(unsigned int i = 1; i < vertices.size(); i++)
	{
		VT position(vertices[i]->x);

		min_compare_x = MIN(position.x, min_compare_x);
		min_compare_y = MIN(position.y, min_compare_y);
		min_compare_z = MIN(position.z, min_compare_z);

		max_compare_x = MAX(position.x, max_compare_x);
		max_compare_y = MAX(position.y, max_compare_y);
		max_compare_z = MAX(position.z, max_compare_z);
	}

	min = VT(min_compare_x, min_compare_y, min_compare_z);
	max = VT(max_compare_x, max_compare_y, max_compare_z);
}

void TRIANGULAR_SURFACE::Read(const char* filename)
{
	ifstream ist(filename, ios::binary);
	if(ist == NULL)
	{
		cout << "Failed to open file: " << filename << endl;
	}

	int number_of_vertices = 0;
	int number_of_triangles = 0;

	// Read vertex positions and normals
	ist.read((char*)&number_of_vertices, sizeof(int));
	cout << "Reading vertices " << number_of_vertices << endl;

	for(int i = 0; i < number_of_vertices; i++)
	{
		T position[3], normal[3], velocity[3];
		ist.read((char*)&position[0], sizeof(T));
		ist.read((char*)&position[1], sizeof(T));
		ist.read((char*)&position[2], sizeof(T));
		ist.read((char*)&normal[0], sizeof(T));
		ist.read((char*)&normal[1], sizeof(T));
		ist.read((char*)&normal[2], sizeof(T));
		ist.read((char*)&velocity[0], sizeof(T));
		ist.read((char*)&velocity[1], sizeof(T));
		ist.read((char*)&velocity[2], sizeof(T));
		AddVertex(position, normal);
	}

	// Read Triangles
	ist.read((char*)number_of_triangles, sizeof(int));
	cout << "Reading triangles " << number_of_triangles << endl;
	for(int i = 0; i < number_of_triangles; i++)
	{
		int index[3];
		ist.read((char*)&index[0], sizeof(int));
		ist.read((char*)&index[1], sizeof(int));
		ist.read((char*)&index[2], sizeof(int));
		AddTriangle(index[0], index[1], index[2]);
	}

	ist.close();
}

void TRIANGULAR_SURFACE::ReadOBJ(const char* filename, const T angle, const VT axis, const VT scale, const VT translation)
{
	ifstream file(filename);

	if(file == 0)
	{
		cout << filename << " does not exist. Program terminated." << endl;
		exit(-1);
	}

	vector<VT> normals_temp;
	vector<VT> texture_coordinates_temp;

	char c[255];

	bool flag_v(false);
	bool flag_vt(false);
	bool flag_vn(false);
	bool flag_f(false);

	bool is_first_vertex(true);

	QUATERNION quaternion(-angle*(T)0.0174532925, axis);

	while(1)
	{
		file >> c;

		if(file.eof() != 0)
		{
			break;
		}

		if(strcmp(c, "#") == 0)
		{
			file.getline(c, 255);		// Comments (less than 255 characters)
		}
		else if(strcmp(c, "v") == 0)	// Vertices
		{
			if(flag_v == false)
			{
				flag_v = true;
			}
			VT vertex_pos;
			file >> vertex_pos.x >> vertex_pos.y >> vertex_pos.z;

			if(angle != (T)0)
			{
				vertex_pos = quaternion.Rotate(vertex_pos);
			}

			vertex_pos = vertex_pos + translation;

			vertex_pos.x *= scale.x;
			vertex_pos.y *= scale.y;
			vertex_pos.z *= scale.z;

			VERTEX* vertex_temp = new VERTEX(vertex_pos.values);
			vertex_temp->is_mc_vertex = true;
			vertices.push_back(vertex_temp);

			if(is_first_vertex == true)
			{
				is_first_vertex =false;
				bounding_box.Initialize(vertex_pos, vertex_pos);
			}
			else
			{
				bounding_box.Extend(vertex_pos);
			}
		}
		else if(strcmp(c,"vt") == 0) // texture coordinate
		{	
			if(flag_vt == false)
			{
				flag_vt = true;
			}

			VT texure_vt;
			file >> texure_vt.x >> texure_vt.y;

			texture_coordinates_temp.push_back(texure_vt);			
		}
		else if(strcmp(c,"vn") == 0) // vertex normal
		{	
			if(flag_vn==false)
			{
				flag_vn=true;
			}

			VT vertex_n;
			if(flag_vn)
			{
				file >> vertex_n.x >> vertex_n.y >> vertex_n.z;
			}

			if(angle!=(T)0)
			{
				vertex_n = quaternion.Rotate(vertex_n);
			}

			normals_temp.push_back(vertex_n);
		}
		else if(strcmp(c,"f") == 0)
		{
			if(flag_f == false)
			{
				flag_f = true;
			}

			int v[3],vt[3],vn[3];
			if(flag_vt == true && flag_vn == true)
			{
				for(int i=0;i<3;i++)
				{
					file>>v[i];file.get(c,2);
					file>>vt[i];file.get(c,2);
					file>>vn[i];

					v[i]--;
					vt[i]--;
					vn[i]--;
				}
			}
			else if(flag_vt==false && flag_vn==true)
			{
				for(int i=0;i<3;i++)
				{
					file>>v[i];file.get(c,2);file.get(c,2);
					file>>vn[i];
					v[i]--;
					vn[i]--;
				}
			}
			else if(flag_vt==false && flag_vn==false)
			{
				for(int i=0;i<3;i++)
				{
					file>>v[i];					
					v[i]--;
				}
			}			

			// add triangle
			TRIANGLE* triangle=TRIANGULAR_SURFACE::AddTriangle(vertices[v[0]], vertices[v[1]], vertices[v[2]]);
			triangle->vertex_indices.i = v[0];
			triangle->vertex_indices.j = v[1];
			triangle->vertex_indices.k = v[2];

			if(flag_vt==true)
			{
				// set texture coordinate for vertex
				triangle->uv[0].x = texture_coordinates_temp[vt[0]].x; triangle->uv[0].y = texture_coordinates_temp[vt[0]].y;
				triangle->uv[1].x = texture_coordinates_temp[vt[1]].x; triangle->uv[1].y = texture_coordinates_temp[vt[1]].y;
				triangle->uv[2].x = texture_coordinates_temp[vt[2]].x; triangle->uv[2].y = texture_coordinates_temp[vt[2]].y; 
			}

			if(flag_vn==true)
			{
				T n[3];
				ARRAY_VECTOR3::set<T>(n,(T)0);
				ARRAY_VECTOR3::add<T>(n,normals_temp[vn[0]].values,n);
				ARRAY_VECTOR3::add<T>(n,normals_temp[vn[1]].values,n);
				ARRAY_VECTOR3::add<T>(n,normals_temp[vn[2]].values,n);
				ARRAY_VECTOR3::normalize<T>(n);
				triangle->SetNormal(n);
			}
			else
			{
				triangle->DetermineNormal();
			}
		}
	}

	TRIANGULAR_SURFACE::DetVertexNormals();

	file.clear();
	file.close();
}

void TRIANGULAR_SURFACE::ReadSMF(const char* filename)
{
	cout << "# Reading " << filename << endl;
	
	ifstream file(filename);
	char c[255];
	int frame = 0;
	int flag = 0;
	
	while(1)
	{
		frame++;
		if(frame > 20000)
		{
			cout << ".";
			frame = 0;
		}
		file >> c;
		if(file.eof() != 0)
		{
			break;
		}
		if(strcmp(c, "#") == 0)
		{
			file.getline(c, 255);
		}
		else if(strcmp(c, "v") == 0)
		{
			if(flag == 0)
			{
				cout << "v";
				flag++;
			}
			T x[3];
			file >> x[0] >> x[1] >> x[2];
			TRIANGULAR_SURFACE::AddVertex(x);
		}
		else if(strcmp(c, "f"))
		{
			if(flag == 1)
			{
				cout << "f";
				flag++;
			}
			int v[3];
			file >> v[0] >> v[1] >> v[2];
			v[0]--;
			v[1]--;
			v[2]--;
			TRIANGLE* triangle = AddTriangle(vertices[v[0]], vertices[v[1]], vertices[v[2]]);
			triangle->DetermineNormal();
		}
	}

	DetVertexNormals();
	file.close();

	cout << "# Finished" << endl;
}

void TRIANGULAR_SURFACE::Write(const char* filename)
{
	ofstream ost(filename, ios::binary);
	if(!ost.is_open())
	{
		cout << "Failed to write " << filename << endl;
	}
	int number_of_vertices = (int)vertices.size();
	int number_of_triangles = (int)triangles.size();

	// Write vertex positions and normals
	cout << "Writing vertices " << number_of_vertices << endl;
	ost.write((char*)number_of_vertices, sizeof(int));
	for(int i = 0; i < number_of_vertices; i++)
	{
		VERTEX* v = vertices[i];
		v->index = i;
		ost.write((char*)&v->x[0], sizeof(T));
		ost.write((char*)&v->x[1], sizeof(T));
		ost.write((char*)&v->x[2], sizeof(T));
		ost.write((char*)&v->n[0], sizeof(T));
		ost.write((char*)&v->n[1], sizeof(T));
		ost.write((char*)&v->n[2], sizeof(T));
		ost.write((char*)&v->velocity[0], sizeof(T));
		ost.write((char*)&v->velocity[1], sizeof(T));
		ost.write((char*)&v->velocity[2], sizeof(T));
	}

	// Write Triangles
	cout << "Writing Triangles " << number_of_triangles << endl;
	ost.write((char*)&number_of_triangles, sizeof(int));
	list<TRIANGLE*>::iterator itr_triangle;
	for(itr_triangle = triangles.begin(); itr_triangle != triangles.end(); itr_triangle++)
	{
		ost.write((char*)&(*itr_triangle)->vertices[0]->index, sizeof(int));
		ost.write((char*)&(*itr_triangle)->vertices[1]->index, sizeof(int));
		ost.write((char*)&(*itr_triangle)->vertices[2]->index, sizeof(int));
	}

	ost.close();
}

void TRIANGULAR_SURFACE::WriteOBJ(const char* filename)
{
	cout << "# Writing " << filename << endl;
	ofstream file(filename);

	int index;
	
	// Vertex
	for(index = 1; index <= (int)vertices.size(); index++)
	{
		TRIANGULAR_SURFACE::vertices[index - 1]->index = index;
		file << "v " << vertices[index - 1]->x[0] << " " << vertices[index - 1]->x[1] << " " << vertices[index - 1]->x[2] << endl;
	}
	
	// vt
	for(index = 1; index <= (int)vertices.size(); index++)
	{
		file << "vt " << "1" << " " << "0" << endl;
	}
	
	// Vertex normal
	for(index = 1; index <= (int)vertices.size(); index++)
	{
		file << "vn " << vertices[index - 1]->n[0] << " " << vertices[index - 1]->n[1] << " " << vertices[index - 1]->n[2] << endl;
	}
	
	// Scaled velocity
	for(index = 1; index <= (int)vertices.size(); index++)
	{
		T scale = (T)0.1;
		file << "vv " << vertices[index - 1]->velocity[0]*scale << " " << vertices[index - 1]->velocity[1]*scale << " " << vertices[index - 1]->velocity[2]*scale << endl;
	}

	TRAVERSE_TRIANGLES
	{
		int index[3] = {(*itr_triangle)->vertices[0]->index, (*itr_triangle)->vertices[1]->index, (*itr_triangle)->vertices[2]->index};
		file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << endl;
	}

	file.close();
}

void TRIANGULAR_SURFACE::WriteOBJ(const char* filename, VT& position, QUATERNION& quat)
{
	cout << "# Writing " << filename << endl;
	ofstream file(filename);

	int index;

	VT vertex_position;
	VT vertex_normal;
	VT vertex_velocity;

	for(index = 1; index <= (int)vertices.size(); index++)
	{
		TRIANGULAR_SURFACE::vertices[index - 1]->index = index;

		VT x(vertices[index - 1]->x);
		vertex_position = quat.Rotate(x);
		vertex_position = vertex_position+position;

		file << "v " << vertex_position.x << " " << vertex_position.y << " " << vertex_position.z << endl;
	}

	for(index = 1; index <= (int)vertices.size(); index++)
	{
		file << "vt " << "1" << " " << "0" << endl;
	}

	for(index = 1; index <= (int)vertices.size(); index++)
	{
		VT x(vertices[index - 1]->n);
		vertex_normal = quat.Rotate(x);
		
		file << "vn " << vertex_normal.x << " " << vertex_normal.y << " " << vertex_normal.z << endl;
	}

	for(index = 1; index <= (int)vertices.size(); index++)
	{
		T scale = (T)0.1;
		VT x(vertices[index - 1]->velocity);
		vertex_velocity = quat.Rotate(x);

		file << "vv " << vertex_velocity.x*scale << " " << vertex_velocity.y*scale << " " << vertex_velocity.z*scale << endl;
	}

	TRAVERSE_TRIANGLES
	{
		int index[3] = {(*itr_triangle)->vertices[0]->index, (*itr_triangle)->vertices[1]->index, (*itr_triangle)->vertices[2]->index};
		file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << endl;
	}

	file.close();
}

void TRIANGULAR_SURFACE::WritePOLY(const char* filename)
{
	cout << "# Writing " << filename << endl;
	ofstream file(filename);

	// # Part 1 - Node List
	file << "# Part 1 - Node List" << endl;
	int node_count = (int)vertices.size();
	int dimension = 3;
	int attribute = 0;
	int boundary_marker = 0;

	file << node_count << " " << dimension << " " << attribute << " " << boundary_marker << endl;

	for(int i = 0; i < node_count; i++)
	{
		file << (i+1) << " " << vertices[i]->x[0] << " " << vertices[i]->x[1] << " " << vertices[i]->x[2] << endl;
	}

	// # Part 2 - Facet List
	file << "# Part 2 - Facet List" << endl;
	int facet_count = (int)triangles.size();

	file << facet_count << " " << boundary_marker << endl;

	TRAVERSE_TRIANGLES
	{
		file << 1 << " " << 0 << " " << 1 << endl;
		file << 3 << " ";

		int index[3] = {(*itr_triangle)->vertex_indices.i, (*itr_triangle)->vertex_indices.j, (*itr_triangle)->vertex_indices.k};
		file << index[0] + 1 << " " << index[1] + 1 << " " << index[2] + 1 << endl;
	}

	// # Part 3 - Hole List
	file << "# Part 3 - Hole List" << endl;
	file << 0 << endl;

	// # Part 4 - Region List
	file << "# Part 4 - Region List" << endl;
	file << 0 << endl;

	file.close();
}

void TRIANGULAR_SURFACE::WriteOBJ(ARRAY<TRIANGULAR_SURFACE*>& arr, const char* filename)
{
	cout << "# Writing " << filename << endl;
	ofstream file(filename);

	ARRAY<int> vertex_count(arr.num_elements);

	int vertex_counter = 0;

	for(int i = 0; i < arr.num_elements; i++)
	{
		TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
		int vertices_size = (int)triangular_surface->vertices.size();
		int index;

		for(index = 1; index <= vertices_size; index++)
		{
			file << "v " << triangular_surface->vertices[index - 1]->x[0] << " " << triangular_surface->vertices[index - 1]->x[1] << " " << triangular_surface->vertices[index - 1]->x[2] << endl;
		}
	}

	for(int i = 0; i < arr.num_elements; i++)
	{
		TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
		int vertices_size = (int)triangular_surface->vertices.size();
		int index;

		for(index = 1; index <= vertices_size; index++)
		{
			file << "vt " << "1" << " " << "0" << endl;
		}
	}

	for(int i = 0; i < arr.num_elements; i++)
	{
		TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
		int vertices_size = (int)triangular_surface->vertices.size();
		int index;

		for(index = 1; index <= vertices_size; index++)
		{
			file << "vn " << triangular_surface->vertices[index - 1]->n[0] << " " << triangular_surface->vertices[index - 1]->n[1] << " " << triangular_surface->vertices[index - 1]->n[2] << endl;
		}
	}

	for(int i = 0; i < arr.num_elements; i++)
	{
		TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
		int vertices_size = (int)triangular_surface->vertices.size();
		int index;

		for(index = 1; index <= vertices_size; index++)
		{
			T scale = (T)0.1;
			file << "vv " << triangular_surface->vertices[index - 1]->velocity[0]*scale << " " << triangular_surface->vertices[index - 1]->velocity[1]*scale << " " << triangular_surface->vertices[index - 1]->velocity[2]*scale << endl;
		}
	}

	for(int i = 0; i < arr.num_elements; i++)
	{
		TRIANGULAR_SURFACE* triangular_surface = arr.values[i];
		list<TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = triangular_surface->triangles.begin(); itr_triangle != triangular_surface->triangles.end(); itr_triangle++)
		{
			int index[3] = {(*itr_triangle)->vertices[0]->index, (*itr_triangle)->vertices[1]->index, (*itr_triangle)->vertices[2]->index};
			file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << endl;
		}
	}

	file.close();
}

void TRIANGULAR_SURFACE::WriteBinMeshFile(ARRAY<TRIANGULAR_SURFACE*>& surface_arr, int i_frame, const char* dir_path)
{
	//file name  .
	string file_name;
	string folder_name(dir_path);//TODO : folder name scripting
	file_name.append(folder_name.c_str());
	file_name.append("\\"+std::string(dir_path)+"_");//TODO : file name scripting
	file_name.append(boost::str(boost::format("%05d")%i_frame));
	file_name.append(".bin");

	cout << "bin:" << dir_path << endl;

	if(_access(folder_name.c_str(), 0)==-1)//folder not exist
	{
		_mkdir(folder_name.c_str());
	}
	FILE* fp = fopen(file_name.c_str(),"wb");	
	if(fp == NULL)
	{
		cout << i_frame << "RF Bin (mesh) file pointer is NULL" << endl; 
		return; 
	}
		

	//HEADER
	unsigned int ID_code = 0xDADADADA;
	fwrite(&ID_code, sizeof(unsigned int), 1, fp);

	unsigned int version = 4;
	fwrite(&version, sizeof(unsigned int), 1, fp);

	unsigned int chunk_code = 0xCCCCCCCC;
	fwrite(&chunk_code, sizeof(unsigned int), 1, fp);


	//1. PARTICLE VTX POSITIONS
	int n_tot_vtx = 0;
	for(int i_sur=0; i_sur < surface_arr.num_elements; i_sur++)//°¢ surface
	{
		n_tot_vtx += (int)surface_arr[i_sur]->vertices.size();
	}		
	fwrite(&n_tot_vtx, sizeof(int), 1, fp);


	int n_vtx;
	float x, y, z;
	for(int i_sur=0; i_sur < surface_arr.num_elements; i_sur++)//°¢ surface	
	{
		TRIANGULAR_SURFACE *surface = surface_arr[i_sur];
		n_vtx = (int)surface->vertices.size();
		
		for(int i_vtx = 0; i_vtx < n_vtx; i_vtx++)
		{
			surface->vertices[i_vtx]->index = i_vtx;
			x = (float)surface->vertices[i_vtx]->x[0];
			y = (float)surface->vertices[i_vtx]->x[1];
			z = (float)surface->vertices[i_vtx]->x[2];
			fwrite(&x, sizeof(float), 1, fp);
			fwrite(&y, sizeof(float), 1, fp);
			fwrite(&z, sizeof(float), 1, fp);
		}
	}	


	//2. PARTICLE FACES
	int n_tot_face = 0;
	for(int i_sur=0; i_sur<surface_arr.num_elements; i_sur++)
	{
		n_tot_face += (int)surface_arr[i_sur]->triangles.size();
	}		               
	fwrite(&n_tot_face, sizeof(int), 1, fp);
	

	int n_accum_vtx=0;
	int i, j, k;
	for(int i_sur=0; i_sur<surface_arr.num_elements; i_sur++)//°¢ surface
	{
		TRIANGULAR_SURFACE *surface = surface_arr[i_sur];
		
		list <TRIANGLE*>::iterator itr_tri;
		for(itr_tri = surface->triangles.begin(); itr_tri != surface->triangles.end(); itr_tri++)		
		{	
			i = (*itr_tri)->vertices[0]->index + n_accum_vtx;
			j = (*itr_tri)->vertices[1]->index + n_accum_vtx;
			k = (*itr_tri)->vertices[2]->index + n_accum_vtx;
			fwrite(&i, sizeof(int), 1, fp);
			fwrite(&j, sizeof(int), 1, fp);
			fwrite(&k, sizeof(int), 1, fp);
		}
		n_accum_vtx += (int)surface->vertices.size();
	}


	//3. FLUID TEXTURE
	/*[unsigned int] ; texture chunk code = 0xCCCCCC00 (**)
		[int] ; number of fluids
		loop for [number of vertices]
	loop for [number of fluids-1] ; version>=3 (***)
		[float] ; texture weight (***)
		endloop
		[float] ; X texture coordinate
		[float] ; Y texture coordinate
		[float] ; Z texture coordinate
		endloop
	//*/	
	
	
	//4. VTX VELOCITY
	unsigned int vel_chunk_code = 0xCCCCCC11;
	fwrite(&vel_chunk_code, sizeof(unsigned int), 1, fp);

	float vel_x, vel_y, vel_z;
	for(int i_sur = 0; i_sur < surface_arr.num_elements; i_sur++)// surface
	{
		TRIANGULAR_SURFACE *surface = surface_arr[i_sur];
		n_vtx = (int)surface->vertices.size();

		for(int i_vtx = 0; i_vtx < n_vtx; i_vtx++)
		{			
			vel_x = (float)surface->vertices[i_vtx]->velocity[0];
			vel_y = (float)surface->vertices[i_vtx]->velocity[1];
			vel_z = (float)surface->vertices[i_vtx]->velocity[2];
			fwrite(&vel_x, sizeof(float), 1, fp);
			fwrite(&vel_y, sizeof(float), 1, fp);
			fwrite(&vel_z, sizeof(float), 1, fp);
		}			
			
	}
	
	//end of file mark
	unsigned int eof_mark = 0xDEDEDEDE;
	fwrite(&eof_mark, sizeof(unsigned int), 1, fp);
	
	fflush(fp);
	fclose(fp);
}

void TRIANGULAR_SURFACE::GetGeometryInfo(vector<VT>& polygon_vertices, vector<VI>& faces)
{
	int index;
	
	for(index = 1; index <= (int)TRIANGULAR_SURFACE::vertices.size(); index++)
	{
		TRIANGULAR_SURFACE::vertices[index - 1]->index = index;
		polygon_vertices.push_back(VT(vertices[index - 1]->x[0], vertices[index - 1]->x[1], vertices[index - 1]->x[2]));
	}

	TRAVERSE_TRIANGLES
	{
		faces.push_back(VI((*itr_triangle)->vertices[0]->index, (*itr_triangle)->vertices[1]->index, (*itr_triangle)->vertices[2]->index));
	}
}











	





		



	








		


















	


	

































	

