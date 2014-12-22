#pragma once

// Headers those use a lot
#include <iostream>
#include <fstream>
#include <assert.h>
#include <limits>

#include "MACROS.h"
#include "VECTOR_2D.h"
#include "VECTOR_3D.h"

using namespace std;

// We want to avoid tons of overhead when we do the generic programming
#ifdef USE_FLOAT_T
typedef float T;
#else
typedef double T;
#endif

// Determine whether 2D or 3D
//#ifdef 2D_SIMULATION
//typedef VECTOR_2D<T> VT;
//typedef VECTOR_2D<int> VI;
//#else
typedef VECTOR_3D<T> VT;
typedef VECTOR_3D<int> VI;
//#endif



