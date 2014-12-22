#include "stdafx.h"
#include "POISSON_SOLVER_3D.h"

int  POISSON_SOLVER_3D::smoother_id = 0;
bool POISSON_SOLVER_3D::use_detailed_log = false;
T*	 POISSON_SOLVER_3D::euclidean_res_accums = 0;
T*   POISSON_SOLVER_3D::infinity_res_accums = 0;
int* POISSON_SOLVER_3D::num_frames_ids = 0;