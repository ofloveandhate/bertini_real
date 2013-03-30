#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <mpfr.h>
#include <mpf2mpfr.h>


#ifndef _CURVE_SELFCONJ_H
#define _CURVE_SELFCONJ_H

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"

#include "partitionParse.h"

#include "lintolinSolver.h"
#include "linprodtodetjacSolver.h"
#include "detjactodetjacSolver.h"

//the main function for computing cell decom for a curve.  only for use on a self-conjugate component.  
void computeCurveSelfConj(char * inputFile,
													witness_set_d,
													vec_d,
													curveDecomp_d*,
													int num_vars,
													int num_var_gps,
													unsigned int currentSeed);

// gets the jacobian (homogeneous) of the functions at the point current_values.  returns mat_d jacobian.  primarily for testing.
void get_jacobian(point_d current_values,
									int MPType,
									int num_var_gps,
									prog_t SLP,
									tracker_config_t T,
									mat_d jacobian);


//read the file "deg.out" and takes the sum of the numbers appearing there. used for determining the number of lintolin solves to perform to get the critical points WRT the projection and coordinate axes (bounding box).
int get_sum_degrees(char filename[], int num_funcs);



#endif
