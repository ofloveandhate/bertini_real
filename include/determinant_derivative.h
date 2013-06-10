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


#ifndef _DETERMINANT_DERIVATIVE_H
#define _DETERMINANT_DERIVATIVE_H

#include "solver_linprodtodetjac.h"
#include "solver_detjactodetjac.h"
#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"

#include "missing_bertini_headers.h"


//the amount by which to perturb for numerical derivative
#define PERTURBATION_VALUE 1e-5
#define PERTURBATION_VALUE_mp 1e-8

#define TOL_DOUBLE_PRECISION 1e-13
#define LARGECHANGE_DOUBLEPRECISION 1e14


#define TOL_MP 1e-40
#define LARGECHANGE_MP 1e50

int take_determinant_d(comp_d determinant, mat_d source_matrix);
int take_determinant_mp(comp_mp determinant, mat_mp source_matrix);

//patch_eval_data_d patch,

int detjac_numerical_derivative_d(mat_d Jv, point_d current_variable_values, comp_d pathVars, vec_d projection,
																	int num_variables,
																	prog_t *SLP,
																	mat_d n_minusone_randomizer_matrix);

//																	 patch_eval_data_mp patch,

int detjac_numerical_derivative_mp(mat_mp Jv, //  the returned value
																	 point_mp current_variable_values, comp_mp pathVars, vec_mp projection,
																	 int num_variables,
																	 prog_t *SLP,
																	 mat_mp n_minusone_randomizer_matrix); // inputs





#endif



