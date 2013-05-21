
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


#ifndef NULLSPACE_H_
#define NULLSPACE_H_

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"
#include "witnessSet.h"
#include "partitionParse.h"
#include "multilintolinSolver.h"
#include "solverNullspace.h"
#include "output.h"









int compute_crit_nullspace(witness_set *W_crit_real, // the returned value
													 witness_set W,
													 mat_mp n_minusone_randomizer_matrix,
													 vec_mp *random_complex_projection,
													 vec_mp *pi,
													 int *randomized_degrees,
													 int target_dim, // this should also be the number of vectors in the *pi entry
													 program_configuration *program_options,
													 solver_configuration *solve_options);




int increment_subindices(int **function_indices,
												 int **subindices,
												 int * randomized_degrees,
												 int num_inner_indices, //these *should* be implicitly available.  switch to c++ plx.
												 int num_functions);




int increment_function_indices(int **function_indices,
												int **unused_function_indices,
												int num_inner_indices,
												int num_functions);


void nullspace_config_setup(nullspace_config *ns_config,
														int target_dim,
														int *randomized_degrees,
														mat_mp randomizer_matrix,
														witness_set W);
#endif


