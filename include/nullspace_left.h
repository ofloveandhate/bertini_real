
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
#include "solverNullspace_left.h"
#include "output.h"
#include "missing_bertini_headers.h"







/**
 the main function for computing critical sets.
 
 \param W_crit_real			the main returned structure.  
 \param W								input witness_set.
 \param randomizer_matrix	randomizes the system down to the correct number of equations.
 \param pi							the set of projections to use.
 \param randomized_degrees		the degrees of the randomized functions, before differentiation.
 \param target_dim						the dimension of the object to find.
 \param program_options				holds the configuration for the main program.  is a pointer so that it is mutable.
 \param solve_options					holds the configuration for any solvers called.  is a pointer so that it is mutable.
 */
int compute_crit_nullspace(witness_set *W_crit_real, // the returned value
													 witness_set W,
													 mat_mp randomizer_matrix,
													 vec_mp *pi,
													 int *randomized_degrees,
													 int ambient_dim,
													 int target_dim, // this should also be the number of vectors in the *pi entry
													 int target_crit_dim,
													 program_configuration *program_options,
													 solver_configuration *solve_options);



/**
 increments the odometer keeping track of which indices in the set of linears.
 
 \param function_indices					the current functions we are working on building up linears for.
 \param subindices								the current linears we will use.  this is set in this function.
 \param randomized_degrees				the degrees of the pre-differentiated functions we are building up to.
 \param	num_inner_functions				cardinal number of the subset we are considering
 \param num_functions							cardinal number of the set of functions we are indexing into.
 */
int increment_subindices(int **function_indices,
												 int **subindices,
												 int * randomized_degrees,
												 int num_inner_indices, //these *should* be implicitly available.  switch to c++ plx.
												 int num_functions);



/**
 increments the monotonically increasing array of indices which keep track of which functions we are taking linears from
 
 \param function_indices					the array we are incrementing
 \param unused_function_indices		the array of unused functions.  this is the complement of function_indices.
 \param	num_inner_functions				cardinal number of the subset we are considering
 \param num_functions							cardinal number of the set of functions we are indexing into.
 */
int increment_function_indices(int **function_indices,
												int **unused_function_indices,
												int num_inner_indices,
												int num_functions);

/**
 performs the setup for the nullspace_config which is used in the compute_crit_nullspace method, and is passed into the solverNullspace. 
 
 \param ns_config						the data structure we are setting up.
 \param target_dim					the dimension of the object we are detecting.
 \param randomized_degrees	array of integers holding the degree of each equation, *before* differentiation.
 \param randomizer_matrix		input matrix which randomizes the system down to the appropriate number of equations.
 \param W										the input witness_set
 */
void nullspace_config_setup(nullspace_config *ns_config,
														vec_mp *pi, // an array of projections, the number of which is the target dimensions
														int ambient_dim,
														int target_dim,
														int target_crit_codim,
														int max_degree,
														mat_mp randomizer_matrix,
														witness_set W,
														solver_configuration *solve_options);
#endif


