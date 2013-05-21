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
#include "witnessSet.h"
#include "partitionParse.h"

#include "checkSelfConjugate.h"
#include "lintolinSolver.h"
#include "linprodtodetjacSolver.h"
#include "detjactodetjacSolver.h"

#include "nullspace.h"

#include "output.h"
/**
//the main function for computing cell decom for a curve.  only for use on a self-conjugate component.  
 
 \param inputFile the name of the input file.
 \param W	witness_set containing linears, patches, points, etc.  much info and calculation performed from this little guy.
 \param pi the set of projections to use.  in this case, there should be only 1.
 \param C		curve decomposition structure into which to place computed data.
 \param V		vertex set structure into which to place collected data.
 \param num_vars		the total number of variables for the problem.
 \param num_var_gps	the number of variable groups.  should be 1.
 \param options program configuration.
 \param solve_options solver configuration.
 */
void computeCurveSelfConj(char * inputFile,
													witness_set W,
													vec_mp *pi,
													curveDecomp_d *C,
													vertex_set *V,
													int num_vars,
													int num_var_gps,
													program_configuration *options,
													solver_configuration *solve_options);



/**
 the linprodtodetjac method for getting the critical points
 */
int compute_crit_linprodtodetjac(witness_set *W_crit_real, // the returned value
																	witness_set W,
																	mat_mp n_minusone_randomizer_matrix,
																	vec_mp random_complex_projection,
																	vec_mp pi,
																	int num_new_linears,
																	program_configuration *program_options,
																	solver_configuration *solve_options);




/**
 
 */
void check_patch_values(witness_set W);

/**
// checks to see what the determinant of the jacobian is at the points in W.
 */
void check_detjac(witness_set W, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection);

/**
// gets the jacobian (homogeneous) of the functions at the point current_values.  returns mat_d jacobian.  primarily for testing.
 */
void get_jacobian(point_d current_values,
									int MPType,
									int num_var_gps,
									prog_t SLP,
									tracker_config_t T,
									mat_d jacobian);

/**
//read the file "deg.out" and takes the sum of the numbers appearing there. used for determining the number of lintolin solves to perform to get the critical points WRT the projection and coordinate axes (bounding box).
 */
int get_sum_degrees(char filename[], int num_funcs);


void sort_for_membership(char * input_file,
												 witness_set *W_out,
												 witness_set W_in,
												 char *stifle_text);


void sort_for_unique(witness_set *W_out,
										 witness_set W_in,
										 tracker_config_t T);


void sort_for_real(witness_set *W_out,
									 witness_set W_in,
									 tracker_config_t T);

void sort_increasing_by_real(vec_mp *projections_sorted, int **index_tracker, vec_mp projections_input);

void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, int ** randomized_degrees,
																								int num_variables, int num_funcs);
int compare_integers_decreasing(const void * left_in, const void * right_in);



int verify_projection_ok(witness_set W,
												 mat_mp n_minusone_randomizer_matrix,
												 vec_mp projection,
												 solver_configuration *solve_options);


#endif
