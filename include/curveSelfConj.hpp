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

#include "missing_bertini_headers.hpp"



#include "fileops.hpp"
#include "data_type.hpp"
#include "witnessSet.hpp"
extern "C" {
#include "partitionParse.h"
}

#include "checkSelfConjugate.hpp"
#include "solver_lintolin.hpp"
#include "solver_linprodtodetjac.hpp"
#include "solver_detjactodetjac.hpp"

#include "nullspace_left.hpp"

#include "output.hpp"
/**
//the main function for computing cell decom for a curve.  only for use on a self-conjugate component.  
 
 \param inputFile the name of the input file.
 \param W	witness_set containing linears, patches, points, etc.  much info and calculation performed from this little guy.
 \param pi the set of projections to use.  in this case, there should be only 1.
 \param C		curve decomposition structure into which to place computed data.
 \param V		vertex set structure into which to place collected data.
 \param num_vars		the total number of variables for the problem.
 \param options program configuration.
 \param solve_options solver configuration.
 */
void computeCurveSelfConj(boost::filesystem::path inputFile,
													witness_set & W,
													vec_mp *pi,
													curve_decomposition &C,
													vertex_set &V,
													int num_vars,
													BR_configuration & options,
													solver_configuration & solve_options);



int slice_and_dice(witness_set & W_curve,
							witness_set & W_crit_real,
							mat_mp randomizer_matrix,
							vec_mp *pi,
							BR_configuration & program_options,
							solver_configuration solve_options,
							curve_decomposition & C,
							vertex_set & V);


int curve_compute_critical_points(witness_set & W_curve,
																	mat_mp randomizer_matrix,
																	int *randomized_degrees,
																	vec_mp *pi,
																	BR_configuration & program_options,
																	solver_configuration solve_options,
																	witness_set & W_crit_real);
/**
 the linprodtodetjac method for getting the critical points
 */
int compute_crit_linprodtodetjac(witness_set *W_crit_real, // the returned value
																	witness_set & W,
																	mat_mp n_minusone_randomizer_matrix, 
																	vec_mp pi,
																	int num_new_linears,
																	BR_configuration & program_options,
																	solver_configuration & solve_options);




int curve_get_additional_critpts(witness_set *W_crit_real,
							 witness_set & W,
							 mat_mp randomizer_matrix,
							 vec_mp pi,
							 int *randomized_degrees,
							 BR_configuration & program_options,
							 solver_configuration & solve_options);

//
///**
// 
// */
//void check_patch_values(witness_set W);

///**
// checks to see what the determinant of the jacobian is at the points in W.
// */
//void check_detjac(witness_set W, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection);
//
///**
//// gets the jacobian (homogeneous) of the functions at the point current_values.  returns mat_d jacobian.  primarily for testing.
// */
//void get_jacobian(point_d current_values,
//									int MPType,
//									int num_var_gps,
//									prog_t SLP,
//									tracker_config_t T,
//									mat_d jacobian);

/**
 read the file "deg.out" and takes the sum of the numbers appearing there. used for determining the number of lintolin solves to perform to get the critical points WRT the projection and coordinate axes (bounding box).
 */
int get_sum_degrees(char filename[], int num_funcs);


void sort_for_membership(char * input_file,
												 witness_set *W_out,
												 witness_set & W_in,
												 char *stifle_text);









int verify_projection_ok(witness_set & W,
												 vec_mp projection,
												 solver_configuration & solve_options);

int verify_projection_ok(witness_set & W,
												 mat_mp randomizer_matrix,
												 vec_mp projection,
												 solver_configuration & solve_options);

#endif
