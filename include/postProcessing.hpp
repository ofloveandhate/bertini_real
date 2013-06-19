
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

#include <getopt.h>



#ifndef _POST_PROCESSING_H
#define _POST_PROCESSING_H


extern "C" {
#include "polysolve.h"
}

#include "witnessSet.hpp"
#include "fileops.hpp"
#include "data_type.hpp"
#include "solver.hpp"

#include "programConfiguration.hpp"

#include "missing_bertini_headers.hpp"



/**
 conversion function for turning endgame_data_t into post_process_t.
 */
void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG);

/**
 bertini_real's own method of finding SINGULAR solutions among a collection of post_process_t's.
 this singular detection does NOT declare a point singular based on multiplicity -- just on the basis of condition number.
 */
int BRfindSingularSolns(post_process_t *endPoints,
												int num_sols, int num_vars,
												tracker_config_t *T );

/**
 bertini_real's own method of finding FINITE solutions among a collection of post_process_t's.
 */
int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
											tracker_config_t *T );

/**
 determines if a solution is allowable.  based on solver_configuration.
 */
int is_acceptable_solution(post_process_t endPoint, solver_configuration *solve_options);


/**
 bertini_real's version of post-processing.  options are set via the solver_configuration.
 */
void BRpostProcessing(post_process_t *endPoints, witness_set *W_new, int num_pts,
											preproc_data *preProcData, tracker_config_t *T,
											solver_configuration *solve_options);




#endif



