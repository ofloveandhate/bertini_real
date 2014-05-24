#ifndef _POST_PROCESSING_H
#define _POST_PROCESSING_H

/** \file postProcessing.hpp */

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





#include "bertini_headers.hpp"


#include "fileops.hpp"
#include "data_type.hpp"
#include "solver.hpp"

#include "programConfiguration.hpp"




class solver_configuration; // forward declaration

/**
 \brief Conversion function for turning endgame_data_t into post_process_t.
 
 \todo Eliminate the nonsensical conversion of data, as all the necessary data was already present in the endgame_data_t, and for (it appears) legacy reasons, it is converted in to a more useless form.  Bleh.
 
 \param endPoint Result of this function, for post_processing.
 \param EG The source for this function, converting into a post_process_t.
 */
void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG);

/**
 \brief Vertini_real's own method of finding SINGULAR solutions among a collection of post_process_t's.
 
 this singular detection does NOT declare a point singular based on multiplicity -- just on the basis of condition number.
 
 \return number of singular solutions
 \param endPoints pointers to post_process_t's containing the output from Bertini solver loops.
 \param num_sols How many solutions to test; cannot exceed the number of endPoints.
 \param num_vars The number of variables in the solutions.
 \param T Configuration of the Bertini tracker loop.
 */
int BRfindSingularSolns(post_process_t *endPoints,
						int num_sols, int num_vars,
						tracker_config_t *T );

/**
 \brief Bertini_real's own method of finding FINITE solutions among a collection of post_process_t's.
 
 \return number of finite solutions
 \param endPoints pointers to post_process_t's containing the output from Bertini solver loops.
 \param num_sols How many solutions to test; cannot exceed the number of endPoints.
 \param num_vars The number of variables in the solutions.
 \param T Configuration of the Bertini tracker loop.
 */
int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
					  tracker_config_t *T );



/**
 \brief Bertini_real's method for determining which solutions are real.
 
 \return number of real solutions
 \param endPoints pointers to post_process_t's containing the output from Bertini solver loops.
 \param num_sols How many solutions to test; cannot exceed the number of endPoints.
 \param num_vars The number of variables in the solutions.
 \param T Configuration of the Bertini tracker loop.
 */
int BRfindRealSolns(post_process_t *endPoints, int num_sols, int num_vars,
					tracker_config_t *T );


/**
 \brief Extract the solution from a post_process_t
 
 \param veccie The computed value, a pre-initted vec_mp.
 \param endPoint Result of a tracker loop.
 */
void endpoint_to_vec_mp(vec_mp veccie, post_process_t *endPoint);




#endif



