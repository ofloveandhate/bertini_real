#ifndef _POST_PROCESSING_H
#define _POST_PROCESSING_H

/** \file postProcessing.hpp 

\brief For adding metadata to output from the various solvers.
*/

#include <boost/bind.hpp>

#include "bertini1/bertini_headers.hpp"

#include "io/fileops.hpp"
#include "nag/solvers/solver.hpp"
#include "programConfiguration.hpp"




class SolverConfiguration; // forward declaration

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
\brief Read a midpath-style file, get which paths appear to be crossed at the endgame boundary

\param num_paths The number of paths to expect.  Not in header for the file.
\param num_vars The number of variables per point.  Not in header for file.
\param tol The tolerance for infinity-norm of difference in points to be labeled as crossing.
\param [out] num_crossings The number of apparent path crossings
\param filename The name of the file to open.
\param crossed_indices The indices of the crossed paths.

Derived directly from the bertini1 function midpoint_checker
*/
void BRmidpointChecker(int num_paths, int num_vars, double tol, int *num_crossings, boost::filesystem::path const& filename, std::vector<int> & crossed_indices);

/**
 \brief Extract the solution from a post_process_t
 
 \param veccie The computed value, a pre-initted vec_mp.
 \param endPoint Result of a tracker loop.
 */
void endpoint_to_vec_mp(vec_mp veccie, post_process_t *endPoint);




#endif



