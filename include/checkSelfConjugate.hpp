#ifndef _CHECK_SELFCONJUGATE_H
#define _CHECK_SELFCONJUGATE_H


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





#include "programConfiguration.hpp"
#include "fileops.hpp"
#include "bertini_headers.hpp"
#include "bertini_extensions.hpp"




/**
 main function for performing self-conjugacy check.
 
 returns 1 if self-conjugate, 0 else.
 
 \param W input witness set
 \param num_vars the number of variables in the problem, including the homogeneous one.
 \param input_file the name of the bertini input file
 
 \return boolean integer indicating whether the component is self-conjugate.
 
 */
bool checkSelfConjugate(vec_mp test_point,
						BR_configuration & program_options,
						boost::filesystem::path input_file);




/**
 returns the incidence number according to the incidence matrix
 
 \param W the input witness set
 \param num_vars the number of variables in the problem, including homogeneous.
 \param input_file the name of the bertini input file
 */
int get_incidence_number(vec_mp test_point, BR_configuration & program_options, boost::filesystem::path input_file);

/**
 write a single point to "member_points"
 
 \param point_to_write the homogeneous point to write
 */
int write_member_points_singlept(vec_mp point_to_write);

/**
 write a single point, and its complex conjugate, to "member_points"
 
 \param point_to_write the point to write
 */
int write_member_points_sc(vec_mp point_to_write);

/**
 write the input file to feed bertini to perform membership testing
 
 \param outputFile the name of the file to write
 \param funcInput the name of the func_input file
 \param configInput the name of the config file
 \param tracktype bertini track type.  should be 3.
 */
void membership_test_input_file(boost::filesystem::path outputFile,
                                boost::filesystem::path funcInput,
                                boost::filesystem::path configInput,
                                int  tracktype);

/**
 read the incicence_matrix file.  return the incidence number for the member_points
 
 \param component_numbers The returned value of this function.
 */
std::vector<int> read_incidence_matrix();

/**
 read the incicence matrix
 
 returns a logical integer array, indicating whether member_point is on the given_incidence_number.  That is, the returned array (which must be initialized before passing into this function), will be 0 at all incidence number positions for which the point is not a member, and at least 1 for all points it is a member.
 */
void read_incidence_matrix_wrt_number(int *component_numbers, int given_incidence_number);




#endif
