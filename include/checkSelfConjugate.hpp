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



#ifndef _CHECK_SELFCONJUGATE_H
#define _CHECK_SELFCONJUGATE_H


extern "C" {
#include "polysolve.h"
}


#include "witnessSet.hpp"
//extern "C" {
//#include "partitionParse.h"
//}
#include "data_type.hpp"
#include "fileops.hpp"
#include "missing_bertini_headers.hpp"

/**
 main function for performing self-conjugacy check.
 
 returns 1 if self-conjugate, 0 else.
 
 \param W input witness set
 \param num_vars the number of variables in the problem, including the homogeneous one.
 \param input_file the name of the bertini input file
 \param stifle_text text to append to the system call, possibly " < /dev/null"
 
 \return boolean integer indicating whether the component is self-conjugate.
 
 */
int checkSelfConjugate(witness_set & W,
                       int           num_vars,
                      boost::filesystem::path input_file,
											 std::string stifle_text);




/**
 returns the incidence number according to the incidence matrix
 
 \param W the input witness set
 \param num_vars the number of variables in the problem, including homogeneous.
 \param input_file the name of the bertini input file
 \param stifle_text text to append to the command to stifle screen output.
 */
int get_incidence_number(const witness_set & W,
												 int           num_vars,
												 boost::filesystem::path input_file,
												 std::string stifle_text);

/**
 write a single point to "member_points"
 
 \param point_to_write the homogeneous point to write
 \param the format
 */
int write_member_points_singlept(point_mp point_to_write);

/**
  write a single point, and its complex conjugate, to "member_points"
 
 \param point_to_write the point to write
 \param fmt the format of the string.
 */
int write_member_points_sc(point_mp point_to_write);

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
void read_incidence_matrix(int *component_numbers);

/**
 //read the incicence matrix
 
 returns a logical integer array, indicating whether member_point is on the given_incidence_number.  That is, the returned array (which must be initialized before passing into this function), will be 0 at all incidence number positions for which the point is not a member, and at least 1 for all points it is a member.
 */
void read_incidence_matrix_wrt_number(int *component_numbers, int given_incidence_number);




#endif
