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


#include "partitionParse.h"
#include "data_type.h"
#include "polysolve.h"
#include "fileops.h"


//returns 1 if self-conjugate, 0 else.
int checkSelfConjugate(witness_set_d W,
                       int           num_vars,
                       char          *input_file);

// returns the component number according to the incidence matrix
int get_component_number(witness_set_d W,
												 int           num_vars,
												 char          *input_file);

//write a single point to "member_points"
int write_member_points_singlept(point_d point_to_write, char * fmt);

//write a single point, and its complex conjugate, to "member_points"
int write_member_points_sc(point_d, char * fmt);

//write the input file to feed bertini to perform membership testing
void membership_test_input_file(char *outputFile,
                                char *funcInput,
                                char *configInput,
                                int  tracktype);

//read the incicence matrix
void read_incidence_matrix(int *component_numbers);
void read_incidence_matrix_wrt_number(int *component_numbers, int given_incidence_number);




#endif
