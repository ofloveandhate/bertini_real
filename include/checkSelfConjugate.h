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


void membership_test_input_file(char *outputFile,
                                char *funcInput,
                                char *configInput,
                                int  tracktype);

int write_member_points(point_d, char * fmt);


void read_incidence_matrix(int component_numbers[]);


int checkSelfConjugate(witness_set_d W,
                       int           num_vars,
                       char          *input_file);


#endif
