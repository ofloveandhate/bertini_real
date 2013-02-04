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

#include "partitionParse.h"
#include "memory.h"
#include "data_type.h"


#ifndef _CHECK_SELFCONJUGATE_H
#define _CHECK_SELFCONJUGATE_H

void membership_test_input_file(char *outputFile,
                                char *funcInput,
                                char *configInput,
                                int  tracktype);

int write_member_points(point_d coord, int num_vars, char * fmt);

int write_member_points_conjugated(point_d point_to_write, int num_vars, char * fmt);

int read_incidence_matrix();


int checkSelfConjugate(witness_set_d W,
                       int           num_vars,
                       char          *input_file);


#endif
