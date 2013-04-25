
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#ifndef _WITNESS_SET_PARSE
#define _WITNESS_SET_PARSE

#include "data_type.h"
#include "polysolve.h"
#include "fileops.h"

void get_variable_names(witness_set *W);

int witnessSetParse(witness_set *W, char *witness_set_file, int num_vars);

#endif
