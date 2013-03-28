
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#ifndef _WITNESS_SET_PARSE
#define _WITNESS_SET_PARSE

#include "data_type.h"
#include "polysolve.h"
#include "fileops.h"

int witnessSetParse(witness_set_d *W, char *witness_set_file, int num_vars);
//int witnessSetParse_mp(witness_set_mp *W, char *witness_set_file, int num_vars);

#endif
