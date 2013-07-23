#include <dirent.h>

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

#include <boost/filesystem/path.hpp>

#ifndef _BR_OUTPUT_H
#define _BR_OUTPUT_H

extern "C" {
#include "cascade.h"
}
extern "C" {
#include "polysolve.h"
}

#include "witnessSet.hpp"

#include "fileops.hpp"
#include "data_type.hpp"

#include "partitionParse.h"
#include "missing_bertini_headers.hpp"
#include "programConfiguration.hpp"
/**
 
 */
void Output_Main(BR_configuration program_options, witness_set & W, curve_decomposition & C, vertex_set & V);






/**
 
 */
void print_matrix_to_file_mp(FILE *OUT, int digits, mat_mp );
#endif


