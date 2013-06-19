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
void Output_Main(program_configuration program_options, witness_set W, curveDecomp_d C, vertex_set V);

/**
 
 */
void print_vertices(vertex_set V, int num_vars,
										char *outputfile, int MPType);

/**
 
 */
void print_edges(edge *E, int num_edges, int num_vars, char *outputfile, int MPType);


/**
 
 */
void print_curve(curveDecomp_d C, int num_vars, boost::filesystem::path input_deflated_Name, boost::filesystem::path  outputfile, int MPType);

/**
 
 */
void print_matrix_to_file_mp(FILE *OUT, int digits, mat_mp );
#endif


