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


#ifndef _BR_OUTPUT_H
#define _BR_OUTPUT_H

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs


#include "witnessSet.h"

#include "fileops.h"
#include "data_type.h"

#include "partitionParse.h"

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
void print_curve(curveDecomp_d C, int num_vars, char *input_deflated_Name, char *outputfile, int MPType);

/**
 
 */
void print_matrix_to_file_mp(FILE *OUT, int digits, mat_mp );
#endif


