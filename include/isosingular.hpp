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

#include <sstream>


#ifndef _ISOSINGULAR_H
#define _ISOSINGULAR_H


extern "C" {
#include "polysolve.h"
}


#include "data_type.hpp"

extern "C" {
#include "partitionParse.h"
}
#include "fileops.hpp"
#include "missing_bertini_headers.hpp"
// isosingular.c
int isosingular_deflation(int *num_deflations, int **deflation_sequence,
													boost::filesystem::path inputFile,
													boost::filesystem::path point,
													std::string bertini_command, std::string matlab_command,
													int max_deflations,
													int dim, int component_number);


void createMatlabDeflation(FILE *OUT, int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs, FILE *IN, int minorSize, int *degrees, int deflation_number);

void addItems(int *numItems, char ***itemNames, int **itemLines, FILE *IN, int lineNumber);

void parse_names(int *numItems, char ***itemNames, int **itemLines, FILE *IN, char *name, int num_declarations);

void isosingular_deflation_iteration(int *declarations,
																		 boost::filesystem::path inputOutputName,
																		 std::string matlab_command, int nullSpace, int deflation_number);

void stabilization_input_file(boost::filesystem::path outputFile,
															boost::filesystem::path funcInput,
															boost::filesystem::path configInput);


void check_declarations(int *declarations);

#endif

