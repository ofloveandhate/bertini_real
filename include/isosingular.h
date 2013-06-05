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



#ifndef _ISOSINGULAR_H
#define _ISOSINGULAR_H

#include "polysolve.h"
#include "data_type.h"
#include "partitionParse.h"
#include "fileops.h"
#include "missing_bertini_headers.h"
// isosingular.c
int isosingular_deflation(int *num_deflations, int **deflation_sequence, char *inputFile, char *point, char *bertini_command, char *matlab_command, int max_deflations);


void createMatlabDeflation(FILE *OUT, int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs, FILE *IN, int minorSize, int *degrees, int deflation_number);
void addItems(int *numItems, char ***itemNames, int **itemLines, FILE *IN, int lineNumber);
void parse_names(int *numItems, char ***itemNames, int **itemLines, FILE *IN, char *name, int num_declarations);
void isosingular_deflation_iteration(int *declarations, char *inputOutputName, char *matlab_command, int nullSpace, int deflation_number);
void stabilization_input_file(char *outputFile, char *funcInput, char *configInput);


void check_declarations(int *declarations);

#endif

