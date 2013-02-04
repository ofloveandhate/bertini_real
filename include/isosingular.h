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

#include "memory.h"
#include "data_type.h"
#include "partitionParse.h"

#ifndef _ISOSINGULAR_H
#define _ISOSINGULAR_H


// isosingular.c
int isosingular_deflation(int *num_deflations, int **deflation_sequence, char *inputFile, char *point, char *bertini_command, char *matlab_command, int max_deflations);




#endif

