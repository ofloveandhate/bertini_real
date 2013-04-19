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

#include "fileops.h"
#include "data_type.h"

#include "partitionParse.h"


void copyfile(char *INfile,char *OUTfile);

void print_vertices(vertex_d *V, int num_V, int num_vars,char *outputfile, int MPType);
void print_edges(edge_d *E, int num_E, int num_vars,char *input_deflated_Name,char *outputfile, int MPType);
void print_each_edge(edge_d E,int num_vars,FILE *OUT, int MPType);
void Output_Main(char *inputName, char *input_deflated_Name,int *deflation_sequence,int num_vars, curveDecomp_d C, int MPType);

#endif


