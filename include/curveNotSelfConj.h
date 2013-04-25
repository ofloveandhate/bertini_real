



#ifndef _CURVE_NOTSELFCONJ_H
#define _CURVE_NOTSELFCONJ_H



#include "data_type.h"
#include "polysolve.h"
#include "partitionParse.h"
#include "fileops.h"


void  get_random_mat_d(mat_d, int,int);
void 	diag_homotopy_input_file(char*, char*,char*, char*,vec_d, int);
void 	diag_homotopy_start_file(char*, witness_set W, int);
void 	computeCurveNotSelfConj(witness_set, vec_mp, curveDecomp_d*,int, char*);

#endif
