#include "memory.h"
#include "data_type.h"
#include "partitionParse.h"


#ifndef _CURVE_NOTSELFCONJ_H
#define _CURVE_NOTSELFCONJ_H

double	get_comp_rand_d(comp_d);
void   	get_random_mat_d(mat_d, int,int);
void 	diag_homotopy_input_file(char*, char*,char*, char*,vec_d, int);
void 	diag_homotopy_start_file(char*, witness_point_set_d, int);
void 	computeCurveNotSelfConj(witness_set_d, vec_d, curveDecomp_d*,int, char*);

#endif
