



#ifndef _CURVE_NOTSELFCONJ_H
#define _CURVE_NOTSELFCONJ_H


#include "programStartup.h"
#include "data_type.h"
#include "polysolve.h"
#include "partitionParse.h"
#include "fileops.h"


void  get_random_mat_d(mat_d, int,int);
void 	diag_homotopy_input_file(char*, char*,char*, char*,vec_d, int);
void 	diag_homotopy_start_file(char*, witness_set W, int);
void 	computeCurveNotSelfConj(witness_set W,
															vec_mp pi,
															curveDecomp_d *C,
															vertex_set			*V,
															int num_vars,
															char *input_file,
															program_configuration *program_options,
															solver_configuration * solve_options);

#endif
