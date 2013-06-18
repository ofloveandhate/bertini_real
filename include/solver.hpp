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


#ifndef SOLVER_MAIN_HEADER_H
#define SOLVER_MAIN_HEADER_H

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"

#include "determinant_derivative.h"
#include "programConfiguration.h"
#include "postProcessing.h"
#include "witnessSet.h"
#include "missing_bertini_headers.h"


typedef struct
{
	void *asdf;
	
	int (*evaluator_function_d) (point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *); // function handle to evaluator to use
  int (*evaluator_function_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
	
	int (*precision_changer)(void const *ED, int new_prec);
	
} solver_handle_struct;




void generic_track_path_d(int pathNum, endgame_data_t *EG_out,
													point_data_d *Pin,
													FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
													void const *ED_d, void const *ED_mp,
													int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
													int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													int (*change_prec)(void const *, int),
													int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));



#endif