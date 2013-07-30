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

#include <string>
#ifndef SOLVER_MAIN_HEADER_H
#define SOLVER_MAIN_HEADER_H

extern "C" {
#include "cascade.h"
}
extern "C" {
#include "polysolve.h"
}
#include "fileops.hpp"
#include "data_type.hpp"

//#include "programConfiguration.hpp"
//#include "postProcessing.hpp"
#include "witnessSet.hpp"
#include "missing_bertini_headers.hpp"



///////////
//
//    SOLVER CONFIGURATION
//
//////////


typedef struct
{
	
	tracker_config_t T;
	preproc_data PPD;
	
	int allow_multiplicity;
	int allow_singular;
	int allow_infinite;
	int allow_unsuccess;
	
	int verbose_level;
	int show_status_summary;
	
	int use_midpoint_checker;
	double midpoint_tol;
	
	int use_gamma_trick;
	
	int complete_witness_set;
} solver_configuration;



void generic_set_start_pts(point_data_d ** startPts,
													 witness_set & W);

void generic_set_start_pts(point_data_mp ** startPts,
													 witness_set & W);

void generic_setup_patch(patch_eval_data_d *P, const witness_set & W); // for mp type 0
void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W);// for my type 2
void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W, int prec); // for mp type 1

/** reads the tracker_config_t from file. */
void get_tracker_config(solver_configuration *solve_options,int MPType);
void solver_init_config(solver_configuration *options);
void solver_clear_config(solver_configuration *options);


void generic_track_path_d(int pathNum, endgame_data_t *EG_out,
									point_data_d *Pin,
									FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
									void const *ED_d, void const *ED_mp,
									int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
									int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
									int (*change_prec)(void const *, int),
									int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void generic_track_path_mp(int pathNum, endgame_data_t *EG_out,
													 point_data_mp *Pin,
													 FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
													 void const *ED,
													 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													 int (*change_prec)(void const *, int),
													 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

class solver
{
public:
	
	// these virtual functions will need to be programmed into the derived classes.
	
	//	virtual evaluator_d();
	//	virtual evaluator_mp();
	//	virtual change_precision();
	//
	
	int (*evaluator_function_d) (point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *); // function handle to evaluator to use
	int (*evaluator_function_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
	
	int (*precision_changer)(void const *ED, int new_prec);
	
};










#endif