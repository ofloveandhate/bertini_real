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


#ifndef _MULTILINSOLVER_H
#define _MULTILINSOLVER_H

#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"
#include "programConfiguration.hpp"
#include "postProcessing.hpp"
#include "witnessSet.hpp"
#include "missing_bertini_headers.hpp"





// the mp version
typedef struct
{

	int num_linears;
	int num_variables;
	
	
	
	vec_mp *current_linear;						// has current precision
	vec_mp *current_linear_full_prec; // carries full precision data for AMP
	vec_mp *old_linear;								// has current precision
	vec_mp *old_linear_full_prec;			// carries full precision data for AMP
	
	
	int curr_prec;
	
	
  patch_eval_data_mp patch;
  preproc_data preProcData;

	prog_t *SLP;
	
	mpq_t *gamma_rat;
	comp_mp gamma;
	
	mat_mp randomizer_matrix;
	mat_mp randomizer_matrix_full_prec;
	
	
#ifdef printpathmultilintolin
	FILE *FOUT;
	int num_steps;
#endif
	
	int verbose_level;
	
} multilintolin_eval_data_mp;



typedef struct
{
	
	int num_variables;
	int num_linears;
	
	vec_d *current_linear;
	vec_d *old_linear;
	
  patch_eval_data_d patch;

  preproc_data preProcData;
  multilintolin_eval_data_mp *BED_mp;		// used only for AMP

	prog_t *SLP;													// straight-line program
	comp_d gamma;													// \gamma randomizer
	mat_d randomizer_matrix;		// for reducing the number of equations via randomization
	


	// these variables are only used via conditional compilation
	// this property is changed in the Makefile
#ifdef printpathmultilintolin
	FILE *FOUT;
	int num_steps;
#endif
	
	int verbose_level;
	
} multilintolin_eval_data_d;




/** the main function for finding critical conditions WRT a projection
 */
int multilintolin_solver_main(int MPType,
															witness_set & W,
															mat_mp randomizer_matrix_full_prec,
															vec_mp *new_linears_full_prec,
															witness_set *W_new,
															solver_configuration *solve_options);




int multilin_to_lin_solver_d(int MPType,
												witness_set & W,  // should include the old linear
												mat_mp randomizer_matrix_full_prec,
												vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												witness_set *W_new,
												solver_configuration *solve_options);



void multilin_to_lin_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set & W,
												vec_mp *new_linears_full_prec,
												post_process_t *endPoints,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												multilintolin_eval_data_d *ED_d, multilintolin_eval_data_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),solver_configuration *solve_options);



int multilin_to_lin_setup_d(FILE **OUT, boost::filesystem::path outName,
														FILE **midOUT, boost::filesystem::path midName,
														tracker_config_t *T,
														multilintolin_eval_data_d *ED,
														prog_t *dummyProg,
														int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
														int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														boost::filesystem::path preprocFile, boost::filesystem::path degreeFile,
														int findStartPts, boost::filesystem::path pointsIN, boost::filesystem::path pointsOUT,
														mat_mp randomizer_matrix_full_prec,
														witness_set & W,
														solver_configuration *solve_options);

//the new custom evaluator for this solver

int multilin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);





void multilintolin_eval_clear_d(multilintolin_eval_data_d *ED, int clearRegen, int MPType);




void setupmultilintolinEval_d(tracker_config_t *T,
															char preprocFile[], char degreeFile[],
															prog_t *dummyProg,
															multilintolin_eval_data_d *BED,
															mat_mp randomizer_matrix_full_prec,
															witness_set & W,
															solver_configuration *solve_options);



void cp_multilintolin_eval_data_d(multilintolin_eval_data_d *BED, multilintolin_eval_data_d *BED_d_input, multilintolin_eval_data_mp *BED_mp_input, int MPType);

int multilintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);




int change_multilintolin_eval_prec(void const *ED, int new_prec);




int multilin_to_lin_solver_mp(int MPType,
												 witness_set & W,  // includes the initial linear.
												 mat_mp randomizer_matrix_full_prec,  // for randomizing down to N-1 equations.
												 vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												 witness_set *W_new,
												 solver_configuration *solve_options);

void multilin_to_lin_track_mp(trackingStats *trackCount,
												 FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												 witness_set & W,
												 vec_mp *new_linears_full_prec,
												 post_process_t *endPoints,
												 FILE *FAIL,
												 int pathMod, tracker_config_t *T,
												 multilintolin_eval_data_mp *ED_d,
												 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												 int (*change_prec)(void const *, int),
												 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
												 solver_configuration *solve_options);




int multilin_to_lin_setup_mp(FILE **OUT, boost::filesystem::path outName,
														 FILE **midOUT, boost::filesystem::path midName,
														 tracker_config_t *T,
														 multilintolin_eval_data_mp *ED,
														 prog_t *dummyProg,
														 int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
														 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 boost::filesystem::path preprocFile, boost::filesystem::path degreeFile,
														 int findStartPts,
														 boost::filesystem::path pointsIN, boost::filesystem::path pointsOUT,
												mat_mp randomizer_matrix_full_prec,
												witness_set & W,
												solver_configuration *solve_options);

int multilin_to_lin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);

void multilintolin_eval_clear_mp(multilintolin_eval_data_mp *ED, int clearRegen, int MPType);


void setupmultilintolinEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg,
													int squareSize, int patchType, int ssType, int prec,
													void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
													multilintolin_eval_data_mp *BED, int adjustDegrees,
													mat_mp randomizer_matrix_full_prec,
													witness_set & W,
													solver_configuration *solve_options);


int check_issoln_multilintolin_d(endgame_data_t *EG,
														tracker_config_t *T,
														void const *ED);
int check_issoln_multilintolin_mp(endgame_data_t *EG,
														 tracker_config_t *T,
														 void const *ED);



#endif


