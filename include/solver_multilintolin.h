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

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"
#include "programStartup.h"
#include "postProcessing.h"
#include "witnessSet.h"
#include "missing_bertini_headers.h"

enum {MESHSTYLE=0, LISTSTYLE=1};






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
	
	mat_mp n_minusone_randomizer_matrix;
	mat_mp n_minusone_randomizer_matrix_full_prec;
	
	
#ifdef printpathmultilintolin
	FILE *FOUT;
	int num_steps;
#endif
	
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
	mat_d n_minusone_randomizer_matrix;		// for reducing the number of equations via randomization
	


	// these variables are only used via conditional compilation
	// this property is changed in the Makefile
#ifdef printpathmultilintolin
	FILE *FOUT;
	int num_steps;
#endif
	
} multilintolin_eval_data_d;




/** the main function for finding critical conditions WRT a projection
 */
int multilintolin_solver_main(int MPType,
															witness_set W,
															mat_mp n_minusone_randomizer_matrix_full_prec,
															vec_mp *new_linears_full_prec,
															witness_set *W_new,
															solver_configuration *solve_options);




int multilin_to_lin_solver_d(int MPType,
												witness_set W,  // should include the old linear
												mat_mp n_minusone_randomizer_matrix_full_prec,
												vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												witness_set *W_new,
												solver_configuration *solve_options);



void multilin_to_lin_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set W,
												vec_mp *new_linears_full_prec,
												post_process_t *endPoints,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												multilintolin_eval_data_d *ED_d, multilintolin_eval_data_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),solver_configuration *solve_options);

void multilin_to_lin_track_path_d(int pathNum,
														 endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
														 void const *ED_d, void const *ED_mp,
														 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 int (*change_prec)(void const *, int),
														 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));



int multilin_to_lin_setup_d(FILE **OUT, char *outName,
											 FILE **midOUT, char *midName,
											 tracker_config_t *T,
											 multilintolin_eval_data_d *ED,
											 prog_t *dummyProg,  // arg7
											 int **startSub, int **endSub,
											 int **startFunc, int **endFunc,
											 int **startJvsub, int **endJvsub,
											 int **startJv, int **endJv,
											 int ***subFuncsBelow,
											 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
											 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
											 char *preprocFile, char *degreeFile,
											 int findStartPts, char *pointsIN, char *pointsOUT,
											 mat_mp n_minusone_randomizer_matrix_full_prec,
											 witness_set W,
											 solver_configuration *solve_options);

//the new custom evaluator for this solver

int multilin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);



void printmultilintolinRelevantData(multilintolin_eval_data_d *ED_d, multilintolin_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP);


void multilintolin_eval_clear_d(multilintolin_eval_data_d *ED, int clearRegen, int MPType);



void setup_multilin_to_lin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														multilintolin_eval_data_d **BED_copy,
														multilintolin_eval_data_d *ED_d, multilintolin_eval_data_mp *ED_mp);


void clear_multilintolin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, multilintolin_eval_data_d **BED_copy);



void setupmultilintolinEval_d(tracker_config_t *T,
												 char preprocFile[], char degreeFile[],
											 prog_t *dummyProg, int squareSize, int patchType, int ssType, int MPType,
											 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
											 multilintolin_eval_data_d *BED, int adjustDegrees,
												 mat_mp n_minusone_randomizer_matrix_full_prec,
												 witness_set W,
												 solver_configuration *solve_options);


void start_system_eval_data_clear_d(start_system_eval_data_d *SSED);//actually lives in bertini library...  testing if this works.

void patch_eval_data_clear_d(patch_eval_data_d *PED);//another which lives in bertini
void patch_eval_data_clear_mp(patch_eval_data_mp *PED);//another which lives in bertini
void changePatchPrec_mp(int new_prec, patch_eval_data_mp *PED); // in bertini



void cp_multilintolin_eval_data_d(multilintolin_eval_data_d *BED, multilintolin_eval_data_d *BED_d_input, multilintolin_eval_data_mp *BED_mp_input, int MPType);

int multilintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);




int change_multilintolin_eval_prec(void const *ED, int prec);
void change_multilintolin_eval_prec_mp(int new_prec, multilintolin_eval_data_mp *BED);





int multilin_to_lin_solver_mp(int MPType,
												 witness_set W,  // includes the initial linear.
												 mat_mp n_minusone_randomizer_matrix_full_prec,  // for randomizing down to N-1 equations.
												 vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												 witness_set *W_new,
												 solver_configuration *solve_options);

void multilin_to_lin_track_mp(trackingStats *trackCount,
												 FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												 witness_set W,
												 vec_mp *new_linears_full_prec,
												 post_process_t *endPoints,
												 FILE *FAIL,
												 int pathMod, tracker_config_t *T,
												 multilintolin_eval_data_mp *ED_d,
												 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												 int (*change_prec)(void const *, int),
												 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
												 solver_configuration *solve_options);

//												 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),

void multilin_to_lin_track_path_mp(int pathNum, endgame_data_t *EG_out,
															point_data_mp *Pin,
															FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
															void const *ED_mp,
															int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
															int (*change_prec)(void const *, int),
															int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

//int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),

int multilin_to_lin_setup_mp(FILE **OUT, char *outName,
												FILE **midOUT, char *midName,
												tracker_config_t *T,
												multilintolin_eval_data_mp *ED,
												prog_t *dummyProg,
												int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
												int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												char *preprocFile, char *degreeFile,
												int findStartPts, char *pointsIN, char *pointsOUT,
												mat_mp n_minusone_randomizer_matrix_full_prec,
												witness_set W,
												solver_configuration *solve_options);

int multilin_to_lin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);

void multilintolin_eval_clear_mp(multilintolin_eval_data_mp *ED, int clearRegen, int MPType);

void setup_multilin_to_lin_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														 FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														 FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														 FILE ***NONSOLN_copy, FILE *NONSOLN,
														 tracker_config_t **T_copy, tracker_config_t *T,
														 multilintolin_eval_data_mp **BED_copy, multilintolin_eval_data_mp *ED_mp);

void clear_multilintolin_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, multilintolin_eval_data_mp **BED_copy);

void setupmultilintolinEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg,
													int squareSize, int patchType, int ssType, int prec,
													void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
													multilintolin_eval_data_mp *BED, int adjustDegrees,
													mat_mp n_minusone_randomizer_matrix_full_prec,
													witness_set W,
													solver_configuration *solve_options);


//void cp_multilintolin_eval_data_d(multilintolin_eval_data_d *BED, multilintolin_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType);


void cp_multilintolin_eval_data_mp(multilintolin_eval_data_mp *BED, multilintolin_eval_data_mp *BED_mp_input,  int MPType);



int check_issoln_multilintolin_d(endgame_data_t *EG,
														tracker_config_t *T,
														void const *ED);
int check_issoln_multilintolin_mp(endgame_data_t *EG,
														 tracker_config_t *T,
														 void const *ED);



#endif


