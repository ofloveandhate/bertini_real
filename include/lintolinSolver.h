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


#ifndef _LINTOLINSOLVER_H
#define _LINTOLINSOLVER_H

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"
#include "programStartup.h"
#include "postProcessing.h"
#include "witnessSet.h"


// the mp version
typedef struct
{
	
	
  patch_eval_data_mp patch;
  preproc_data preProcData;

	prog_t *SLP;
	
	mpq_t *gamma_rat;
	comp_mp gamma;
	mat_mp n_minusone_randomizer_matrix;
	mat_mp n_minusone_randomizer_matrix_full_prec;
	
	vec_mp current_linear;
	vec_mp current_linear_full_prec;
	vec_mp old_linear;
	vec_mp old_linear_full_prec;
	
	int num_linears;
	int num_variables;
	int curr_prec;
	
	
#ifdef printpathlintolin
	FILE *FOUT;
	int num_steps;
#endif
	
} lintolin_eval_data_mp;

//derived from basic_eval_data_d



// HERE which are the unnecessary fields???  DAB
typedef struct
{
//  square_system_eval_data_d squareSystem;
  patch_eval_data_d patch; // ???
//  start_system_eval_data_d startSystem;// ???
  preproc_data preProcData;
  lintolin_eval_data_mp *BED_mp; // used only for AMP
//  eqData_t *EqD;              // for equation-by-equation  ???
	prog_t *SLP;
	comp_d gamma;
	mat_d n_minusone_randomizer_matrix;
	
	vec_d current_linear;
	vec_d old_linear;
	int num_variables;
	int num_linears;
	
#ifdef printpathlintolin
	FILE *FOUT;
	int num_steps;
#endif
	
} lintolin_eval_data_d;





//typedef struct
//{
//  basic_eval_data_d *BED_d;
//} lintolin_eval_data_d;
////derived from basic_eval_data_d





/** the main function for finding critical conditions WRT a projection
 */
int lintolin_solver_main(int MPType,
													 witness_set W,
													 mat_mp n_minusone_randomizer_matrix_full_prec,
													 vec_mp *new_linears_full_prec,
													 int num_new_linears,
													 witness_set *W_new,
													 solver_configuration *solve_options);




int lintolin_solver_d(int MPType,
												witness_set W,  // should include the old linear
												mat_mp n_minusone_randomizer_matrix_full_prec,
												vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												int num_new_linears,
												witness_set *W_new,
												solver_configuration *solve_options);



void lintolin_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set W,
												vec_mp *new_linears_full_prec,
												int num_new_linears,
												post_process_t *endPoints,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												lintolin_eval_data_d *ED_d, lintolin_eval_data_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),solver_configuration *solve_options);

void lintolin_track_path_d(int pathNum,
														 endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
														 void const *ED_d, void const *ED_mp,
														 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 int (*change_prec)(void const *, int),
														 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));



int lintolin_setup_d(FILE **OUT, char *outName,
											 FILE **midOUT, char *midName,
											 tracker_config_t *T,
											 lintolin_eval_data_d *ED,
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

int lintolin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);



void printlintolinRelevantData(lintolin_eval_data_d *ED_d, lintolin_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP);


void lintolin_eval_clear_d(lintolin_eval_data_d *ED, int clearRegen, int MPType);



void setup_lintolin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														lintolin_eval_data_d **BED_copy,
														lintolin_eval_data_d *ED_d, lintolin_eval_data_mp *ED_mp);


void clear_lintolin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, lintolin_eval_data_d **BED_copy);



void setuplintolinEval_d(tracker_config_t *T,
												 char preprocFile[], char degreeFile[],
											 prog_t *dummyProg, int squareSize, int patchType, int ssType, int MPType,
											 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
											 lintolin_eval_data_d *BED, int adjustDegrees,
												 mat_mp n_minusone_randomizer_matrix_full_prec,
												 witness_set W,
												 solver_configuration *solve_options);


void start_system_eval_data_clear_d(start_system_eval_data_d *SSED);//actually lives in bertini library...  testing if this works.

void patch_eval_data_clear_d(patch_eval_data_d *PED);//another which lives in bertini
void patch_eval_data_clear_mp(patch_eval_data_mp *PED);//another which lives in bertini
void changePatchPrec_mp(int new_prec, patch_eval_data_mp *PED); // in bertini



void cp_lintolin_eval_data_d(lintolin_eval_data_d *BED, lintolin_eval_data_d *BED_d_input, lintolin_eval_data_mp *BED_mp_input, int MPType);

int lintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);




int change_lintolin_eval_prec(void const *ED, int prec);
void change_lintolin_eval_prec_mp(int new_prec, lintolin_eval_data_mp *BED);





int lintolin_solver_mp(int MPType,
												 witness_set W,  // includes the initial linear.
												 mat_mp n_minusone_randomizer_matrix_full_prec,  // for randomizing down to N-1 equations.
												 vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												 int num_new_linears,
												 witness_set *W_new,
												 solver_configuration *solve_options);

void lintolin_track_mp(trackingStats *trackCount,
												 FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												 witness_set W,
												 vec_mp *new_linears_full_prec,
												 int num_new_linears,
												 post_process_t *endPoints,
												 FILE *FAIL,
												 int pathMod, tracker_config_t *T,
												 lintolin_eval_data_mp *ED_d,
												 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												 int (*change_prec)(void const *, int),
												 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
												 solver_configuration *solve_options);

//												 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),

void lintolin_track_path_mp(int pathNum, endgame_data_t *EG_out,
															point_data_mp *Pin,
															FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
															void const *ED_mp,
															int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
															int (*change_prec)(void const *, int),
															int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

//int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),

int lintolin_setup_mp(FILE **OUT, char *outName,
												FILE **midOUT, char *midName,
												tracker_config_t *T,
												lintolin_eval_data_mp *ED,
												prog_t *dummyProg,
												int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
												int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												char *preprocFile, char *degreeFile,
												int findStartPts, char *pointsIN, char *pointsOUT,
												mat_mp n_minusone_randomizer_matrix_full_prec,
												witness_set W,
												solver_configuration *solve_options);

int lintolin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);

void lintolin_eval_clear_mp(lintolin_eval_data_mp *ED, int clearRegen, int MPType);

void setup_lintolin_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														 FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														 FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														 FILE ***NONSOLN_copy, FILE *NONSOLN,
														 tracker_config_t **T_copy, tracker_config_t *T,
														 lintolin_eval_data_mp **BED_copy, lintolin_eval_data_mp *ED_mp);

void clear_lintolin_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, lintolin_eval_data_mp **BED_copy);

void setuplintolinEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg,
													int squareSize, int patchType, int ssType, int prec,
													void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
													lintolin_eval_data_mp *BED, int adjustDegrees,
													mat_mp n_minusone_randomizer_matrix_full_prec,
													witness_set W,
													solver_configuration *solve_options);


//void cp_lintolin_eval_data_d(lintolin_eval_data_d *BED, lintolin_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType);


void cp_lintolin_eval_data_mp(lintolin_eval_data_mp *BED, lintolin_eval_data_mp *BED_mp_input,  int MPType);



int check_issoln_lintolin_d(endgame_data_t *EG,
														tracker_config_t *T,
														void const *ED);
int check_issoln_lintolin_mp(endgame_data_t *EG,
														 tracker_config_t *T,
														 void const *ED);



#endif


