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


#ifndef _LINPROD_TO_DETJAC_H
#define _LINPROD_TO_DETJAC_H

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs

#include "fileops.h"
#include "data_type.h"

//#include "partitionParse.h"

//the amount by which to perturb for numerical derivative
#define PERTURBATION_VALUE 1e-5


// HERE which are the unnecessary fields???  DAB
typedef struct
{
	//  square_system_eval_data_d squareSystem;
  patch_eval_data_d patch; // ???
	//  start_system_eval_data_d startSystem;// ???
  preproc_data preProcData;
  basic_eval_data_mp *BED_mp; // used only for AMP
	//  eqData_t *EqD;              // for equation-by-equation  ???
	prog_t *SLP;
	comp_d gamma;
	mat_d n_minusone_randomizer_matrix;
	
	vec_d *linears;
	int num_linears;
	vec_d projection;
	int num_variables;
} linprodtodetjac_eval_data_d;
//derived from basic_eval_data_d

//typedef struct
//{
//  basic_eval_data_d *BED_d;
//} linprodtodetjac_eval_data_d;
////derived from basic_eval_data_d





/** the main function for finding critical conditions WRT a projection
 */
int linprod_to_detjac_solver_d(int MPType, double parse_time, unsigned int currentSeed,
															 witness_set_d W,  // includes the initial linear.
															 mat_d n_minusone_randomizer_matrix,  // for randomizing down to N-1 equations.
															 vec_d projection,
															 witness_set_d *W_new, // for passing the data back out of this function tree
															 int my_id, int num_processes, int headnode);



void linprod_to_detjac_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set_d W,
												witness_set_d *W_new,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												linprodtodetjac_eval_data_d *ED_d, basic_eval_data_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void linprod_to_detjac_track_path_d(int pathNum,
														 endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
														 void const *ED_d, void const *ED_mp,
														 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 int (*change_prec)(void const *, int),
														 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));



int linprod_to_detjac_setup_d(FILE **OUT, char *outName,
											 FILE **midOUT, char *midName,
											 tracker_config_t *T,
											 linprodtodetjac_eval_data_d *ED,
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
											 mat_d n_minusone_randomizer_matrix,
											 witness_set_d W,
															vec_d projection);

//the new custom evaluator for this solver

int linprod_to_detjac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int take_determinant(comp_d determinant, mat_d source_matrix);
int detjac_numerical_derivative_d(mat_d Jv, point_d current_variable_values, comp_d pathVars, vec_d projection, void const *ED);


void printlinprodtodetjacRelevantData(linprodtodetjac_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP);


void linprodtodetjac_eval_clear_d(linprodtodetjac_eval_data_d *ED, int clearRegen, int MPType);



void setup_linprod_to_detjac_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														linprodtodetjac_eval_data_d **BED_copy,
														linprodtodetjac_eval_data_d *ED_d, basic_eval_data_mp *ED_mp);


void clear_linprodtodetjac_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, linprodtodetjac_eval_data_d **BED_copy);



void setuplinprodtodetjacEval_d(char preprocFile[], char degreeFile[],
												 prog_t *dummyProg, int squareSize, int patchType, int ssType, int MPType,
												 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
												 linprodtodetjac_eval_data_d *BED, int adjustDegrees,
												 mat_d n_minusone_randomizer_matrix,
																witness_set_d W,
																vec_d projection);


void start_system_eval_data_clear_d(start_system_eval_data_d *SSED);//actually lives in bertini library...  testing if this works.

void patch_eval_data_clear_d(patch_eval_data_d *PED);//another which lives in bertini

void cp_linprodtodetjac_eval_data_d(linprodtodetjac_eval_data_d *BED, linprodtodetjac_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType);

int linprodtodetjac_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);


#endif


