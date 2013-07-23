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


#ifndef SOLVER_NULLSPACE_H
#define SOLVER_NULLSPACE_H

#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"

#include "determinant_derivative.hpp"
#include "programConfiguration.hpp"
#include "postProcessing.hpp"
#include "witnessSet.hpp"
#include "missing_bertini_headers.hpp"

class nullspace_config
{
	
public:
	int target_dim;   // r			the dimension of the real set we are looking for
	int ambient_dim;  // k			the dimension of the complex component we are looking IN.  
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	int num_jac_equations;    // the number of equations (after the randomization)
	
	int num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
	int num_x_vars;  // N-k+\ell
	
	int num_randomized_eqns;	// N-k (N-ambient_dim)
	int max_degree;						// the max degree of differentiated (randomized) functions
	int *randomized_degrees; // the degrees of the randomized functions (not derivatives)
	
	vec_mp **starting_linears;	// outer layer should have as many as there are randomized equations (N-k)
															// inside layer has number corresponding to max of randomized_degrees
	
	int num_additional_linears;
	vec_mp *additional_linears_terminal;
	vec_mp *additional_linears_starting;
	
	
	int num_v_linears;   // # is  (N), cause there are N equations in the subsystem. 
	vec_mp *v_linears;
	
	vec_mp v_patch; // length of this should be N-k+\ell
	
	vec_mp *target_projection; // # of these should be \ell
	
	mat_mp randomizer_matrix;  // R, the main randomizer matrix, which was passed in.  randomizes f and Jf down to N-k equations.
	
	mat_mp post_randomizer_matrix;  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations
	
	void clear();
	
	nullspace_config(){
		target_dim = ambient_dim = target_crit_codim = num_jac_equations = -1;
		num_x_vars = num_v_vars = -1;
		num_randomized_eqns = max_degree = -1;
		
		randomized_degrees = NULL;
		
		starting_linears = NULL;
		
		additional_linears_starting = additional_linears_terminal = NULL;
		
		num_v_linears = -1;
		v_linears = NULL;
		
		target_projection = NULL;
	}
	
	void print();
};



//
//void cp_nullspace_config(nullspace_config *ns_out, nullspace_config ns_in){
//	int ii;
//	ns_out->target_dim = ns_in.target_dim;
//	ns_out->ambient_dim = ns_in.ambient_dim;
//	ns_out->target_crit_codim = ns_in.target_crit_codim;
//	
//	ns_out->num_v_vars = ns_in->num_v_vars;
//	ns_out->num_x_vars = ns_in->num_x_vars;
//	
//	ns_out->num_randomized_eqns = ns_in.num_randomized_eqns;
//	ns_out->max_degree = ns_in.max_degree;
//	
//	ns_out->starting_linears = (vec_mp **) br_malloc(ns_in.num_x_vars*sizeof(vec_mp *));
//	for (ii=0; ii<ns_in.num_x_vars; ii++) {
//		ns_out->starting_linears[ii] = (vec_mp *) br_malloc(ns_in.max_degree*sizeof(vec_mp));
//		for (jj=0; jj<ns_in.max_degree; jj++) {
////			init_vec_mp();
//		}
//	}
//	
//	
//}



// the mp version
// this must be defined before the double version, because double has mp.
typedef struct
{
	
  patch_eval_data_mp patch; // patch in x
	
  preproc_data preProcData; // information related to the SLP for system
	
	prog_t *SLP; // the SLP
	
	mpq_t *gamma_rat; // randomizer
	comp_mp gamma;    // randomizer
	
	mat_mp randomizer_matrix;     // randomizer
	mat_mp randomizer_matrix_full_prec;  // randomizer
	
	int num_jac_equations;
	int target_dim;   // r			the dimension of the real set we are looking for
	int ambient_dim;  // k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	
	int num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
	int num_x_vars;  // N-k+\ell
	
	int num_randomized_eqns;	// N-k (N-ambient_dim)
	int max_degree;						// the max degree of differentiated (randomized) functions
	int *randomized_degrees;

	
	int num_additional_linears;
	vec_mp *additional_linears_terminal;
	vec_mp *additional_linears_terminal_full_prec;

	vec_mp *additional_linears_starting;
	vec_mp *additional_linears_starting_full_prec;
	
	
	mat_mp post_randomizer_matrix;  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations
	mat_mp post_randomizer_matrix_full_prec;  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations
	
	
	
	
	vec_mp **starting_linears; // outer layer should have as many as there are randomized equations
	// inside layer has number corresponding to randomized_degrees
	vec_mp **starting_linears_full_prec; // outer layer should have as many as there are randomized equations
	// inside layer has number corresponding to randomized_degrees
	
	
	int num_v_linears;
	vec_mp *v_linears;         // should be as many in here as there are randomized equations
	vec_mp *v_linears_full_prec;         // should be as many in here as there are randomized equations

	vec_mp v_patch;
	vec_mp v_patch_full_prec;
	
	
	mat_mp jac_with_proj;
	mat_mp jac_with_proj_full_prec;
	
	
	comp_mp perturbation;
	comp_mp perturbation_full_prec;
	comp_mp half;
	comp_mp half_full_prec;
////	vec_mp *source_projection;
	vec_mp *target_projection; // # of these should be target_dim (for now)
//
////	vec_mp *source_projection_full_prec;
	vec_mp *target_projection_full_prec; // # of these should be target_dim (for now)
	
	int num_variables;
	
#ifdef printpathlinprod
	FILE *FOUT;
	int num_steps;
#endif
	
	int curr_prec;
} nullspacejac_eval_data_mp;



typedef struct
{
	
	nullspacejac_eval_data_mp *BED_mp; // used only for AMP
	
	
  patch_eval_data_d patch;
	preproc_data preProcData;
	prog_t *SLP;
	
	comp_d gamma;

	mat_d randomizer_matrix;     // randomizer
	
	int num_jac_equations;
	int target_dim;   // r			the dimension of the real set we are looking for
	int ambient_dim;  // k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	
	int num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
	int num_x_vars;  // N-k+\ell
	
	int num_randomized_eqns;	// N-k (N-ambient_dim)
	int max_degree;						// the max degree of differentiated (randomized) functions
	int *randomized_degrees;
	
	
	int num_additional_linears;
	vec_d *additional_linears_terminal;
	vec_d *additional_linears_starting;
	
	mat_d post_randomizer_matrix;  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations

	
	
	
	vec_d **starting_linears; // outer layer should have as many as there are randomized equations
	// inside layer has number corresponding to randomized_degrees
	vec_d **starting_linears_full_prec; // outer layer should have as many as there are randomized equations
	// inside layer has number corresponding to randomized_degrees
	
	
	int num_v_linears;
	vec_d *v_linears;         // should be as many in here as there are randomized equations
	
	vec_d v_patch;
	
	mat_d jac_with_proj;
	vec_d *target_projection; // 
	

	int num_variables;

	comp_d perturbation;
	comp_d half;
	
#ifdef printpathlinprod
	FILE *FOUT;
	int num_steps;
#endif
	
	
} nullspacejac_eval_data_d;







/** the main function for finding critical conditions WRT a projection
 */

int nullspacejac_solver_main(int										MPType,
														 witness_set						& W, // carries with it the start points, and the linears.
														 witness_set						*W_new, // new data goes in here
														 nullspace_config				*ns_config,
														 solver_configuration		*solve_options);



int nullspacejac_solver_d(int											MPType, //, double parse_time, unsigned int currentSeed
													witness_set							& W,  // includes the initial linear.
													witness_set							*W_new,
													nullspace_config				*ns_config,
													solver_configuration		*solve_options);



void nullspacejac_track_d(trackingStats *trackCount,
													FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
													witness_set & W,
													post_process_t *endPoints,  // for holding the produced data.
													FILE *FAIL,
													int pathMod, tracker_config_t *T,
													nullspacejac_eval_data_d *ED_d,
													nullspacejac_eval_data_mp *ED_mp,
													int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
													int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													int (*change_prec)(void const *, int),
													int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
													solver_configuration *solve_options);

void nullspacejac_track_path_d(int pathNum, endgame_data_t *EG_out,
															 point_data_d *Pin,
															 FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
															 void const *ED_d, void const *ED_mp,
															 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
															 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
															 int (*change_prec)(void const *, int),
															 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));



int nullspacejac_setup_d(FILE **OUT, boost::filesystem::path outName,
												 FILE **midOUT, boost::filesystem::path midName,
												 tracker_config_t *T,
												 nullspacejac_eval_data_d *ED,
												 prog_t *dummyProg,
												 int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
												 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												 boost::filesystem::path preprocFile, boost::filesystem::path degreeFile,
												 int findStartPts, boost::filesystem::path pointsIN, boost::filesystem::path pointsOUT,
												 witness_set & W,
												 nullspace_config				*ns_config,
												 solver_configuration *solve_options);

//the new custom evaluator for this solver

int nullspacejac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED);




void printnullspacejacRelevantData(nullspacejac_eval_data_d *ED_d, nullspacejac_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP);


void nullspacejac_eval_clear_d(nullspacejac_eval_data_d *ED, int clearRegen, int MPType);



void setup_nullspacejac_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
															FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
															FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
															FILE ***NONSOLN_copy, FILE *NONSOLN,
															tracker_config_t **T_copy, tracker_config_t *T,
															nullspacejac_eval_data_d **BED_copy, nullspacejac_eval_data_d *ED_d, nullspacejac_eval_data_mp *ED_mp);


void clear_nullspacejac_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, nullspacejac_eval_data_d **BED_copy);



void setupnullspacejacEval_d(tracker_config_t *T,char preprocFile[], char degreeFile[], prog_t *dummyProg,
														 int squareSize, int patchType, int ssType, int MPType,
														 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,// what are these supposed to point to?
														 nullspacejac_eval_data_d *BED, int adjustDegrees,
														 witness_set & W,
														 nullspace_config *ns_config,
														 solver_configuration *solve_options);






void cp_nullspacejac_eval_data_d(nullspacejac_eval_data_d *BED, nullspacejac_eval_data_d *BED_d_input, nullspacejac_eval_data_mp *BED_mp_input, int MPType);


int nullspacejac_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);



int change_nullspacejac_eval_prec(void const *ED, int new_prec);

int nullspacejac_solver_mp(int MPType,  
													 witness_set & W,  // includes the initial linears.
													 witness_set *W_new,
													 nullspace_config *ns_config,
													 solver_configuration *solve_options);



void nullspacejac_track_mp(trackingStats *trackCount,
													 FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
													 witness_set & W,
													 post_process_t *endPoints,
													 FILE *FAIL,
													 int pathMod, tracker_config_t *T,
													 nullspacejac_eval_data_mp *ED,
													 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													 int (*change_prec)(void const *, int),
													 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
													 solver_configuration *solve_options);


void nullspacejac_track_path_mp(int pathNum, endgame_data_t *EG_out,
																point_data_mp *Pin,
																FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
																void const *ED,
																int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
																int (*change_prec)(void const *, int),
																int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));


int nullspacejac_setup_mp(FILE **OUT, boost::filesystem::path outName,
													FILE **midOUT, boost::filesystem::path midName,
													tracker_config_t *T,
													nullspacejac_eval_data_mp *ED,
													prog_t *dummyProg,
													int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
													int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
													int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													boost::filesystem::path preprocFile, boost::filesystem::path degreeFile,
													int findStartPts,
													boost::filesystem::path pointsIN, boost::filesystem::path pointsOUT,
													witness_set & W,
													nullspace_config *ns_config,
													solver_configuration *solve_optionss);


int nullspacejac_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);



void nullspacejac_eval_clear_mp(nullspacejac_eval_data_mp *ED, int clearRegen, int MPType);


void setup_nullspacejac_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
															 FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
															 FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
															 FILE ***NONSOLN_copy, FILE *NONSOLN,
															 tracker_config_t **T_copy, tracker_config_t *T,
															 nullspacejac_eval_data_mp **BED_copy, nullspacejac_eval_data_mp *ED_mp);


void clear_nullspacejac_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, nullspacejac_eval_data_mp **BED_copy);


void setupnullspacejacEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg,
															int squareSize, int patchType, int ssType, int prec,
															void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
															nullspacejac_eval_data_mp *BED, int adjustDegrees,
															witness_set & W,
															nullspace_config				*ns_config,
															solver_configuration *solve_options);

void cp_nullspacejac_eval_data_mp(nullspacejac_eval_data_mp *BED, nullspacejac_eval_data_mp *BED_mp_input, int MPType);



int check_issoln_nullspacejac_d(endgame_data_t *EG,
																tracker_config_t *T,
																void const *ED);
int check_issoln_nullspacejac_mp(endgame_data_t *EG,
																 tracker_config_t *T,
																 void const *ED);

int check_isstart_nullspacejac_d(point_d testpoint,
																 tracker_config_t *T,
																 void const *ED);


void check_nullspace_evaluator(point_mp current_values,
															 void const *ED);



#endif


