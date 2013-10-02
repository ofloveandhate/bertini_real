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

#include "missing_bertini_headers.hpp"



#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"
#include "programConfiguration.hpp"
#include "postProcessing.hpp"







class multilin_config
{
	
public:
	mat_mp randomizer_matrix;  ///< R, the main randomizer matrix, which was passed in.  randomizes f and Jf down to N-k equations.

	SLP_global_pointers SLP_memory;
	prog_t * SLP;
	
	int MPType;
	
	bool have_rand;
	bool have_mem;
	
	
	multilin_config(solver_configuration & solve_options,
									const witness_set & W)
	{
		init();
		
		set_memory(solve_options);
		
		make_randomizer(solve_options, W);
	}
	
	
	
	multilin_config(solver_configuration & solve_options,
									mat_mp _random)
	{
		init();
		
		set_memory(solve_options);
		
		set_randomizer(_random);
	}
	
	
	multilin_config(solver_configuration & solve_options)
	{
		init();
		set_memory(solve_options);
	}
	
	
	multilin_config(mat_mp _random)
	{
		init();
		set_randomizer(_random);
	}
	
	void make_randomizer(const solver_configuration & solve_options, const witness_set & W)
	{
		std::vector< int > randomized_degrees;
		//get the matrix and the degrees of the resulting randomized functions.
		
		
		std::cout << "making randomizer of size	" <<  W.num_variables-W.num_linears-W.num_patches << "x" << solve_options.PPD.num_funcs << std::endl;
		make_randomization_matrix_based_on_degrees(this->randomizer_matrix, randomized_degrees, W.num_variables-W.num_linears-W.num_patches, solve_options.PPD.num_funcs);
		
		have_rand = true;
	}

	
	void set_memory(solver_configuration & solve_options)
	{
		
//TODO: should i assume here that the input file is already parsed??
		this->MPType = solve_options.T.MPType;
		solve_options.T.numVars = setupProg(SLP, solve_options.T.Precision, solve_options.T.MPType);
		//make randomizer matrix here
		SLP_memory.capture_globals();
		SLP_memory.set_globals_null();
		have_mem = true;
		
		
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
	}
	
	
	
	void set_randomizer(mat_mp _random)
	{
		mat_cp_mp(randomizer_matrix, _random);
		have_rand = true;
	}
	
	
	
	
	void copy(const multilin_config & other)
	{
		
		//ah yes, this problem.
		//		this->SLP = other.SLP;// this needs to be a deep copy
	}
	
	void init()
	{
		init_mat_mp2(randomizer_matrix,1,1,1024); randomizer_matrix->rows = randomizer_matrix->cols = 1;

		SLP = new prog_t;
		have_mem = false;
		have_rand = false;
	}
	
	
	void clear()
	{
		
		
		clear_mat_mp(randomizer_matrix);
		
		SLP_memory.set_globals_to_this();
		clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
		delete SLP;
	}
	
	
	
	multilin_config()
	{
		init();
	}
	
	multilin_config(const multilin_config & other)
	{
		init();
		copy(other);
	}
	
	multilin_config & operator=(const multilin_config & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	
	~multilin_config(){
		clear();
	}
};



















// the mp version
// this must be defined before the double version, because double has mp.
class multilintolin_eval_data_mp : public solver_mp
{
public:
	
	
	int num_linears;
	
	
	
	vec_mp *current_linear;						// has current precision
	vec_mp *current_linear_full_prec; // carries full precision data for AMP
	vec_mp *old_linear;								// has current precision
	vec_mp *old_linear_full_prec;			// carries full precision data for AMP
	
	
	
	// default initializer
	multilintolin_eval_data_mp() : solver_mp(){
		init();
	}
	
	multilintolin_eval_data_mp(int mp) : solver_mp(mp){
		this->MPType = mp;
		init();
	}
	
	
	multilintolin_eval_data_mp(const multilintolin_eval_data_mp & other) : solver_mp(other)
	{
		init();
		copy(other);
	}
	
	multilintolin_eval_data_mp & operator=(const multilintolin_eval_data_mp & other)
	{
		
		init();
		copy(other);
		return *this;
	}
	
	
	
	~multilintolin_eval_data_mp(){
		clear();
		// no need to reset the counters.
	}
	

	
	void print()
	{
		
	}
	

	
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config){return 0;}
	
	int receive(parallelism_config & mpi_config){return 0;}
	
	
	
	int setup(const multilin_config & config,
						const witness_set & W,
						vec_mp * target_linears,
						solver_configuration & solve_options);
	
	
	
protected:
	
	
	void clear()
	{
		solver_mp::clear();
		
		
		for (int ii=0; ii<num_linears; ii++) {
			clear_vec_mp(current_linear[ii]);
			clear_vec_mp(old_linear[ii]);
		}
		free(current_linear);
		free(old_linear);
		
		
		if (this->MPType == 2) {
			for (int ii=0; ii<num_linears; ii++) {
				clear_vec_mp(current_linear_full_prec[ii]);
				clear_vec_mp(old_linear_full_prec[ii]);
			}
			free(current_linear_full_prec);
			free(old_linear_full_prec);
		}
		
		
		
	} // re: clear
	
	void init();
	
	
	void copy(const multilintolin_eval_data_mp & other)
	{
		solver_mp::copy(other);
		
		
		if (this->num_linears==0) {
			current_linear = (vec_mp *) br_malloc(other.num_linears*sizeof(vec_mp));
			old_linear = (vec_mp *) br_malloc(other.num_linears*sizeof(vec_mp));
		}
		else
		{
			current_linear = (vec_mp *) br_realloc(current_linear, other.num_linears*sizeof(vec_mp));
			old_linear = (vec_mp *) br_realloc(old_linear, other.num_linears*sizeof(vec_mp));
		}
		
		for (int ii=0; ii<other.num_linears; ii++) {
			init_vec_mp(current_linear[ii],0);
			init_vec_mp(old_linear[ii],0);
			vec_cp_mp(current_linear[ii],other.current_linear[ii]);
			vec_cp_mp(old_linear[ii],other.old_linear[ii]);
		}
		
		
		
		if (this->MPType==2) {
			if (this->num_linears==0) {
				current_linear_full_prec = (vec_mp *) br_malloc(other.num_linears*sizeof(vec_mp));
				old_linear_full_prec = (vec_mp *) br_malloc(num_linears*sizeof(vec_mp));
			}
			else
			{
				current_linear_full_prec = (vec_mp *) br_realloc(current_linear_full_prec, other.num_linears*sizeof(vec_mp));
				old_linear_full_prec = (vec_mp *) br_realloc(old_linear_full_prec, other.num_linears*sizeof(vec_mp));
			}
			
			for (int ii=0; ii<other.num_linears; ii++) {
				init_vec_mp2(current_linear_full_prec[ii],0,1024);
				init_vec_mp2(old_linear_full_prec[ii],0,1024);
				vec_cp_mp(current_linear_full_prec[ii],other.current_linear_full_prec[ii]);
				vec_cp_mp(old_linear_full_prec[ii],other.old_linear_full_prec[ii]);
			}
		}

		this->num_linears= other.num_linears;
	} // re: copy
	
	
	
}; // re: class nullspace_eval_data_mp














// the mp version
// this must be defined before the double version, because double has mp.
class multilintolin_eval_data_d : public solver_d
{
public:
	
	multilintolin_eval_data_mp * BED_mp;
	int num_linears;
	
	
	
	vec_d *current_linear;						
	vec_d *old_linear;								
	
	
	
	// default initializer
	multilintolin_eval_data_d() : solver_d(){
		init();
	}
	
	multilintolin_eval_data_d(int mp) : solver_d(mp)
	{
		this->MPType = mp;
		init();
	}
	
	
	multilintolin_eval_data_d(const multilintolin_eval_data_d & other) : solver_d(other)
	{
		init();
		copy(other);
	}
	
	multilintolin_eval_data_d & operator=(const multilintolin_eval_data_d & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	
	
	~multilintolin_eval_data_d(){
		clear();
		// no need to reset the counters.
	}
	
	
	
	void print()
	{
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config){return 0;}
	
	int receive(parallelism_config & mpi_config){return 0;}
	
	
	
	int setup(const multilin_config & config,
						const witness_set & W,
						vec_mp * target_linears,
						solver_configuration & solve_options);
	
	
	
protected:
	void init();
	
	void clear()
	{
		solver_d::clear();
		
		for (int ii=0; ii<num_linears; ii++) {
			clear_vec_d(current_linear[ii]);
			clear_vec_d(old_linear[ii]);
		}
		free(current_linear);
		free(old_linear);
		
		if (this->MPType==2)
		{
			delete this->BED_mp;
		}
	}
};



int multilin_solver_master_entry_point(const witness_set						&W, // carries with it the start points, and the linears.
																			 witness_set							*W_new, // new data goes in here
																			 vec_mp * new_linears,
																			 const multilin_config &		config,
																			 solver_configuration		& solve_options);




//the new custom evaluator for this solver

int multilin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);




int multilin_to_lin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);



int check_issoln_multilintolin_d(endgame_data_t *EG,
														tracker_config_t *T,
														void const *ED);
int check_issoln_multilintolin_mp(endgame_data_t *EG,
														 tracker_config_t *T,
														 void const *ED);



int multilintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);




int change_multilintolin_eval_prec(void const *ED, int new_prec);






































//
//
//
//
///** the main function for finding critical conditions WRT a projection
// */
//int multilintolin_solver_main(int MPType,
//															const witness_set & W,
//															mat_mp randomizer_matrix_full_prec,
//															vec_mp *new_linears_full_prec,
//															witness_set *W_new,
//															solver_configuration & solve_options);
//
//
//
//
//int multilin_to_lin_solver_d(int MPType,
//														 const witness_set & W,  // should include the old linear
//														 mat_mp randomizer_matrix_full_prec,
//														 vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
//														 witness_set *W_new,
//														 solver_configuration & solve_options);
//
//
//
//int multilin_to_lin_solver_mp(int MPType,
//															const witness_set & W,  // includes the initial linear.
//															mat_mp randomizer_matrix_full_prec,  // for randomizing down to N-1 equations.
//															vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
//															witness_set *W_new,
//															solver_configuration & solve_options);
//



void multilin_to_lin_track_d(trackingStats *trackCount,
														 FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
														 const witness_set & W,
														 vec_mp *new_linears_full_prec,
														 post_process_t *endPoints,
														 FILE *FAIL,
														 int pathMod,
														 multilintolin_eval_data_d *ED_d, multilintolin_eval_data_mp *ED_mp,
														 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 int (*change_prec)(void const *, int),
														 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),solver_configuration & solve_options);



int multilin_to_lin_setup_d(FILE **OUT, boost::filesystem::path outName,
														FILE **midOUT, boost::filesystem::path midName,
														multilintolin_eval_data_d *ED,
														prog_t *dummyProg,
														int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
														int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														boost::filesystem::path preprocFile, boost::filesystem::path degreeFile,
														int findStartPts, boost::filesystem::path pointsIN, boost::filesystem::path pointsOUT,
														mat_mp randomizer_matrix_full_prec,
														const witness_set & W,
														solver_configuration & solve_options);





void multilintolin_eval_clear_d(multilintolin_eval_data_d *ED, int clearRegen, int MPType);




void setupmultilintolinEval_d(tracker_config_t *T,
															char preprocFile[], char degreeFile[],
															prog_t *dummyProg,
															multilintolin_eval_data_d *BED,
															mat_mp randomizer_matrix_full_prec,
															const witness_set & W,
															solver_configuration & solve_options);



void cp_multilintolin_eval_data_d(multilintolin_eval_data_d *BED, multilintolin_eval_data_d *BED_d_input, multilintolin_eval_data_mp *BED_mp_input, int MPType);






void multilin_to_lin_track_mp(trackingStats *trackCount,
															FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
															const witness_set & W,
															vec_mp *new_linears_full_prec,
															post_process_t *endPoints,
															FILE *FAIL,
															int pathMod,
															multilintolin_eval_data_mp *ED_d,
															int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
															int (*change_prec)(void const *, int),
															int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
															solver_configuration & solve_options);




int multilin_to_lin_setup_mp(FILE **OUT, boost::filesystem::path outName,
														 FILE **midOUT, boost::filesystem::path midName,
														 multilintolin_eval_data_mp *ED,
														 prog_t *dummyProg,
														 int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
														 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 boost::filesystem::path preprocFile, boost::filesystem::path degreeFile,
														 int findStartPts,
														 boost::filesystem::path pointsIN, boost::filesystem::path pointsOUT,
														 mat_mp randomizer_matrix_full_prec,
														 const witness_set & W,
														 solver_configuration & solve_options);



void multilintolin_eval_clear_mp(multilintolin_eval_data_mp *ED, int clearRegen, int MPType);


void setupmultilintolinEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg,
															 int squareSize, int patchType, int ssType, int prec,
															 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
															 multilintolin_eval_data_mp *BED, int adjustDegrees,
															 mat_mp randomizer_matrix_full_prec,
															 const witness_set & W,
															 solver_configuration & solve_options);

#endif


