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

#include "missing_bertini_headers.hpp"

#include "fileops.hpp"
#include "data_type.hpp"

#include "programConfiguration.hpp"
#include "witnessSet.hpp"


///////////
//
//    SOLVER CONFIGURATION
//
//////////


class solver_configuration : public prog_config
{
public:
	
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
};








/**
 conversion function for turning endgame_data_t into post_process_t.
 */
void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG);

/**
 bertini_real's own method of finding SINGULAR solutions among a collection of post_process_t's.
 this singular detection does NOT declare a point singular based on multiplicity -- just on the basis of condition number.
 */
int BRfindSingularSolns(post_process_t *endPoints,
												int num_sols, int num_vars,
												tracker_config_t *T );

/**
 bertini_real's own method of finding FINITE solutions among a collection of post_process_t's.
 */
int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
											tracker_config_t *T );

/**
 determines if a solution is allowable.  based on solver_configuration.
 */
int is_acceptable_solution(post_process_t endPoint, solver_configuration & solve_options);


/**
 bertini_real's version of post-processing.  options are set via the solver_configuration.
 */
void BRpostProcessing(post_process_t *endPoints, witness_set *W_new, int num_pts,
											preproc_data *preProcData, tracker_config_t *T,
											solver_configuration & solve_options);











void call_for_help(int solver_type, parallelism_config & mpi_config);



class solver
{
public:
	
	boost::filesystem::path preproc_file;
	boost::filesystem::path function_file;
	
	int num_steps; ///< the number of evaluations made using this evaluator
	int verbose_level;  ///< how verbose to be
	int MPType; ///< the multiple precision type for solve
	preproc_data preProcData; ///< information related to the SLP for system
	prog_t *SLP; ///< the SLP
	int num_variables; ///< the grand total number of variables to be computed on in the problem at hand 
	
	// these function handles are required for all solvers.
	int (*evaluator_function_d) (point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *); // function handle to evaluator to use
	int (*evaluator_function_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
	
	int (*precision_changer)(void const *ED, int new_prec);
	
	int (*dehomogenizer)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *);
	
	int (*is_solution_checker_d)(endgame_data_t *, tracker_config_t *, void const *);
	int (*is_solution_checker_mp)(endgame_data_t *, tracker_config_t *, void const *);
	
	solver(){
		init();
	}
	
	
	solver(int new_mp_type){
		init();
		this->MPType = new_mp_type;
	}
	
	void copy(const solver & other){
		cp_preproc_data(&this->preProcData, other.preProcData);
		cp_prog_t(this->SLP, other.SLP);
		
		this->is_solution_checker_d = other.is_solution_checker_d;
		this->is_solution_checker_mp = other.is_solution_checker_mp;
		
		this->evaluator_function_mp = other.evaluator_function_mp;
		this->evaluator_function_d = other.evaluator_function_d;
		this->precision_changer = other.precision_changer;
		this->dehomogenizer = other.dehomogenizer;
		
		this->num_variables = other.num_variables;
		this->MPType = other.MPType;
		this->verbose_level = other.verbose_level;
		this->num_steps = other.num_steps;
	}
	
	int send();
	
	int receive();
	
	void setup(prog_t * _SLP)
	{
		setupPreProcData(const_cast<char *>(preproc_file.c_str()), &this->preProcData);
		this->SLP = _SLP; // assign the pointer
	}
	
private:
	void init()
	{
		
		this->preproc_file = "preproc_data";
		this->function_file = "func_input";
		
		this->MPType = 2;
		SLP = NULL;
		num_variables = 0;
		
		this->num_steps = 0;
		this->verbose_level = 0;
		
		// initialize the function handles.
		is_solution_checker_d = NULL;
		is_solution_checker_mp = NULL;
		evaluator_function_d = NULL;
		evaluator_function_mp = NULL;
		precision_changer = NULL;
		dehomogenizer = NULL;
	}
	// these virtual functions will need to be programmed into the derived classes.
	
	

};  // re: generic solver BASE class









class solver_mp : public solver
{
	
public:
	
	patch_eval_data_mp patch; ///< patch in x
  
	comp_mp gamma;    ///< randomizer
	mpq_t *gamma_rat; ///< randomizer
	mat_mp randomizer_matrix;     ///< randomizer
	mat_mp randomizer_matrix_full_prec;  ///< randomizer

	int curr_prec;

	
	solver_mp() : solver(){
		init();
	}; // re: default constructor
	
	
	solver_mp(int new_mp_type) : solver(new_mp_type){
		this->MPType = new_mp_type;
		init();
	}; // re: default constructor
	

	void init(){
		init_mp(gamma);
		init_mat_mp(randomizer_matrix,0,0);
		
		if (this->MPType == 2 ) {
			gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
			init_rat(gamma_rat);
			init_mat_mp2(randomizer_matrix_full_prec,0,0,1024);
		}
	}
	
	~solver_mp(){
		
		clear_mat_mp(randomizer_matrix);
		clear_mp(gamma);
		
		if (MPType==2) {
			clear_mat_mp(randomizer_matrix_full_prec);
			clear_rat(gamma_rat);
		}
	}// re: default destructor
	
	
	solver_mp & operator=( const solver_mp & other)
	{
		solver::operator= (other);
		
		copy(other);
		return *this;
	}  // re: assigment
	
	
	solver_mp(const solver_mp & other) : solver(other){
		copy(other);
	} // re: copy
	
	int send();
	
	
	int receive();
	
	void setup(prog_t * _SLP)
	{
		solver::setup(_SLP);
	}
	
private:
	
	void copy(const solver_mp & other){
		cp_patch_mp(&this->patch, other.patch);
		
		set_mp(this->gamma, other.gamma);
		set_rat(this->gamma_rat, other.gamma_rat);
		
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		mat_cp_mp(this->randomizer_matrix_full_prec, other.randomizer_matrix_full_prec);
		
		this->curr_prec = other.curr_prec;
	}
	

};








class solver_d : public solver
{
	
public:
	
	patch_eval_data_d patch; ///< patch in x
  
	comp_d gamma;    ///< randomizer
	mat_d randomizer_matrix;     ///< randomizer

	
	solver_d() : solver(){
		init_d(gamma);
		init_mat_d(randomizer_matrix,0,0);

	}; // re: default constructor
	
	
	solver_d(int new_mp_type) : solver(new_mp_type){
		init_d(gamma);
		init_mat_d(randomizer_matrix,0,0);
	}; // re: default constructor
	
	
	
	~solver_d(){
		clear_mat_d(randomizer_matrix);
		clear_d(gamma);
	}// re: default destructor
	
	
	solver_d & operator=( const solver_d & other)
	{
		solver::operator= (other);

		copy(other);
		return *this;
	}  // re: assigment
	
	
	solver_d(const solver_d & other) : solver(other){
		copy(other);
	} // re: copy
	
	int send();
	
	
	int receive();
	
	void setup(prog_t * _SLP)
	{
		solver::setup(_SLP);
	}
	
private:
	
	void copy(const solver_d & other){
		
		cp_patch_d(&this->patch, other.patch);
		set_d(this->gamma, other.gamma);
		mat_cp_d(this->randomizer_matrix, other.randomizer_matrix);
	}
	

	
};






/**
 reads in projection from file if user specified, creates one otherwise.
 --
 // currently defaults to create a random real projection with homogeneous value 0;
 */
void get_projection(vec_mp *pi,
										BR_configuration program_options,
										solver_configuration solve_options,
										int num_vars,
										int num_projections);









void generic_set_start_pts(point_data_d ** startPts,
													 witness_set & W);

void generic_set_start_pts(point_data_mp ** startPts,
													 witness_set & W);

void generic_setup_patch(patch_eval_data_d *P, const witness_set & W); // for mp type 0
void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W);// for my type 2
void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W, int prec); // for mp type 1

int generic_setup_files(FILE ** OUT, boost::filesystem::path outname,
												FILE ** MIDOUT, boost::filesystem::path midname);


/** reads the tracker_config_t from file. */
void get_tracker_config(solver_configuration &solve_options,int MPType);
void solver_init_config(solver_configuration &options);
void solver_clear_config(solver_configuration &options);


void generic_tracker_loop(trackingStats *trackCount,
													FILE * OUT, FILE * midOUT,
													witness_set & W,  // was the startpts file pointer.
													post_process_t *endPoints,
													solver_d * ED_d, solver_mp * ED_mp,
													solver_configuration & solve_options);


void generic_tracker_loop_worker(trackingStats *trackCount,
																 FILE * OUT, FILE * MIDOUT,
																 solver_d * ED_d, solver_mp * ED_mp,
																 solver_configuration & solve_options);



void generic_track_path(int pathNum, endgame_data_t *EG_out,
												point_data_d *Pin, point_data_mp *Pin_mp,
												FILE *OUT, FILE *MIDOUT,
												tracker_config_t *T,
												void const *ED_d, void const *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));




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




#endif