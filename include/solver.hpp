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


#include "postProcessing.hpp"
///////////
//
//    SOLVER CONFIGURATION
//
//////////


extern _comp_d  **mem_d;
extern _comp_mp **mem_mp;
extern int *size_d;  // size of mem_d
extern int *size_mp;  // size of mem_mp
extern int *mem_needs_init_d; // determine if mem_d has been initialized
extern int *mem_needs_init_mp; // determine if mem_mp has been initialized


class SLP_global_pointers
{
public:
	
	void capture_globals()
	{
		
		local_mem_d = mem_d;
		local_mem_mp = mem_mp;
		local_size_d = size_d;  // size of mem_d
		local_size_mp = size_mp;  // size of mem_mp
		local_mem_needs_init_d = mem_needs_init_d; // determine if mem_d has been initialized
		local_mem_needs_init_mp = mem_needs_init_mp; // determine if mem_mp has been initialized
	}
	
	
	void set_globals_to_this()
	{
		mem_d							= local_mem_d;
		mem_mp						= local_mem_mp;
		size_d						= local_size_d;  // size of mem_d
		size_mp						= local_size_mp;  // size of mem_mp
		mem_needs_init_d	= local_mem_needs_init_d; // determine if mem_d has been initialized
		mem_needs_init_mp = local_mem_needs_init_mp; // determine if mem_mp has been initialized
	}
	
	void set_globals_null()
	{
		mem_d							= NULL;
		mem_mp						= NULL;
		size_d						= NULL;  // size of mem_d
		size_mp						= NULL;  // size of mem_mp
		mem_needs_init_d	= NULL; // determine if mem_d has been initialized
		mem_needs_init_mp = NULL; // determine if mem_mp has been initialized
	}
	
	
	
	
	SLP_global_pointers()
	{
		init();
	}
	
	
	
	SLP_global_pointers & operator=( const SLP_global_pointers & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	SLP_global_pointers(const SLP_global_pointers & other)
	{
		init();
		copy(other);
	} // re: copy
	
	

	
	_comp_d  **local_mem_d;
	_comp_mp **local_mem_mp;
	int *local_size_d;  // size of mem_d
	int *local_size_mp;  // size of mem_mp
	int *local_mem_needs_init_d; // determine if mem_d has been initialized
	int *local_mem_needs_init_mp; // determine if mem_mp has been initialized
	
	protected:
	
	void init()
	{
		local_mem_d = NULL;
		local_mem_mp = NULL;
		local_size_d = NULL;  // size of mem_d
		local_size_mp = NULL;  // size of mem_mp
		local_mem_needs_init_d = NULL; // determine if mem_d has been initialized
		local_mem_needs_init_mp = NULL; // determine if mem_mp has been initialized
	}
	
	void copy(const SLP_global_pointers & other)
	{
		this->local_mem_d = other.local_mem_d;
		this->local_mem_mp = other.local_mem_mp;
		this->local_size_d = other.local_size_d;  // size of mem_d
		this->local_size_mp = other.local_size_mp;  // size of mem_mp
		this->local_mem_needs_init_d = other.local_mem_needs_init_d; // determine if mem_d has been initialized
		this->local_mem_needs_init_mp = other.local_mem_needs_init_mp; // determine if mem_mp has been initialized
	}
	
};



class solver_configuration : public prog_config
{
public:
	
	tracker_config_t T;
	tracker_config_t T_orig;
	preproc_data PPD;
	
	int allow_multiplicity;
	int allow_singular;
	int allow_infinite;
	int allow_unsuccess;
	
	int path_number_modulus;
	
	int verbose_level;
	int show_status_summary;
	
	int use_midpoint_checker;
	double midpoint_tol;
	
	int use_gamma_trick;
	
	int complete_witness_set;
	
	void reset_tracker_config()
	{
		cp_tracker_config_t(&T,&T_orig);
	}
	
	void backup_tracker_config()
	{
		cp_tracker_config_t(&T_orig,&T);
	}
	
	
	~solver_configuration(){
//		tracker_config_clear(&this->T);
//		tracker_config_clear(&this->T_orig);
//		preproc_data_clear(&this->PPD);
	}
	
	
	void get_PPD()
	{
		parse_preproc_data("preproc_data", &this->PPD);
	}

};













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
	bool have_SLP;
	
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
	
	virtual ~solver()
	{
		freeEvalProg(this->MPType);
	}
	
	virtual int send(parallelism_config & mpi_config);
	
	virtual int receive(parallelism_config & mpi_config);
	
	void setup(prog_t * _SLP)
	{
		setupPreProcData(const_cast<char *>(preproc_file.c_str()), &this->preProcData);
		this->SLP = _SLP; // assign the pointer
		this->have_SLP = true;
	}
	
	void setup()
	{
		setupPreProcData(const_cast<char *>(preproc_file.c_str()), &this->preProcData);
	}
	
protected:
	
	
	
	void init()
	{
		
		this->preproc_file = "preproc_data";
		this->function_file = "func_input";
		
		this->MPType = 2;
		SLP = NULL;
		have_SLP = false;
		
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
	// these referenced functions will need to be programmed into the derived classes.
	
	
	void copy(const solver & other){
		cp_preproc_data(&this->preProcData, other.preProcData);
		
		this->SLP = other.SLP;

		
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
		this->MPType = 2;
		init();
	}; // re: default constructor
	
	
	solver_mp(int new_mp_type) : solver(new_mp_type){
		this->MPType = new_mp_type;
		init();
	}; // re: default constructor
	

	void init(){
		solver::init();
		
		init_mp(gamma);
		init_mat_mp(randomizer_matrix,0,0);
		
		if (this->MPType == 2 ) {
			gamma_rat = (mpq_t *)br_malloc(2 * sizeof(mpq_t));
			init_rat(gamma_rat);
			init_mat_mp2(randomizer_matrix_full_prec,0,0,1024);
		}
	}
	
	virtual ~solver_mp(){
		
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
	
	virtual int send(parallelism_config & mpi_config);
	
	
	virtual int receive(parallelism_config & mpi_config);
	
	void setup(prog_t * _SLP)
	{
		solver::setup(_SLP);
	}
	
	void setup()
	{
		solver::setup();
	}
protected:
	
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
		init();
	}; // re: default constructor
	
	
	solver_d(int new_mp_type) : solver(new_mp_type){
		init();
	}; // re: default constructor
	
	
	
	virtual ~solver_d(){
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
	
	virtual int send(parallelism_config & mpi_config);
	
	
	virtual int receive(parallelism_config & mpi_config);
	
	void setup(prog_t * _SLP)
	{
		solver::setup(_SLP);
	}
	
	void setup()
	{
		solver::setup();
	}
protected:
	
	
	void init()
	{
		solver::init();
		
		init_d(gamma);
		init_mat_d(randomizer_matrix,0,0);
	}
	
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
										const solver_configuration & solve_options,
										int num_vars,
										int num_projections);




void adjust_tracker_AMP(tracker_config_t * T, int num_variables);




void generic_set_start_pts(point_data_d ** startPts,
													 const witness_set & W);

void generic_set_start_pts(point_data_mp ** startPts,
													 const witness_set & W);

void generic_setup_patch(patch_eval_data_d *P, const witness_set & W); // for mp type 0
void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W);// for my type 2
void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W, int prec); // for mp type 1

int generic_setup_files(FILE ** OUT, boost::filesystem::path outname,
												FILE ** MIDOUT, boost::filesystem::path midname);


/** reads the tracker_config_t from file. */
void get_tracker_config(solver_configuration &solve_options,int MPType);
void solver_init_config(solver_configuration &options);
void solver_clear_config(solver_configuration &options);

void generic_solver_master(witness_set * W_new, witness_set & W,
													 solver_d * ED_d, solver_mp * ED_mp,
													 solver_configuration & solve_options);

void generic_tracker_loop(trackingStats *trackCount,
													FILE * OUT, FILE * midOUT,
													witness_set & W,  // was the startpts file pointer.
													post_process_t *endPoints,
													solver_d * ED_d, solver_mp * ED_mp,
													solver_configuration & solve_options);

void generic_tracker_loop_master(trackingStats *trackCount,
																 FILE * OUT, FILE * MIDOUT,
																 witness_set & W,  // was the startpts file pointer.
																 post_process_t *endPoints,
																 solver_d * ED_d, solver_mp * ED_mp,
																 solver_configuration & solve_options);

void send_start_points(int next_worker, int num_packets,
											point_data_d *startPts_d,
											point_data_mp *startPts_mp,
											int & next_index,
											solver_configuration & solve_options);

int receive_endpoints(trackingStats *trackCount,
											endgame_data_t *EG_receives, int & max_incoming,
											int & solution_counter,
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