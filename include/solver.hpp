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

#include "programConfiguration.hpp"
//#include "postProcessing.hpp"
#include "witnessSet.hpp"
#include "missing_bertini_headers.hpp"
//#include "programConfiguration.hpp"


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


class solver_mp : public solver
{
	
public:
	
	patch_eval_data_mp patch; // patch in x
  preproc_data preProcData; // information related to the SLP for system
	prog_t *SLP; // the SLP
	mpq_t *gamma_rat; // randomizer
	comp_mp gamma;    // randomizer
	mat_mp randomizer_matrix;     // randomizer
	mat_mp randomizer_matrix_full_prec;  // randomizer
	int num_variables;
	FILE *FOUT;
	int num_steps;
	int curr_prec;
	int verbose_level;
	
	
	
	
	
	
	solver_mp() : solver(){
		SLP = NULL;
		
		init_mp(gamma);
		init_rat(gamma_rat);
		
		init_mat_mp(randomizer_matrix,0,0);
		init_mat_mp2(randomizer_matrix_full_prec,0,0,1024);
		
		num_variables = -1;
		
		
		std::cout << "initialized solver_mp" << std::endl;
	}; // re: default constructor
	
	
	
	
	~solver_mp(){
		clear_mat_mp(randomizer_matrix);
		clear_mat_mp(randomizer_matrix_full_prec);
		
		clear_mp(gamma);
		clear_rat(gamma_rat);
	}// re: default destructor
	
	
	solver_mp & operator=( const solver_mp & other)
	{
		cp_patch_mp(&this->patch, other.patch);
		
		cp_preproc_data(&this->preProcData, other.preProcData);
		
		cp_prog_t(this->SLP, other.SLP);
		
		set_mp(this->gamma, other.gamma);
		set_rat(this->gamma_rat, other.gamma_rat);
		
		
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		mat_cp_mp(this->randomizer_matrix_full_prec, other.randomizer_matrix_full_prec);
		
		this->num_variables = other.num_variables;
		
		this->FOUT = other.FOUT;
		this->num_steps = other.num_steps;
		this->curr_prec = other.curr_prec;
		this->verbose_level = other.verbose_level;
		return *this;
	}  // re: assigment
	
	
	solver_mp(const solver_mp & other){
		cp_patch_mp(&this->patch, other.patch);
		
		cp_preproc_data(&this->preProcData, other.preProcData);
		
		cp_prog_t(this->SLP, other.SLP);
		
		set_mp(this->gamma, other.gamma);
		set_rat(this->gamma_rat, other.gamma_rat);
		
		
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		mat_cp_mp(this->randomizer_matrix_full_prec, other.randomizer_matrix_full_prec);
		
		this->num_variables = other.num_variables;
		
		this->FOUT = other.FOUT;
		this->num_steps = other.num_steps;
		this->curr_prec = other.curr_prec;
		this->verbose_level = other.verbose_level;
	} // re: copy
	
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




#endif