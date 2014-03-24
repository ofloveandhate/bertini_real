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

#include "bertini_headers.hpp"



#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"
#include "programConfiguration.hpp"
#include "postProcessing.hpp"







class multilin_config
{
	
public:
	
	bool have_randomizer;
	system_randomizer *randomizer;
	
	
	SLP_global_pointers SLP_memory;
	prog_t * SLP;
	
	int MPType;
	
	bool have_mem;
	
	
	multilin_config(solver_configuration & solve_options,
					const witness_set & W)
	{
		init();
		
		set_memory(solve_options);
		
		make_randomizer(solve_options, W);
	}
	
	
	
	multilin_config(solver_configuration & solve_options,
					system_randomizer * _random)
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
	
	
	multilin_config(system_randomizer * _random)
	{
		init();
		set_randomizer(_random);
	}
	
	void make_randomizer(const solver_configuration & solve_options, const witness_set & W)
	{
		
		randomizer = new system_randomizer;
		have_randomizer = true;
		this->randomizer->setup(W.num_variables-W.num_linears-W.num_patches, solve_options.PPD.num_funcs);
		
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
	
	
	
	void set_randomizer(system_randomizer * _random)
	{
		randomizer = _random;
	}
	
	
	
	
	void copy(const multilin_config & other)
	{
		
		//ah yes, this problem.
		//		this->SLP = other.SLP;// this needs to be a deep copy
	}
	
	void init()
	{

		SLP = new prog_t;
		have_mem = false;
		
		have_randomizer = false;
	}
	
	
	void clear()
	{
		
		if (have_randomizer) {
			delete randomizer;
		}
		
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
		std::cout << "instantiating multilin eval d without declaring mp type"	 << std::endl;
		mypause();
		init();
	}
	
	multilintolin_eval_data_mp(int mp) : solver_mp(mp){
		this->MPType = mp;
		init();
	}
	
	
	multilintolin_eval_data_mp(const multilintolin_eval_data_mp & other) : solver_mp(other)
	{
		this->MPType = other.MPType;
		init();
		copy(other);
	}
	
	multilintolin_eval_data_mp & operator=(const multilintolin_eval_data_mp & other)
	{
		this->MPType = other.MPType;
		copy(other);
		return *this;
	}
	
	
	
	virtual ~multilintolin_eval_data_mp(){
		multilintolin_eval_data_mp::clear();
	}
	
	
	
	virtual void print()
	{
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(const multilin_config & config,
			  const witness_set & W,
			  vec_mp * target_linears,
			  solver_configuration & solve_options);
	
	
	
protected:
	
	
	void clear()
	{
		
		
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
		std::cout << "instantiating multilin eval d without declaring mp type"	 << std::endl;
		mypause();
		init();
	}
	
	multilintolin_eval_data_d(int mp) : solver_d(mp)
	{
		this->MPType = mp;
		
		init();
	}
	
	
	multilintolin_eval_data_d(const multilintolin_eval_data_d & other) : solver_d(other)
	{
		this->MPType = other.MPType;
		init();
		copy(other);
	}
	
	multilintolin_eval_data_d & operator=(const multilintolin_eval_data_d & other)
	{
		this->MPType = other.MPType;
		copy(other);
		return *this;
	}
	
	
	
	virtual ~multilintolin_eval_data_d(){
		clear();
	}
	
	
	
	virtual void print()
	{
		solver_d::print();
		
		std::cout << "multilintolin evaluator data (double):" << std::endl;
		for (int ii=0; ii<num_linears; ii++) {
			
			std::stringstream name;
			name << "old_linear_" << ii;
			print_point_to_screen_matlab(old_linear[ii],name.str());
		}
		
		for (int ii=0; ii<num_linears; ii++) {
			
			std::stringstream name;
			name << "current_linear_" << ii;
			print_point_to_screen_matlab(current_linear[ii],name.str());
		}
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(const multilin_config & config,
			  const witness_set & W,
			  vec_mp * target_linears,
			  solver_configuration & solve_options);
	
	
	
protected:
	void init();
	
	void clear()
	{
		
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



int multilin_solver_master_entry_point(const witness_set & W, // carries with it the start points, and the linears.
									   solver_output & solve_out, // new data goes in here
									   vec_mp * new_linears,
									   const multilin_config &		config,
									   solver_configuration		& solve_options);



int multilin_slave_entry_point(solver_configuration & solve_options);


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









#endif


