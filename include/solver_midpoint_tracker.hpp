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


#ifndef SOLVER_MIDPOINT_H
#define SOLVER_MIDPOINT_H

#include "missing_bertini_headers.hpp"



#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"

#include "determinant_derivative.hpp"
#include "programConfiguration.hpp"
#include "postProcessing.hpp"
#include "witnessSet.hpp"

class midpoint_config
{
	
public:
	
	mat_mp randomizer_matrix;
	
	
	int num_x_vars;
	int num_y_vars;
	int num_z_vars;
	
	prog_t *SLP_crit;
	
	vec_mp *pi;
	int num_projections;
	
	//patch already lives in the base class.
	
	comp_mp v_target;
	comp_mp u_target;
	
	mat_mp randomizer_matrix_crit;
	
	comp_mp crit_val_left;
	comp_mp crit_val_right;
	
	midpoint_config(){
		
	}
	
	void add_projection(vec_mp proj)
	{
		if (num_projections==0) {
			pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else
		{
			pi = (vec_mp *) br_realloc(pi,(num_projections+1)*sizeof(vec_mp));
		}
		
		init_vec_mp2(pi[num_projections],0,1024);
		vec_cp_mp(pi[num_projections],proj);
		
		num_projections++;
	}
	
	void copy(const midpoint_config & ns_config){
		
	}
	
	
	
	void init(){
		num_projections=0;
		pi = NULL;
		
		init_mp2(v_target,1024);
		init_mp2(u_target,1024);
		
		init_mp2(crit_val_left,1024);
		init_mp2(crit_val_right,1024);
		
		
		init_mat_mp2(randomizer_matrix,0,0,1024);
		init_mat_mp2(randomizer_matrix_crit,0,0,1024);
		
		this->SLP_crit = NULL;
	}
	
	void clear()
	{
		clear_mp(v_target);
		clear_mp(u_target);
		clear_mp(crit_val_left);
		clear_mp(crit_val_right);
		
		clear_mat_mp(randomizer_matrix);
		clear_mat_mp(randomizer_matrix_crit);
		
		for (int ii=0; ii<num_projections; ii++) {
			clear_vec_mp(pi[ii]);
		}
		free(pi);
		
//		clear_prog_t(SLP_crit);
	}
	

	
	void print()
	{
		
	}
};






// the mp version
// this must be defined before the double version, because double has mp.
class midpoint_eval_data_mp : public solver_mp
{
public:
	
	int num_x_vars;
	int num_y_vars;
	int num_z_vars;
	
	prog_t *SLP_crit;
	
	vec_mp *pi;
	int num_projections;
	
	//patch already lives in the base class.
	
	comp_mp v_target;
	comp_mp u_target;
	
	mat_mp randomizer_matrix_crit;
	
	comp_mp crit_val_left;
	comp_mp crit_val_right;
	
	
	// default initializer
	midpoint_eval_data_mp() : solver_mp(){
		
		reset_counters();
		
		init();
	}
	
	midpoint_eval_data_mp(int mp) : solver_mp(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	~midpoint_eval_data_mp(){
		
		midpoint_eval_data_mp::clear();
		// no need to reset the counters.
	}
	
	void print()
	{

	};
	
	void reset_counters()
	{

	} // re: reset_counters
	
	
	midpoint_eval_data_mp & operator=(const midpoint_eval_data_mp & other)
	{
		clear(); // this is wasteful, but safe for now
		reset_counters();// this is wasteful, but safe for now
		
		copy(other);
		return *this;
	}
	
	midpoint_eval_data_mp(const midpoint_eval_data_mp & other)
	{
		
		
		solver_mp();
		midpoint_eval_data_mp();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(prog_t * _SLP,
						midpoint_config *ns_config,
						witness_set & W,
						solver_configuration & solve_options);
	
	
	
private:
	
	
	void clear()
	{
				
	} // re: clear
	
	void init();
	
	void copy(const midpoint_eval_data_mp & other)
	{
		
		clear();
		reset_counters();
		
	} // re: copy
	
	
	
	
}; // re: class nullspace_eval_data_mp










// the double version
// this must be defined after the mp version, because double has mp.
class midpoint_eval_data_d : public solver_d
{
public:
	
	midpoint_eval_data_mp *BED_mp; // used only for AMP
	
	int num_x_vars;
	int num_y_vars;
	int num_z_vars;
	
	prog_t *SLP_crit;
	
	vec_d *pi;
	int num_projections;
	
	//patch already lives in the base class.
	
	comp_d v_target;
	comp_d u_target;
	
	mat_d randomizer_matrix_crit;
	
	comp_d crit_val_left;
	comp_d crit_val_right;
	
	
	// default initializer
	midpoint_eval_data_d() : solver_d(){
		
		reset_counters();
		
		init();
	}
	
	midpoint_eval_data_d(int mp) : solver_d(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	~midpoint_eval_data_d(){
		
		midpoint_eval_data_d::clear();
		// no need to reset the counters.
	}
	
	void print()
	{
		
	};
	
	void reset_counters()
	{
		
	} // re: reset_counters
	
	
	midpoint_eval_data_d & operator=(const midpoint_eval_data_d & other)
	{
		clear(); // this is wasteful, but safe for now
		reset_counters();// this is wasteful, but safe for now
		
		copy(other);
		return *this;
	}
	
	midpoint_eval_data_d(const midpoint_eval_data_d & other)
	{
		init();
		solver_d();
		midpoint_eval_data_d();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(prog_t * _SLP,
						midpoint_config *ns_config,
						witness_set & W,
						solver_configuration & solve_options);
	
	
private:
	
	
	void clear()
	{
		
		
		
		
		if (this->MPType==2) {
			
			delete this->BED_mp;
		}
	} // re: clear
	
	void init();
	
	void copy(const midpoint_eval_data_d & other)
	{
		
		clear();
		reset_counters();
		
				
	} // re: copy
	
};







/** the main function for finding critical conditions WRT a projection
 */

int midpoint_solver_master_entry_point(int										MPType,
																			 witness_set						& W, // carries with it the start points, and the linears.
																			 witness_set						*W_new, // new data goes in here
																			 midpoint_config				*ns_config,
																			 solver_configuration		& solve_options);



int midpoint_solver(int MPType,
												witness_set & W,
												witness_set *W_new,
												midpoint_config				*ns_config,
												solver_configuration & solve_options);






int midpoint_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED);

int midpoint_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);





int midpoint_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);



int change_midpoint_eval_prec(void const *ED, int new_prec);



int check_issoln_midpoint_d(endgame_data_t *EG,
																tracker_config_t *T,
																void const *ED);
int check_issoln_midpoint_mp(endgame_data_t *EG,
																 tracker_config_t *T,
																 void const *ED);

int check_isstart_midpoint_d(point_d testpoint,
																 tracker_config_t *T,
																 void const *ED);


void check_midpoint_evaluator(point_mp current_values,
															 void const *ED);




void midpoint_slave_entry_point(solver_configuration & solve_options);










#endif


