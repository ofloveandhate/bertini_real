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




#include "bertini_headers.hpp"



#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"

#include "programConfiguration.hpp"
#include "postProcessing.hpp"


#ifndef SOLVER_NULLSPACE_H
#define SOLVER_NULLSPACE_H




class nullspace_config
{
	
public:
	int target_dim;   ///< r			the dimension of the real set we are looking for
	int ambient_dim;  ///< k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	int num_jac_equations;    // the number of equations (after the randomization)
	int num_projections;
	
	int num_v_vars;  ///< N   number of variables in original problem statement (including homogenizing variables)
	int num_natural_vars;  ///< N-k+\ell
	int num_synth_vars;
	
	
	int max_degree;						///< the max degree of differentiated (randomized) functions
	
	system_randomizer * randomizer;
	
	bool randomize_lower;
	mat_mp lower_randomizer;
	
	vec_mp **starting_linears;	///< outer layer should have as many as there are randomized equations (N-k)
								///< inside layer has number corresponding to max of randomized_degrees
	
	int num_additional_linears; ///<
	vec_mp *additional_linears_terminal; ///<
	vec_mp *additional_linears_starting; ///<
	
	
	int num_v_linears;   ///< # is  (N), cause there are N equations in the subsystem.
	vec_mp *v_linears;    ///<
	
	vec_mp v_patch; ///< length of this should be N-k+\ell
	
	vec_mp *target_projection; ///< # of these should be \ell

	bool numerical_derivative;
	
	void clear();
	
	nullspace_config(){
		init();
	}
	
	void print();
	
	~nullspace_config()
	{
		clear();
	}
	
	
private:
	
	void init()
	{
		target_dim = ambient_dim = target_crit_codim = num_jac_equations = -1;
		num_natural_vars = num_v_vars = num_synth_vars = -1;
		max_degree = -1;
		numerical_derivative = false;
		
		starting_linears = NULL;
		
		additional_linears_starting = additional_linears_terminal = NULL;
		
		num_v_linears = -1;
		v_linears = NULL;
		
		target_projection = NULL;
	}
};






// the mp version
// this must be defined before the double version, because double has mp.
class nullspacejac_eval_data_mp : public solver_mp
{
public:
	
	prog_deriv_t * SLP_derivative;
	
	
	int num_jac_equations;
	int target_dim;   // r			the dimension of the real set we are looking for
	int ambient_dim;  // k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	
	int num_natural_vars;  // N   number of variables in original problem statement (including homogenizing variables)
	int num_v_vars;  // (N-k) + (k-\ell+1)
	int num_synth_vars;
	
	
	int max_degree;						// the max degree of differentiated (randomized) functions

	
	int num_additional_linears;
	vec_mp *additional_linears_terminal;
	vec_mp *additional_linears_terminal_full_prec;
	
	vec_mp *additional_linears_starting;
	vec_mp *additional_linears_starting_full_prec;
	
	
	
	
	
	
	vec_mp **starting_linears; // outer layer should have as many as there are randomized equations
							   // inside layer has number corresponding to randomized_degrees
	vec_mp **starting_linears_full_prec; // outer layer should have as many as there are randomized equations
										 // inside layer has number corresponding to randomized_degrees
	
	
	int num_v_linears;
	vec_mp *v_linears;         // should be as many in here as there are randomized equations
	vec_mp *v_linears_full_prec;         // should be as many in here as there are randomized equations
	
	vec_mp v_patch;
	vec_mp v_patch_full_prec;
	
	bool randomize_lower;
	mat_mp lower_randomizer;
	mat_mp lower_randomizer_full_prec;
	
	mat_mp jac_with_proj;
	mat_mp jac_with_proj_full_prec;
	
	
	
	int num_projections;
	vec_mp *target_projection; // # of these should be target_dim (for now)
	vec_mp *target_projection_full_prec; // # of these should be target_dim (for now)
	
	
	// default initializer
	nullspacejac_eval_data_mp() : solver_mp(){
		
		reset_counters();
		
		init();
	}
	
	nullspacejac_eval_data_mp(int mp) : solver_mp(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	virtual ~nullspacejac_eval_data_mp(){
		
		nullspacejac_eval_data_mp::clear();
		// no need to reset the counters.
	}
	
	void print();
    
    
    
    
	void reset_counters()
	{
		num_jac_equations = 0;
		target_dim = 0;   // r			the dimension of the real set we are looking for
		ambient_dim = 0;  // k			the dimension of the complex component we are looking IN.
		target_crit_codim = 0;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		num_v_vars = -771;  // N   number of variables in original problem statement (including homogenizing variables)
		num_natural_vars = -791;  // N-k+\ell
		num_synth_vars = -1341;
		
		
		max_degree = 0;						// the max degree of differentiated (randomized) functions
		
		
		num_projections = 0;
		num_v_linears = 0;
		num_additional_linears = 0;
	} // re: reset_counters
	
	
	nullspacejac_eval_data_mp & operator=(const nullspacejac_eval_data_mp & other)
	{
		clear(); // this is wasteful, but safe for now
		reset_counters();// this is wasteful, but safe for now
		
		copy(other);
		return *this;
	}
	
	nullspacejac_eval_data_mp(const nullspacejac_eval_data_mp & other)
	{
		
		
		solver_mp();
		nullspacejac_eval_data_mp();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(prog_t * _SLP,
			  nullspace_config *ns_config,
			  witness_set & W,
			  solver_configuration & solve_options);
	
	
	
private:
	
	
	void clear()
	{
		
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_mp(additional_linears_terminal[ii]);
			}
			free(additional_linears_terminal);
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_mp(additional_linears_starting[ii]);
			}
			free(additional_linears_starting);
		}
		
		
		
		if (num_jac_equations>0) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					clear_vec_mp(starting_linears[ii][jj]);
				}
				free(starting_linears[ii]);
			}
			free(starting_linears);
		}
		
		if (num_v_linears>0) {
			for (int ii=0; ii<num_v_linears; ii++) {
				clear_vec_mp(v_linears[ii]);
			}
			free(v_linears);
		}
		
		
		clear_vec_mp(v_patch);
		
		clear_mat_mp(jac_with_proj);
		
		clear_mat_mp(lower_randomizer);
		

		
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_mp(target_projection[ii]);
			}
			free(target_projection);
		}
		
		
		
		if (this->MPType==2) {
			
			clear_mat_mp(lower_randomizer_full_prec);
			
			if (num_additional_linears>0) {
				for (int ii=0; ii<num_additional_linears; ii++) {
					clear_vec_mp(additional_linears_terminal_full_prec[ii]);
				}
				free(additional_linears_terminal_full_prec);
				
				for (int ii=0; ii<num_additional_linears; ii++) {
					clear_vec_mp(additional_linears_starting_full_prec[ii]);
				}
				free(additional_linears_starting_full_prec);
			}
			
			
			
			
			for (int ii=0; ii<num_jac_equations; ii++){
				for (int jj=0; jj<max_degree; jj++) {
					clear_vec_mp(starting_linears_full_prec[ii][jj]);
				}
				free(starting_linears_full_prec[ii]);
			}
			free(starting_linears_full_prec);
			
			for (int ii=0; ii<num_v_linears; ii++)
				clear_vec_mp(v_linears_full_prec[ii]);
			
			free(v_linears_full_prec);
			
			clear_vec_mp(v_patch_full_prec);
			clear_mat_mp(jac_with_proj_full_prec);

			
			
			if (num_projections>0) {
				for (int ii=0; ii<num_projections; ii++) {
					clear_vec_mp(target_projection_full_prec[ii]);
				}
				free(target_projection_full_prec);
			}
		}
		
		clear_deriv(SLP_derivative);
		delete this->SLP_derivative;
		
	} // re: clear
	
	void init();
	
	void copy(const nullspacejac_eval_data_mp & other)
	{
		
		clear();
		reset_counters();
		
		
		
		
		if (other.num_additional_linears>0) {
			
			
			this->additional_linears_terminal						= (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			this->additional_linears_terminal_full_prec = (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting						= (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting_full_prec = (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			
			for (int ii=0; ii<other.num_additional_linears; ii++) {
				vec_cp_mp(this->additional_linears_terminal[ii],						other.additional_linears_terminal[ii]);
				vec_cp_mp(this->additional_linears_terminal_full_prec[ii],	other.additional_linears_terminal_full_prec[ii]);
				
				vec_cp_mp(this->additional_linears_starting[ii],						other.additional_linears_starting[ii]);
				vec_cp_mp(this->additional_linears_starting_full_prec[ii],	other.additional_linears_starting_full_prec[ii]);
			}
		}
		else {} // other.num_additional_linears == 0
		
		this->num_additional_linears = other.num_additional_linears;
		
		if (other.num_jac_equations) {
			this->starting_linears_full_prec = (vec_mp **) br_malloc(other.num_jac_equations*sizeof(vec_mp *));
			this->starting_linears = (vec_mp **) br_malloc(other.num_jac_equations*sizeof(vec_mp *));
			
			for (int ii=0; ii<other.num_jac_equations; ii++) {
				this->starting_linears_full_prec[ii] = (vec_mp *) br_malloc(other.max_degree*sizeof(vec_mp ));
				this->starting_linears[ii] = (vec_mp *) br_malloc(other.max_degree*sizeof(vec_mp ));
				for (int jj=0; jj<other.max_degree; jj++) {
					init_vec_mp(this->starting_linears[ii][jj],0);
					init_vec_mp2(this->starting_linears_full_prec[ii][jj],0,1024);
					vec_cp_mp(this->starting_linears[ii][jj],other.starting_linears[ii][jj]);
					vec_cp_mp(this->starting_linears_full_prec[ii][jj],other.starting_linears_full_prec[ii][jj]);
				}
			}
		}
		else{}
		
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->max_degree = other.max_degree;
		
		
		if (other.num_v_linears) {
			this->v_linears = (vec_mp *) br_malloc(other.num_v_linears*sizeof(vec_mp));
			this->v_linears_full_prec = (vec_mp *) br_malloc(other.num_v_linears*sizeof(vec_mp));
			for (int ii=0; ii<num_v_linears; ii++) {
				init_vec_mp(this->v_linears[ii],0);
				init_vec_mp2(this->v_linears_full_prec[ii],0,1024);
				vec_cp_mp(this->v_linears[ii],other.v_linears[ii]);
				vec_cp_mp(this->v_linears_full_prec[ii],other.v_linears_full_prec[ii]);
			}
		}
		else{}
		
		vec_cp_mp(this->v_patch, other.v_patch);
		vec_cp_mp(this->v_patch_full_prec, other.v_patch_full_prec);
		
		mat_cp_mp(this->jac_with_proj, other.jac_with_proj);
		mat_cp_mp(this->jac_with_proj_full_prec, other.jac_with_proj_full_prec);
		

		
		
		
		
		if (other.num_projections>0) {
			this->target_projection = (vec_mp *) br_malloc(other.num_projections*sizeof(vec_mp));
			this->target_projection_full_prec = (vec_mp *) br_malloc(other.num_projections*sizeof(vec_mp));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_mp(this->target_projection[ii],0);
				init_vec_mp2(this->target_projection[ii],0,1024);
				vec_cp_mp(this->target_projection[ii],other.target_projection[ii]);
				vec_cp_mp(this->target_projection_full_prec[ii],other.target_projection_full_prec[ii]);
			}
		}
		else{}
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->target_dim = other.target_dim;   // r			the dimension of the real set we are looking for
		this->ambient_dim = other.ambient_dim;  // k			the dimension of the complex component we are looking IN.
		this->target_crit_codim = other.target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		this->num_v_vars = other.num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
		this->num_natural_vars = other.num_natural_vars;  // N-k+\ell
		this->num_synth_vars = other.num_synth_vars;
		
		this->max_degree = other.max_degree;						// the max degree of differentiated (randomized) functions
		
		
		this->num_projections = other.num_projections;
		this->num_v_linears = other.num_v_linears;
		this->num_additional_linears = other.num_additional_linears;
		
		
		
	} // re: copy
	
	
	
	
}; // re: class nullspace_eval_data_mp










// the double version
// this must be defined after the mp version, because double has mp.
class nullspacejac_eval_data_d : public solver_d
{
public:
	
	nullspacejac_eval_data_mp *BED_mp; // used only for AMP
	
	
	prog_deriv_t * SLP_derivative;
	
	
	int num_jac_equations;
	int target_dim;   // r			the dimension of the real set we are looking for
	int ambient_dim;  // k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	
	int num_natural_vars;  // N   number of variables in original problem statement (including homogenizing variables)
	int num_v_vars;  // (N-k) + (k-\ell+1)
	int num_synth_vars;
	
	int max_degree;						// the max degree of differentiated (randomized) functions
	
	int num_additional_linears;
	vec_d *additional_linears_terminal;
	
	vec_d *additional_linears_starting;
	
	
	
	vec_d **starting_linears; // outer layer should have as many as there are randomized equations
							  // inside layer has number corresponding to randomized_degrees
	
	
	int num_v_linears;
	vec_d *v_linears;         // should be as many in here as there are randomized equations
	
	vec_d v_patch;
	
	bool randomize_lower;
	mat_d lower_randomizer;
	
	mat_d jac_with_proj;
	
	
	
	int num_projections;
	vec_d *target_projection; // # of these should be target_dim (for now)
	
	
	// default initializer
	nullspacejac_eval_data_d() : solver_d(){
		
		reset_counters();
		
		init();
	}
	
	nullspacejac_eval_data_d(int mp) : solver_d(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	virtual ~nullspacejac_eval_data_d(){
		
		nullspacejac_eval_data_d::clear();
		// no need to reset the counters.
	}
	
	void print();
	
	
	void reset_counters()
	{
		num_jac_equations = 0;
		target_dim = 0;   // r			the dimension of the real set we are looking for
		ambient_dim = 0;  // k			the dimension of the complex component we are looking IN.
		target_crit_codim = 0;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		num_synth_vars = -1631;
		num_v_vars = -941;  // N   number of variables in original problem statement (including homogenizing variables)
		num_natural_vars = -746;  // N-k+\ell
		
		max_degree = -112;						// the max degree of differentiated (randomized) functions
		
		
		num_projections = 0;
		num_v_linears = 0;
		num_additional_linears = 0;
	} // re: reset_counters
	
	
	nullspacejac_eval_data_d & operator=(const nullspacejac_eval_data_d & other)
	{
		clear(); // this is wasteful, but safe for now
		reset_counters();// this is wasteful, but safe for now
		
		copy(other);
		return *this;
	}
	
	nullspacejac_eval_data_d(const nullspacejac_eval_data_d & other)
	{
		init();
		solver_d();
		nullspacejac_eval_data_d();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(prog_t * _SLP,
			  nullspace_config *ns_config,
			  witness_set & W,
			  solver_configuration & solve_options);
	
	
private:
	
	
	void clear()
	{
		
		
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_d(additional_linears_terminal[ii]);
			}
			free(additional_linears_terminal);
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_d(additional_linears_starting[ii]);
			}
			free(additional_linears_starting);
		}
		
		
		
		if (num_jac_equations>0) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					clear_vec_d(starting_linears[ii][jj]);
				}
				free(starting_linears[ii]);
			}
			free(starting_linears);
		}
		
		if (num_v_linears>0) {
			for (int ii=0; ii<num_v_linears; ii++) {
				clear_vec_d(v_linears[ii]);
			}
			free(v_linears);
		}
		
		
		clear_vec_d(v_patch);
		clear_mat_d(lower_randomizer);
		clear_mat_d(jac_with_proj);
		

		
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_d(target_projection[ii]);
			}
			free(target_projection);
		}
		
		
		
		if (this->MPType==2) {
			delete this->BED_mp;
		}
		else{
			clear_deriv(SLP_derivative);
			delete this->SLP_derivative;
		}
	} // re: clear
	
	void init();
	
	void copy(const nullspacejac_eval_data_d & other)
	{
		
		clear();
		reset_counters();
		
		
		
		
		if (other.num_additional_linears>0) {
			
			
			this->additional_linears_terminal						= (vec_d *) br_malloc(other.num_additional_linears*sizeof(vec_d));
			this->additional_linears_starting						= (vec_d *) br_malloc(other.num_additional_linears*sizeof(vec_d));
			
			for (int ii=0; ii<other.num_additional_linears; ii++) {
				vec_cp_d(this->additional_linears_terminal[ii],						other.additional_linears_terminal[ii]);
				
				vec_cp_d(this->additional_linears_starting[ii],						other.additional_linears_starting[ii]);
			}
		}
		else {} // other.num_additional_linears == 0
		
		this->num_additional_linears = other.num_additional_linears;
		
		
		
		
		if (other.num_jac_equations) {
			this->starting_linears = (vec_d **) br_malloc(other.num_jac_equations*sizeof(vec_d *));
			
			for (int ii=0; ii<other.num_jac_equations; ii++) {
				this->starting_linears[ii] = (vec_d *) br_malloc(other.max_degree*sizeof(vec_d ));
				for (int jj=0; jj<other.max_degree; jj++) {
					init_vec_d(this->starting_linears[ii][jj],0);
					vec_cp_d(this->starting_linears[ii][jj],other.starting_linears[ii][jj]);
				}
			}
		}
		else{}
		
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->max_degree = other.max_degree;
		
		
		if (other.num_v_linears) {
			this->v_linears = (vec_d *) br_malloc(other.num_v_linears*sizeof(vec_d));
			for (int ii=0; ii<num_v_linears; ii++) {
				init_vec_d(this->v_linears[ii],0);
				vec_cp_d(this->v_linears[ii],other.v_linears[ii]);
			}
		}
		else{}
		
		vec_cp_d(this->v_patch, other.v_patch);
		
		mat_cp_d(this->jac_with_proj, other.jac_with_proj);
		
		
		
		
		if (other.num_projections>0) {
			this->target_projection = (vec_d *) br_malloc(other.num_projections*sizeof(vec_d));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_d(this->target_projection[ii],0);
				vec_cp_d(this->target_projection[ii],other.target_projection[ii]);
			}
		}
		else{}
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->target_dim = other.target_dim;   // r			the dimension of the real set we are looking for
		this->ambient_dim = other.ambient_dim;  // k			the dimension of the complex component we are looking IN.
		this->target_crit_codim = other.target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		this->num_v_vars = other.num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
		this->num_natural_vars = other.num_natural_vars;  // N-k+\ell
		this->num_synth_vars = other.num_synth_vars;
		
		this->max_degree = other.max_degree;						// the max degree of differentiated (randomized) functions
		
		
		this->num_projections = other.num_projections;
		this->num_v_linears = other.num_v_linears;
		this->num_additional_linears = other.num_additional_linears;
		
		
		
	} // re: copy
	
};








/** the main function for finding critical conditions WRT a projection
 */

int nullspacejac_solver_master_entry_point(int										MPType,
										   witness_set & W, // carries with it the start points, and the linears.
										   solver_output & solve_out, // new data goes in here
										   nullspace_config				*ns_config,
										   solver_configuration		& solve_options);










int nullspacejac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED);

int nullspacejac_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);




int nullspacejac_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);



int change_nullspacejac_eval_prec(void const *ED, int new_prec);



int check_issoln_nullspacejac_d(endgame_data_t *EG,
								tracker_config_t *T,
								void const *ED);
int check_issoln_nullspacejac_mp(endgame_data_t *EG,
								 tracker_config_t *T,
								 void const *ED);

int check_isstart_nullspacejac_d(vec_d testpoint,
								 tracker_config_t *T,
								 void const *ED);

int check_isstart_nullspacejac_mp(vec_mp testpoint,
								  tracker_config_t *T,
								  void const *ED);

void check_nullspace_evaluator(point_mp current_values,
							   void const *ED);




void nullspace_slave_entry_point(solver_configuration & solve_options);







#endif


