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


#ifndef _SPHERESOLVER_H
#define _SPHERESOLVER_H

#include "missing_bertini_headers.hpp"



#include "solver.hpp"
#include "fileops.hpp"
#include "data_type.hpp"
#include "programConfiguration.hpp"
#include "postProcessing.hpp"







class sphere_config
{
	
public:
	mat_mp randomizer_matrix;  ///< R, the main randomizer matrix, which was passed in.  randomizes f and Jf down to N-k equations.
	
	SLP_global_pointers SLP_memory;
	prog_t * SLP;
	
	int MPType;
	
	bool have_rand;
	bool have_mem;
	
	vec_mp *starting_linear;
	
	vec_mp center;
	comp_mp radius;
	
	sphere_config(solver_configuration & solve_options,
				  const witness_set & W)
	{
		init();
		
		set_memory(solve_options);
		
		make_randomizer(solve_options, W);
	}
	
	
	
	sphere_config(solver_configuration & solve_options,
				  mat_mp _random)
	{
		init();
		
		set_memory(solve_options);
		
		set_randomizer(_random);
	}
	
	
	sphere_config(solver_configuration & solve_options)
	{
		init();
		set_memory(solve_options);
	}
	
	
	sphere_config(mat_mp _random)
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
	
	
	void set_memory(solver_configuration & solve_options);
	
	
	
	void set_randomizer(mat_mp _random)
	{
		mat_cp_mp(randomizer_matrix, _random);
		have_rand = true;
	}
	
	
	void set_center(vec_mp new_center)
	{
		vec_cp_mp(center, new_center);
	}
	
	void set_radius(comp_mp new_radius)
	{
		set_mp(this->radius, new_radius);
	}
	
	void copy(const sphere_config & other)
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
		
		starting_linear = (vec_mp *) br_malloc(2*sizeof(vec_mp));
		init_vec_mp2(starting_linear[0],0,1024);
		init_vec_mp2(starting_linear[1],0,1024);
		
		init_vec_mp2(center,0,1024);
		init_mp2(radius,1024);
		
		
	}
	
	
	void clear()
	{
		clear_vec_mp(center);
		clear_mp(radius);
		
		
		clear_vec_mp(starting_linear[0]);
		clear_vec_mp(starting_linear[1]);
		free(starting_linear);
		
		
		clear_mat_mp(randomizer_matrix);
		
		SLP_memory.set_globals_to_this();
		clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
		delete SLP;
	}
	
	
	
	sphere_config()
	{
		init();
	}
	
	sphere_config(const sphere_config & other)
	{
		init();
		copy(other);
	}
	
	sphere_config & operator=(const sphere_config & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	
	~sphere_config(){
		clear();
	}
};



















// the mp version
// this must be defined before the double version, because double has mp.
class sphere_eval_data_mp : public solver_mp
{
public:
	
    int num_natural_vars;
    
	comp_mp two, two_full_prec;
	
	//there had better be two starting linears.
	vec_mp *starting_linear;					// has current precision
	vec_mp *starting_linear_full_prec;			// carries full precision data for AMP
	
	
	
	vec_mp *static_linear;
	vec_mp *static_linear_full_prec;
	int num_static_linears;
	
	
	vec_mp center, center_full_prec;
	comp_mp radius, radius_full_prec;
	
	
	// default initializer
	sphere_eval_data_mp() : solver_mp(){
		init();
	}
	
	sphere_eval_data_mp(int mp) : solver_mp(mp){
		this->MPType = mp;
		init();
	}
	
	
	sphere_eval_data_mp(const sphere_eval_data_mp & other) : solver_mp(other)
	{
		init();
		copy(other);
	}
	
	sphere_eval_data_mp & operator=(const sphere_eval_data_mp & other)
	{
		
		init();
		copy(other);
		return *this;
	}
	
	
	
	virtual ~sphere_eval_data_mp(){
		sphere_eval_data_mp::clear();
	}
	
	
	
	virtual void print()
	{
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(const sphere_config & config,
			  const witness_set & W,
			  solver_configuration & solve_options);
	
	
	
protected:
	
	
	void clear()
	{
		clear_mp(two);
		clear_vec_mp(center);
		clear_mp(radius);
		
		for (int ii=0; ii<2; ii++) {
			clear_vec_mp(starting_linear[ii]);
		}
		free(starting_linear);
		
		if (num_static_linears>0) {
			for (int ii=0; ii<num_static_linears; ii++) {
				clear_vec_mp(static_linear[ii]);
			}
			
			free(static_linear);
		}
		
		
		
		
		if (this->MPType == 2) {
			clear_vec_mp(center_full_prec);
			clear_mp(radius_full_prec);
			clear_mp(two_full_prec);
			
			for (int ii=0; ii<2; ii++) {
				clear_vec_mp(starting_linear_full_prec[ii]);
			}
			free(starting_linear_full_prec);
			
			if (num_static_linears>0) {
				for (int ii=0; ii<num_static_linears; ii++) {
					clear_vec_mp(static_linear_full_prec[ii]);
				}
				
				free(static_linear_full_prec);
			}
			
			
		}
		
		
		
	} // re: clear
	
	void init();
	
	
	void copy(const sphere_eval_data_mp & other)
	{
		solver_mp::copy(other);
		
		
		if (this->num_static_linears==0) {
			static_linear = (vec_mp *) br_malloc(other.num_static_linears*sizeof(vec_mp));
		}
		else
		{
			static_linear = (vec_mp *) br_realloc(static_linear, other.num_static_linears*sizeof(vec_mp));
		}
		
		for (int ii=0; ii<other.num_static_linears; ii++) {
			init_vec_mp(static_linear[ii],0);
			vec_cp_mp(static_linear[ii],other.static_linear[ii]);
		}
		
		for (int ii=0; ii<2; ii++) {
			vec_cp_mp(starting_linear[ii], other.starting_linear[ii]);
		}
		
		set_mp(this->radius, other.radius);
		vec_cp_mp(this->center, other.center);
		
		if (this->MPType==2) {
			
			set_mp(this->radius_full_prec, other.radius_full_prec);
			vec_cp_mp(this->center_full_prec, other.center_full_prec);
			
			if (this->num_static_linears==0) {
				static_linear_full_prec = (vec_mp *) br_malloc(other.num_static_linears*sizeof(vec_mp));
			}
			else
			{
				static_linear_full_prec = (vec_mp *) br_realloc(static_linear_full_prec, other.num_static_linears*sizeof(vec_mp));
			}
			
			for (int ii=0; ii<other.num_static_linears; ii++) {
				init_vec_mp2(static_linear_full_prec[ii],0,1024);
				vec_cp_mp(static_linear_full_prec[ii],other.static_linear_full_prec[ii]);
			}
			
			for (int ii=0; ii<2; ii++) {
				vec_cp_mp(starting_linear_full_prec[ii], other.starting_linear_full_prec[ii]);
			}
		}
		
		this->num_static_linears = other.num_static_linears;
	} // re: copy
	
	
	
}; // re: class nullspace_eval_data_mp














// the mp version
// this must be defined before the double version, because double has mp.
class sphere_eval_data_d : public solver_d
{
public:
	
    int num_natural_vars;
    
	sphere_eval_data_mp * BED_mp;
	
	vec_d *starting_linear;								// has current precision
	
	comp_d two;
	
	vec_d *static_linear;
	int num_static_linears;
	
	
	vec_d center;
	comp_d radius;
	
	
	
	// default initializer
	sphere_eval_data_d() : solver_d(){
		init();
	}
	
	sphere_eval_data_d(int mp) : solver_d(mp)
	{
		this->MPType = mp;
		init();
	}
	
	
	sphere_eval_data_d(const sphere_eval_data_d & other) : solver_d(other)
	{
		init();
		copy(other);
	}
	
	sphere_eval_data_d & operator=(const sphere_eval_data_d & other)
	{
		copy(other);
		return *this;
	}
	
	
	
	virtual ~sphere_eval_data_d(){
		sphere_eval_data_d::clear();
	}
	
	
	
	virtual void print()
	{
		solver_d::print();
		
		std::cout << "sphere evaluator data (double):" << std::endl;
		for (int ii=0; ii<num_static_linears; ii++) {
			
			std::stringstream name;
			name << "static_linear_" << ii;
			print_point_to_screen_matlab(static_linear[ii],name.str());
		}
		
		for (int ii=0; ii<2; ii++) {
			
			std::stringstream name;
			name << "starting_linear_" << ii;
			print_point_to_screen_matlab(starting_linear[ii],name.str());
		}
		
		print_point_to_screen_matlab(center,"center");
		
		print_comp_matlab(radius,"radius");
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(const sphere_config & config,
			  const witness_set & W,
			  solver_configuration & solve_options);
	
	
	
protected:
	void init();
	
	void clear()
	{
		
		clear_vec_d(center);
		
		for (int ii=0; ii<2; ii++) {
			clear_vec_d(starting_linear[ii]);
		}
		free(starting_linear);
		
		if (num_static_linears>0) {
			for (int ii=0; ii<num_static_linears; ii++) {
				clear_vec_d(static_linear[ii]);
			}
			
			free(static_linear);
		}
		
		
		
		if (this->MPType==2)
		{
			delete this->BED_mp;
		}
	}
	
	
	
	void copy(const sphere_eval_data_d & other)
	{
		solver_d::copy(other);
		
		
		if (this->num_static_linears==0) {
			static_linear = (vec_d *) br_malloc(other.num_static_linears*sizeof(vec_d));
		}
		else
		{
			static_linear = (vec_d *) br_realloc(static_linear, other.num_static_linears*sizeof(vec_d));
		}
		
		for (int ii=0; ii<other.num_static_linears; ii++) {
			init_vec_d(static_linear[ii],0);
			vec_cp_d(static_linear[ii],other.static_linear[ii]);
		}
		
		for (int ii=0; ii<2; ii++) {
			vec_cp_d(starting_linear[ii], other.starting_linear[ii]);
		}
		
		set_d(this->radius, other.radius);
		vec_cp_d(this->center, other.center);
	}
	
};



int sphere_solver_master_entry_point(const witness_set						&W, // carries with it the start points, and the linears.
									 solver_output & solve_out, // new data goes in here
									 const sphere_config &		config,
									 solver_configuration		& solve_options);


int sphere_slave_entry_point(solver_configuration & solve_options);

//the new custom evaluator for this solver

int sphere_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);




int sphere_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);



int check_issoln_sphere_d(endgame_data_t *EG,
						  tracker_config_t *T,
						  void const *ED);
int check_issoln_sphere_mp(endgame_data_t *EG,
						   tracker_config_t *T,
						   void const *ED);



int sphere_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);




int change_sphere_eval_prec(void const *ED, int new_prec);










#endif


