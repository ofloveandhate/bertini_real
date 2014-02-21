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

class midpoint_config;
class midpoint_eval_data_mp;
class midpoint_eval_data_d;

class complete_system
{
	friend class midpoint_config;
	friend class midpoint_eval_data_mp;
	friend class midpoint_eval_data_d;
protected:
	
	boost::filesystem::path input_filename;
	SLP_global_pointers memory;
	prog_t * SLP;
	mat_mp randomizer_matrix;
	int num_variables;
	
	bool have_SLP;
	
	int MPType;
	
public:
	
	void bcast_send(parallelism_config & mpi_config)
	{
		std::cout << "bcasting " << input_filename << std::endl;
		
		
		int * buffer = new int[4];
		buffer[0] = randomizer_matrix->rows;
		buffer[1] = randomizer_matrix->cols;
		buffer[2] = MPType;
		buffer[3] = num_variables;
		MPI_Bcast(buffer, 4, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
		
		delete [] buffer;
		
		bcast_mat_mp(randomizer_matrix,0,0);
		memory.set_globals_to_this();
		bcast_prog_t(SLP, this->MPType, 0, 0); // last two arguments are: myid, headnode
		memory.set_globals_null();
	}
	
	
	
	void bcast_receive(parallelism_config & mpi_config)
	{
		
		if (have_SLP) {
			clear_mat_mp(randomizer_matrix);
			memory.set_globals_to_this();
			clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
		}
		else
		{
			init_mat_mp2(randomizer_matrix,1,1,1024);
			randomizer_matrix->rows = 1;
			randomizer_matrix->cols = 1;
			SLP = new prog_t;
		}
		
		int * buffer = new int[4];
		
		MPI_Bcast(buffer, 4, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
		
		change_size_mat_mp(randomizer_matrix,buffer[0],buffer[1]);

		randomizer_matrix->rows = buffer[0];
		randomizer_matrix->cols = buffer[1];
		MPType = buffer[2];
		num_variables = buffer[3];
		
		
		bcast_mat_mp(randomizer_matrix,1,0);
		
		
		delete [] buffer;
		
		bcast_prog_t(SLP, MPType, 1, 0); // last two arguments are: myid, headnode
		initEvalProg(MPType);
		memory.capture_globals();
		memory.set_globals_null();
		
		have_SLP = true;
	}
	
	
	
	complete_system()
	{
		init();
	}
	
	complete_system(const decomposition & D, tracker_config_t * T)
	{
		init();
		get_system(D,T);
	}
	
	
	~complete_system()
	{
		clear();
	}
	
	
	complete_system & operator=( const complete_system & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	complete_system(const complete_system & other)
	{
		init();
		copy(other);
	} // re: copy
	
	
	
	void init()
	{
		num_variables = 0;
		input_filename = "unset_input_filename_complete_class";
		MPType = -1;
		
		have_SLP = false;
	}
	
	
	void copy(const complete_system & other)
	{

		this->num_variables = other.num_variables;
		
		this->MPType = other.MPType;
		
		this->input_filename = other.input_filename;
		
		if (this->have_SLP) {
			this->memory.set_globals_to_this();
			clearProg(this->SLP, this->MPType, 1); // 1 means call freeprogeval()
			
			if (!other.have_SLP) {
				delete this->SLP;
				clear_mat_mp(this->randomizer_matrix);
				this->have_SLP = false;
			}
			
		}
		
		if (other.have_SLP) {
			
			if (!this->have_SLP)
			{
				this->SLP = new prog_t;
				init_mat_mp2(this->randomizer_matrix,1,1,1024);
				this->randomizer_matrix->rows = 1; this->randomizer_matrix->cols = 1;
			}
			
			mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
			
			cp_prog_t(this->SLP, other.SLP);
			this->memory.set_globals_null();
			initEvalProg(this->MPType);
			this->memory.capture_globals();
			
			this->have_SLP = true;
		}
		
		
	}
	
	
	void clear()
	{
		std::cout << "clearing complete_system, name " << input_filename << std::endl;
		
		if (have_SLP) {
			clear_mat_mp(randomizer_matrix);
			memory.set_globals_to_this();
			clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
			delete SLP;
			have_SLP = false;
		}
		
		
	}
	
	void get_system(const decomposition & D, tracker_config_t * T)
	{
		if (have_SLP) {
			memory.set_globals_to_this();
			clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
		}
		else{
			SLP = new prog_t;
			init_mat_mp2(randomizer_matrix,1,1,1024);
			randomizer_matrix->rows = 1; randomizer_matrix->cols = 1;
		}
		
		
		
		
		int blabla;  // i would like to move this.
		parse_input_file(D.input_filename, &blabla);
		num_variables = D.num_variables;
		input_filename = D.input_filename;
		this->MPType = T->MPType;
		
		
		
		
		//	// setup a straight-line program, using the file(s) created by the parser
		int numVars = setupProg(SLP, T->Precision, T->MPType);
		if (num_variables != numVars) {
			std::cout << "numvars is incorrect..." << std::endl;
			mypause();
			br_exit(-57189); // this should be a throw or something
		}
		
		preproc_data PPD;
		parse_preproc_data("preproc_data", &PPD);
		
		//create the array of integers
		std::vector<int> randomized_degrees;
		
		//get the matrix and the degrees of the resulting randomized functions.
		make_randomization_matrix_based_on_degrees(randomizer_matrix, randomized_degrees, D.num_variables-D.num_patches-D.dimension, PPD.num_funcs);
		
		
		
		memory.capture_globals();
		memory.set_globals_null();
		
		
		
		preproc_data_clear(&PPD);
		
		have_SLP = true;
		std::cout << "have system " << input_filename << " in memory" << std::endl;
	}
};

//
class surface_decomposition; // forward declaration
//
class midpoint_config
{
public:
	friend class complete_system;
	
	std::map<std::string, complete_system> systems;
	
	
    
    int MPType;
    
    
	
	comp_mp v_target; // set during the loop in connect the dots
	comp_mp u_target;// set during the loop in connect the dots
	
	comp_mp crit_val_left;// set during the loop in connect the dots
	comp_mp crit_val_right;// set during the loop in connect the dots
	
	std::string system_name_mid;
	std::string system_name_bottom;
	std::string system_name_top;
	
	
	vec_mp *pi;
	int num_projections;
    
    
	
	midpoint_config(){
		init();
	}
	
	~midpoint_config(){
		clear();
	}
	
	
	midpoint_config & operator=(const midpoint_config & other)
	{
		init();
		
		copy(other);
		return *this;
	}
	
	midpoint_config(const midpoint_config & other)
	{
		init();
		
		copy(other);
		
	}
	
	

	
	
    void print()
	{

		print_comp_matlab(u_target,"u_target");
		print_comp_matlab(v_target,"v_target");
		
		print_comp_matlab(crit_val_right,"crit_val_right");
		print_comp_matlab(crit_val_left,"crit_val_left");
		

	}
	
	void setup(const surface_decomposition & surf,
               solver_configuration & solve_options);
    
    void bcast_send(parallelism_config & mpi_config);
    void bcast_receive(parallelism_config & mpi_config);

	
	
    void add_projection(vec_mp proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else {
			this->pi = (vec_mp *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_mp));
		}
		
		init_vec_mp2(pi[num_projections],0,1024);
		vec_cp_mp(pi[num_projections], proj);
		
		num_projections++;
		
	}
	
	
	int num_top_vars()
	{
		return systems[system_name_top].num_variables;
	}
	
	
	
	int num_bottom_vars()
	{
		return systems[system_name_bottom].num_variables;
	}
	
	
	int num_mid_vars()
	{
		return systems[system_name_mid].num_variables;
	}
	
	
	
	
	
	
private:
    
	void copy(const midpoint_config & other){
		
		this->MPType = other.MPType;
		
		this->systems = other.systems;
		

		
		system_name_mid = other.system_name_mid;
		system_name_top = other.system_name_top;
		system_name_bottom = other.system_name_bottom;
		

		
		set_mp(v_target,other.v_target);
		set_mp(u_target,other.u_target);
		set_mp(crit_val_left,other.crit_val_left);
		set_mp(crit_val_right,other.crit_val_right);
		

        
        for (int ii=0; ii<other.num_projections; ii++) {
            add_projection(other.pi[ii]);
        }
        
	}
	
	
	
	
	void init();
	
	void clear()
	{
		
		clear_mp(v_target);
		clear_mp(u_target);
		clear_mp(crit_val_left);
		clear_mp(crit_val_right);
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_mp(pi[ii]);
			}
			free(pi);
		}
		
	}
	
    
	
    
	
};






// the mp version
// this must be defined before the double version, because double has mp.
class midpoint_eval_data_mp : public solver_mp
{
public:
	
	int num_mid_vars;
	int num_bottom_vars;
	int num_top_vars;
	
    
	SLP_global_pointers top_memory;
	SLP_global_pointers mid_memory;
	SLP_global_pointers bottom_memory;
	
	// these are all merely pointers, and should only be assigned to memory set by the SLP creation routine inside of setupProg(), by the midpoint_config::setup() call
	
	prog_t *SLP_top;
	prog_t *SLP_bottom;
	prog_t *SLP_mid; // these are to be freed by the midpoint_setup object.
	
	
	mat_mp randomizer_matrix_bottom;
	mat_mp randomizer_matrix_top;
	
	
	
	vec_mp *pi;
	int num_projections;
	
	//patch already lives in the base class.
	
	comp_mp v_target;
	comp_mp u_target;
	
    
	
	comp_mp crit_val_left;
	comp_mp crit_val_right;
	
	
	comp_mp u_start;
	comp_mp v_start;
	
	comp_mp one;
	comp_mp zero;
	comp_mp half;
	
	
	vec_mp *pi_full_prec;
	comp_mp v_target_full_prec;
	comp_mp u_target_full_prec;
	
	mat_mp randomizer_matrix_bottom_full_prec;
	mat_mp randomizer_matrix_top_full_prec;
	
	comp_mp crit_val_left_full_prec;
	comp_mp crit_val_right_full_prec;
	
	comp_mp half_full_prec;
	comp_mp u_start_full_prec;
	comp_mp v_start_full_prec;
	comp_mp one_full_prec;
	comp_mp zero_full_prec;
	
	
	
	
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
		
		clear();
		// no need to reset the counters.
	}
	
	void print()
	{
        
	};
	
	void reset_counters()
	{
		this->num_projections = 0;
	} // re: reset_counters
	
	
	midpoint_eval_data_mp & operator=(const midpoint_eval_data_mp & other)
	{
		init();
		
		copy(other);
		return *this;
	}
	
	midpoint_eval_data_mp(const midpoint_eval_data_mp & other)
	{
		init();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	int send(parallelism_config & mpi_config);
	
	int receive(parallelism_config & mpi_config);
	
	
	
	int setup(midpoint_config & md_config,
              const witness_set & W,
              solver_configuration & solve_options);
	
	
	
	
	void add_projection(vec_mp proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else {
			this->pi = (vec_mp *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_mp));
		}
		
		init_vec_mp(pi[num_projections],0);
		vec_cp_mp(pi[num_projections], proj);
		
		if (this->MPType==2){
			if (this->num_projections==0) {
				this->pi_full_prec = (vec_mp *) br_malloc(sizeof(vec_mp));
			}
			else {
				this->pi_full_prec = (vec_mp *) br_realloc(pi_full_prec,(this->num_projections+1) *sizeof(vec_mp));
			}
			
			init_vec_mp2(pi_full_prec[num_projections],0,1024);
			vec_cp_mp(pi_full_prec[num_projections], proj);
		}
		num_projections++;
		
	}
	
	
	
	
	
	void clear()
	{
        
		if (this->num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_mp(pi[ii]);
			}
			free(pi);
		}
		
		clear_mp(v_target);
		clear_mp(u_target);
		clear_mat_mp(randomizer_matrix_top);
		clear_mat_mp(randomizer_matrix_bottom);
		clear_mp(crit_val_left);
		clear_mp(crit_val_right);
		
		clear_mp(u_start);
		clear_mp(v_start);
		
		clear_mp(half);
		clear_mp(one);
		clear_mp(zero);
		
		
		if (this->MPType==2){
            
			if (this->num_projections>0) {
				for (int ii=0; ii<num_projections; ii++) {
					clear_vec_mp(pi_full_prec[ii]);
				}
				free(pi_full_prec);
			}
			
			clear_mp(v_target_full_prec);
			clear_mp(u_target_full_prec);
			clear_mat_mp(randomizer_matrix_top_full_prec);
			clear_mat_mp(randomizer_matrix_bottom_full_prec);
			clear_mp(crit_val_left_full_prec);
			clear_mp(crit_val_right_full_prec);
			
			clear_mp(u_start_full_prec);
			clear_mp(v_start_full_prec);
			
			clear_mp(half_full_prec);
			clear_mp(one_full_prec);
			clear_mp(zero_full_prec);
			
		}
        
		
		
	} // re: clear
	
	void init(); // in the .cpp file.  must be there, not here.
    
	
	void copy(const midpoint_eval_data_mp & other)
	{
		solver_mp::copy(other);
		
		this->num_mid_vars = other.num_mid_vars;
		this->num_top_vars = other.num_top_vars;
		this->num_bottom_vars = other.num_bottom_vars;
		
		this->bottom_memory = other.bottom_memory;
		this->top_memory = other.top_memory;
		this->mid_memory = other.mid_memory;
		
		this->SLP_mid = other.SLP_mid;
		this->SLP_top = other.SLP_top;
		this->SLP_bottom = other.SLP_bottom;
		
		for (int ii=0; ii<other.num_projections; ii++) {
			if (other.MPType==2) {
				add_projection(other.pi_full_prec[ii]);
			}
			else{
				add_projection(other.pi[ii]);
			}
		}
		
		//patch already lives in the base class.
		
		set_mp(v_target,other.v_target);
		set_mp(u_target,other.u_target);
		
		mat_cp_mp(randomizer_matrix_top, other.randomizer_matrix_top);
		mat_cp_mp(randomizer_matrix_bottom, other.randomizer_matrix_bottom);
		
		
		set_mp(crit_val_left,other.crit_val_left);
		set_mp(crit_val_right,other.crit_val_right);
		
		set_mp(u_start, other.u_start);
		set_mp(v_start, other.v_start);
		
		
		if (this->MPType==2) {
			set_mp(v_target_full_prec,other.v_target_full_prec);
			set_mp(u_target_full_prec,other.u_target_full_prec);
			
			mat_cp_mp(randomizer_matrix_top_full_prec,other.randomizer_matrix_top_full_prec);
			mat_cp_mp(randomizer_matrix_bottom_full_prec,other.randomizer_matrix_bottom_full_prec);
			
			
			set_mp(crit_val_left_full_prec,other.crit_val_left_full_prec);
			set_mp(crit_val_right_full_prec,other.crit_val_right_full_prec);
			
			set_mp(u_start_full_prec, other.u_start_full_prec);
			set_mp(v_start_full_prec, other.v_start_full_prec);
		}
        
	} // re: copy
	
	
	
	
}; // re: class nullspace_eval_data_mp










// the double version
// this must be defined after the mp version, because double has mp.
class midpoint_eval_data_d : public solver_d
{
public:
	
	midpoint_eval_data_mp *BED_mp; // used only for AMP
	
	int num_mid_vars;
	int num_top_vars;
	int num_bottom_vars;
	
	SLP_global_pointers bottom_memory;
	SLP_global_pointers top_memory;
	SLP_global_pointers mid_memory;
	
	prog_t *SLP_bottom;
	prog_t *SLP_top;
	prog_t *SLP_mid;
	
	vec_d *pi;
	int num_projections;
	
	mat_d randomizer_matrix_top;
	mat_d randomizer_matrix_bottom;
	
	
	//patch already lives in the base class.
	
	comp_d v_target;
	comp_d u_target;
	
	
	
	comp_d crit_val_left;
	comp_d crit_val_right;
	
	comp_d half;
	comp_d u_start;
	comp_d v_start;
	comp_d one;
	comp_d zero;
	
	
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
		std::cout << "num_variables " << num_variables << std::endl;
		std::cout << "num_mid_vars " << num_mid_vars << std::endl;
		std::cout << "SLP_top "      << SLP_top->size << std::endl;
		std::cout << "SLP_bottom " << SLP_bottom->size << std::endl;
		std::cout << "SLP_mid " << SLP_mid->size << std::endl;
        
		print_comp_matlab(u_target,"u_target");
		print_comp_matlab(v_target,"v_target");
		
		print_comp_matlab(crit_val_left,"crit_val_left");
		print_comp_matlab(crit_val_right,"crit_val_right");
		
		
//		std::cout << num_projections << " projections: " << std::endl;
//		for (int ii=0; ii<2; ii++) {
//			print_point_to_screen_matlab(this->pi[ii],"pi");
//		}
		
	};
	
	
	
	void reset_counters()
	{
		
	} // re: reset_counters
	
	
	midpoint_eval_data_d & operator=(const midpoint_eval_data_d & other)
	{
		this->MPType = other.MPType;
		init();
		
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
	
	
	
	int setup(midpoint_config & md_config,
              const witness_set & W,
              solver_configuration & solve_options);
	
	
	void add_projection(vec_mp proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_d *) br_malloc(sizeof(vec_d));
		}
		else {
			this->pi = (vec_d *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_d));
		}
		
		init_vec_d(pi[num_projections],0);
		vec_mp_to_d(pi[num_projections], proj);
		
		num_projections++;
	}
	
	
	void add_projection(vec_d proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_d *) br_malloc(sizeof(vec_d));
		}
		else {
			this->pi = (vec_d *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_d));
		}
		
		init_vec_d(pi[num_projections],0);
		vec_cp_d(pi[num_projections], proj);
		
		num_projections++;
	}
private:
	
	
	void clear()
	{
        
		//yeah, there's some stuff to clear here.
		if (this->MPType==2) {
			delete this->BED_mp;
		}
		
		
		for (int ii=0; ii<num_projections; ii++) {
			clear_vec_d(pi[ii]);
		}
		free(pi);
		
		
		clear_mat_d(randomizer_matrix_top);
		clear_mat_d(randomizer_matrix_bottom);
	} // re: clear
	
	void init();
	
	void copy(const midpoint_eval_data_d & other)
	{
		
		
		solver_d::copy(other);
		
		this->num_mid_vars = other.num_mid_vars;
		this->num_top_vars = other.num_top_vars;
		this->num_bottom_vars = other.num_bottom_vars;
		
		this->top_memory = other.top_memory;
		this->bottom_memory = other.bottom_memory;
		this->mid_memory = other.mid_memory;
		
		this->SLP_bottom = other.SLP_bottom;
		this->SLP_mid = other.SLP_mid;
		this->SLP_top = other.SLP_top;
        
		
		for (int ii=0; ii<other.num_projections; ii++) {
			add_projection(other.pi[ii]);
		}
		
		//patch already lives in the base class.
		
		set_d(v_target,other.v_target);
		set_d(u_target,other.u_target);
		
		mat_cp_d(randomizer_matrix_top,other.randomizer_matrix_top);
		mat_cp_d(randomizer_matrix_bottom,other.randomizer_matrix_bottom);
		
		
		set_d(crit_val_left,other.crit_val_left);
		set_d(crit_val_right,other.crit_val_right);
		
		set_d(u_start, other.u_start);
		set_d(v_start, other.v_start);
	}
};







/** the main function for finding critical conditions WRT a projection
 */

int midpoint_solver_master_entry_point(const witness_set & W, // carries with it the start points, and the linears.
                                       solver_output & solve_out, // new data goes in here
                                       midpoint_config & md_config,
                                       solver_configuration		& solve_options);






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


