#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

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
#include <signal.h>

#include <set>
#include <iostream>
#include <ios>
#include <string>
#include <fstream>

#include <vector>
#include <sstream>
#include <map>


#include "bertini_headers.hpp"

#include <boost/timer/timer.hpp>
#include "boost/filesystem.hpp"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/bind.hpp>

#include "fileops.hpp"




#define DEFAULT_MAX_PREC 1024


class parallelism_config; // a forward declaration
class BR_configuration;   // a forward declaration

/*** low-level data types. ***/

enum {SUCCESSFUL=0, CRITICAL_FAILURE=-10, TOLERABLE_FAILURE=-1};


enum {SYSTEM_CRIT = -1600, SYSTEM_SPHERE};

//The following lets us use words instead of numbers to indicate vertex type.
enum {UNSET= 100, CRITICAL, SEMICRITICAL, MIDPOINT, ISOLATED, NEW, CURVE_SAMPLE_POINT, SURFACE_SAMPLE_POINT, REMOVED, PROBLEMATIC};



// enum for worker mode choice
enum {NULLSPACE = 3000, LINPRODTODETJAC, DETJACTODETJAC, LINTOLIN, MULTILIN, MIDPOINT_SOLVER, SPHERE_SOLVER};

enum {TERMINATE = 2000, INITIAL_STATE};

enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};

enum {INACTIVE = 500, ACTIVE};

enum {VEC_MP = 4000, VEC_D, MAT_MP, MAT_D, COMP_MP, COMP_D, VEC_RAT, MAT_RAT, COMP_RAT, INDICES,
	DECOMPOSITION, CURVE, SURFACE, EDGE, CELL, FACE, UNUSED,
	VERTEX_SET, WITNESS_SET, VERTEX};

std::string enum_lookup(int flag);

void *br_malloc(size_t size);

void *br_realloc(void *ptr, size_t size);

void br_exit(int errorCode);

void deliberate_segfault();




template <typename key_type, typename value_type>
value_type map_lookup_with_default(const  std::map <key_type,value_type> & mc_mapperson, const key_type & lookup_key, const value_type& default_value ) {
	typename std::map<key_type,value_type>::const_iterator it = mc_mapperson.find( lookup_key );
	if ( it == mc_mapperson.end() ) {
		return default_value;
	}
	else {
		return it->second;
	}
}















bool is_identity(mat_d M);
bool is_identity(mat_mp M);

//function prototypes for bertini_real data clearing etc.

void norm_of_difference(mpf_t result, vec_mp left, vec_mp right);

void dehomogenize(point_d *result, point_d dehom_me);
void dehomogenize(point_d *result, point_d dehom_me, int num_variables);

void dehomogenize(point_mp *result, point_mp dehom_me);
void dehomogenize(point_mp *result, point_mp dehom_me, int num_variables);

void nonconj_transpose(mat_d Res, mat_d M);
void nonconj_transpose(mat_mp Res, mat_mp M);

void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);

void dot_product_mindim(comp_d result, vec_d left, vec_d right);
void dot_product_mindim(comp_mp result, vec_mp left, vec_mp right);

void projection_value_homogeneous_input(comp_d result, vec_d input, vec_d projection);
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection);


int isSamePoint_inhomogeneous_input(point_d left, point_d right);
int isSamePoint_inhomogeneous_input(point_mp left, point_mp right);



int isSamePoint_homogeneous_input(point_d left, point_d right);
int isSamePoint_homogeneous_input(point_mp left, point_mp right);


void real_threshold(comp_mp blabla, double threshold);
void real_threshold(vec_mp blabla, double threshold);
void real_threshold(mat_mp blabla, double threshold);


void print_path_retVal_message(int retVal);

/**
 retrieves the number of variables from the PPD by taking the sum of the sizes, plus the sum of the types.
 */
int get_num_vars_PPD(preproc_data PPD);


void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp PED_input);
void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d PED_input);
void cp_preproc_data(preproc_data *PPD, const preproc_data & PPD_input);

void clear_post_process_t(post_process_t * endPoint, int num_vars);


void print_tracker(const tracker_config_t * T);


int sort_increasing_by_real(vec_mp projections_sorted, std::vector< int > & index_tracker, vec_mp projections_input);

void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, std::vector< int > & randomized_degrees,
																								int num_variables, int num_funcs);
int compare_integers_decreasing(const void * left_in, const void * right_in);
int compare_integers_increasing(const void * left_in, const void * right_in);

void send_patch_mp   (patch_eval_data_mp * patch);
void receive_patch_mp(patch_eval_data_mp * patch);


void send_patch_d   (patch_eval_data_d * patch);
void receive_patch_d(patch_eval_data_d * patch);


void send_preproc_data(preproc_data *PPD);
void receive_preproc_data(preproc_data *PPD);

//
//void send_vec_mp(vec_mp b, int target);
//void receive_vec_mp(vec_mp b, int source);
//
//void send_vec_d(vec_d b, int target);
//void receive_vec_d(vec_d b, int source);



void send_mat_d(mat_d A, int target);
void receive_mat_d(mat_d A, int source);



void send_mat_mp(mat_mp A, int target);
void receive_mat_mp(mat_mp A, int source);

void send_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int target);
void receive_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int source);

void send_vec_d(vec_d b, int target);
void receive_vec_d(vec_d b, int source);


void send_vec_mp(vec_mp b, int target);
void receive_vec_mp(vec_mp b, int source);

void send_vec_rat(mpq_t ***b, int size, int target);
void receive_vec_rat(mpq_t ***b, int size, int source);


void send_comp_d(comp_d c, int target);
void receive_comp_d(comp_d c, int source);


void send_comp_num_d(comp_d *c, int num, int target);
void receive_comp_num_d(comp_d *c, int num, int source);


void send_comp_mp(comp_mp c, int target);
void receive_comp_mp(comp_mp c, int source);


void send_comp_num_mp(comp_mp *c, int num, int target);
void receive_comp_num_mp(comp_mp *c, int num, int source);

void send_comp_rat(mpq_t c[2], int target);
void receive_comp_rat(mpq_t c[2], int source);

void send_comp_num_rat(mpq_t c[][2], int num, int target);
void receive_comp_num_rat(mpq_t c[][2], int num, int source);






void print_point_to_screen_matlab(const vec_d M, std::string name);
void print_point_to_screen_matlab(const vec_mp M, std::string name);
void print_matrix_to_screen_matlab(const mat_d M, std::string name);
void print_matrix_to_screen_matlab(const mat_mp M, std::string name);

void print_comp_matlab(const comp_mp M,std::string name);
void print_comp_matlab(const comp_d M,std::string name);





class mpi_data_class
{
	
};


namespace bertini {
	class complex : public mpi_data_class
	{
		
	};
}

class function
{
	std::string func;  //symbolic representation of function (straight from input file).
										 // this class is woefully incomplete.
};



class system_randomizer
{
private:
	bool square_indicator;
	bool setup_indicator;
	
	mat_mp randomizer_matrix_full_prec;
	mat_mp randomizer_matrix_mp;
	mat_d randomizer_matrix_d;
	
	int num_randomized_funcs;
	int num_original_funcs;
	
	int max_base_degree;
	int max_degree_deficiency;
	std::vector<int> randomized_degrees;
	std::vector<int> original_degrees;
	
	std::vector<std::vector<int>> structure_matrix;
	
	
	vec_mp integer_coeffs_mp;
	mat_mp single_row_input_mp;
	vec_mp temp_homogenizer_mp;
	vec_mp temp_funcs_mp;
	mat_mp temp_jac_mp;
	mat_mp temp_mat_mp;
	vec_mp temp_vec_mp;
	
	
	vec_d integer_coeffs_d;
	mat_d single_row_input_d;
	vec_d temp_homogenizer_d;
	vec_d temp_funcs_d;
	mat_d temp_jac_d;
	mat_d temp_mat_d;
	vec_d temp_vec_d;
	
	
public:
	
	
	
	friend std::ostream & operator<<(std::ostream &os, system_randomizer & s)
	{
		os << "square: " << s.square_indicator << ", is_ready: " << s.setup_indicator << std::endl;
		os << "num_rand: " <<  s.num_randomized_funcs << ", num_orig " << s.num_original_funcs << std::endl;
		os << "max_base_deg: " << s.max_base_degree << ", max_deficiency " << s.max_degree_deficiency << std::endl;
		
		print_matrix_to_screen_matlab(s.randomizer_matrix_d,"rand_mat");
		
		for (auto iter = s.structure_matrix.begin(); iter!=s.structure_matrix.end(); ++iter) {
			for (auto jter = iter->begin(); jter!=iter->end(); ++jter) {
				os << *jter << " ";
			}
			os << std::endl;
		}
		os << std::endl << std::endl;
		
		print_point_to_screen_matlab(s.integer_coeffs_d,"int_coeffs");
		
		std::cout << "randomized" << std::endl;
		for (auto iter = s.randomized_degrees.begin(); iter!=s.randomized_degrees.end(); ++iter) {
			std::cout << *iter << " ";
		}
		os << std::endl << std::endl;
		
		std::cout << "original" << std::endl;
		for (auto iter = s.original_degrees.begin(); iter!=s.original_degrees.end(); ++iter) {
			std::cout << *iter << " ";
		}
		os << std::endl;
		
		//TODO: add mode output here
		return os;
	}
	
	
	
	system_randomizer()
	{
		init();
	}
	
	
	
	system_randomizer & operator=( const system_randomizer & other)
	{
		copy(other);
		return *this;
	}
	
	system_randomizer(const system_randomizer & other)
	{
		init();
		copy(other);
	} // re: copy
	
	~system_randomizer()
	{
		clear_mat_mp(randomizer_matrix_full_prec);
		clear_mat_mp(randomizer_matrix_mp);
		clear_mat_d(randomizer_matrix_d);
		
		clear_vec_mp(temp_homogenizer_mp);
		clear_vec_d(temp_homogenizer_d);
		
		
		
		
		clear_vec_mp(integer_coeffs_mp);
		clear_vec_d(integer_coeffs_d);
		
		
		
		clear_mat_d(single_row_input_d);
		clear_mat_mp(single_row_input_mp);
		
		
		
		clear_vec_d(temp_funcs_d);
		clear_vec_mp(temp_funcs_mp);
		
		clear_mat_d(temp_jac_d);
		clear_mat_mp(temp_jac_mp);
		
		
		clear_mat_d(temp_mat_d);
		clear_mat_mp(temp_mat_mp);
		
		clear_vec_d(temp_vec_d);
		clear_vec_mp(temp_vec_mp);
		
		
	}
	
	int num_base_funcs()
	{
		return num_original_funcs;
	}
	
	int num_rand_funcs()
	{
		return num_randomized_funcs;
	}
	
	int max_degree()
	{
		return max_base_degree;
	}
	
	int max_deficiency()
	{
		return max_degree_deficiency;
	}
	
	
	bool is_square()
	{
		return square_indicator;
	}
	
	bool is_ready()
	{
		return setup_indicator;
	}
	
	
	int base_degree(unsigned int loc)
	{
		if (loc>=original_degrees.size()) {
			br_exit(123312);
			return -1;
		}
		else{
			return original_degrees[loc];
		}
	}
	
	mat_mp * get_mat_full_prec()
	{
		return &randomizer_matrix_full_prec;
	}
	
	mat_d * get_mat_d()
	{
		return &randomizer_matrix_d;
	}
	
	mat_mp * get_mat_mp()
	{
		return &randomizer_matrix_mp;
	}
	
	
	void change_prec(int new_prec);
	
	
	
	void randomize(vec_d randomized_func_vals, mat_d randomized_jacobian,
				   vec_d func_vals, mat_d jacobian_vals,
				   comp_d hom_var);
	
	
	
	void randomize(vec_mp randomized_func_vals, mat_mp randomized_jacobian,
				   vec_mp func_vals, mat_mp jacobian_vals,
				   comp_mp hom_var);
	
	
	/**
	 parses "deg.out" for the degrees of the functions, and sets up the internals of this class object, for randomizing a system.
	 */
	void setup(int num_desired_rows, int num_funcs);
	
	
	void setup_temps();
	
	void send(int target, parallelism_config & mpi_config);
	
	void receive(int source, parallelism_config & mpi_config);
	
	void bcast_send(parallelism_config & mpi_config);
	
	void bcast_receive(parallelism_config & mpi_config);
	
	
protected:
	
	void init()
	{
		
		max_degree_deficiency = -1232;
		max_base_degree = -1321;
		
		setup_indicator = false;
		
		init_mat_mp2(randomizer_matrix_full_prec,0,0,1024);
		init_mat_mp(randomizer_matrix_mp,0,0);
		init_mat_d(randomizer_matrix_d,0,0);
		
		
		init_vec_mp(integer_coeffs_mp,0);
		init_vec_d(integer_coeffs_d,0);
		
		
		
		init_mat_d(single_row_input_d,0,0);
		init_mat_mp(single_row_input_mp,0,0);
		
		
		
		init_vec_d(temp_homogenizer_d,0);
		init_vec_mp(temp_homogenizer_mp,0);
		
		init_vec_d(temp_funcs_d,0);
		init_vec_mp(temp_funcs_mp,0);
		
		init_mat_d(temp_jac_d,0,0);
		init_mat_mp(temp_jac_mp,0,0);
		
		
		init_mat_d(temp_mat_d,0,0);
		init_mat_mp(temp_mat_mp,0,0);
		
		init_vec_d(temp_vec_d,0);
		init_vec_mp(temp_vec_mp,0);
		
		
	}
	
	void copy(const system_randomizer & other)
	{
		std::cout << "copy code for system_randomizer class is incomplete" << std::endl;
		br_exit(-1);
		
		max_degree_deficiency = other.max_degree_deficiency;
		this->square_indicator = other.square_indicator;
		
		mat_cp_d(randomizer_matrix_d, other.randomizer_matrix_d);
		if (other.randomizer_matrix_mp->rows>0 && other. randomizer_matrix_mp->cols>0) {
			change_prec_mat_mp(randomizer_matrix_mp, mpf_get_prec(other.randomizer_matrix_mp->entry[0][0].r));
		}
		
		mat_cp_mp(randomizer_matrix_mp, other.randomizer_matrix_mp);
		mat_cp_mp(randomizer_matrix_full_prec, other.randomizer_matrix_full_prec);
	}
	
	

};//re: randomizer class






class point_holder
{
	
public:
	
	vec_mp *pts_mp;
	int num_points;
	
	
	void copy_points(const point_holder & other) {
		
        for (int ii=0; ii<other.num_points; ii++)
			add_point(other.pts_mp[ii]);
    }
	
	inline bool has_no_points() const
	{
		if (num_points==0) {
			return true;
		}
		else{
			return false;
		}
	}
	
	inline bool has_points() const
	{
		if (num_points==0) {
			return false;
		}
		else{
			return true;
		}
	}
	
	void reset_points()
	{
		for (int ii =0; ii<this->num_points; ii++)
			clear_vec_mp(this->pts_mp[ii]);
		
		if (this->num_points>0) {
			free(this->pts_mp);
		}
		
		this->num_points = 0;
		this->pts_mp = NULL;
	}
	
	void add_point(vec_mp new_point);
	
	
	
	point_holder(){
		init();
	}
	
	
	~point_holder(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	point_holder& operator=( const point_holder & other) {
		
		reset_points();
		
		copy(other);
		return *this;
	}
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	point_holder(const point_holder & other){
		init();
		copy(other);
	}
	
	
	
	void copy(const point_holder & other)
	{
		copy_points(other);
	}
	
	
	
private:
	
	
	
	void init()
	{
		this->pts_mp = NULL; // potential data loss if used improperly, which i may.
		this->num_points = 0;
	}
	
	void clear(){
		reset_points();
	}
};





class patch_holder
{

public:
	
	vec_mp *patch_mp;
	int num_patches;
	
	
	void copy_patches(const patch_holder & other) {
		
        for (int ii=0; ii<other.num_patches; ii++)
			add_patch(other.patch_mp[ii]);
    }
	
	
	
	void reset_patches()
	{
		for (int ii =0; ii<this->num_patches; ii++)
			clear_vec_mp(this->patch_mp[ii]);
		
		if (this->num_patches>0) {
			free(this->patch_mp);
		}
		
		this->num_patches = 0;
		this->patch_mp = NULL;
	}
	
	void add_patch(vec_mp new_patch);
	
	
	
	patch_holder(){
		init();
	}
	
	
	~patch_holder(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	patch_holder& operator=( const patch_holder& other) {
		
		reset_patches();
		
		copy(other);
		return *this;
	}
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	patch_holder(const patch_holder & other){
		init();
		copy(other);
	}
	
	
	
	void copy(const patch_holder & other)
	{
		copy_patches(other);
	}
	
	
	
private:
	
	
	
	void init()
	{
		this->patch_mp = NULL;
		this->num_patches = 0;
	}
	
	void clear(){
		reset_patches();
	}
};



class linear_holder
{
	
public:
	
	vec_mp *L_mp;
	int num_linears;
	
	
	void copy_linears(const linear_holder & other) {
        for (int ii=0; ii<other.num_linears; ii++)
			add_linear(other.L_mp[ii]);
    }
	
	
	
	void reset_linears()
	{
		for (int ii =0; ii<num_linears; ii++)
			clear_vec_mp(L_mp[ii]);
		
		if (this->num_linears>0) {
			free(L_mp);
		}
		
		num_linears = 0;
		L_mp = NULL;
	}
	
	void add_linear(vec_mp new_linear);
	
	
	
	linear_holder(){
		init();
	}
	
	
	~linear_holder(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	linear_holder& operator=( const linear_holder& other) {
		reset_linears();
		copy(other);
		return *this;
	}
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	linear_holder(const linear_holder & other){
		init();
		copy(other);
	}
	
	
	
	void copy(const linear_holder & other)
	{
		copy_linears(other);
	}
	
	
	
private:
	
	
	
	void init()
	{
		this->L_mp = NULL;
		this->num_linears = 0;
	}
	
	void clear(){
		reset_linears();
	}
};


class name_holder
{
public:
	
	std::vector< std::string > variable_names;
	
	
	void reset_names()
	{
		variable_names.resize(0);
	}
	
	void cp_names(const name_holder & nomnom)
	{
		
		this->variable_names = nomnom.variable_names;
		
	}
	
	void get_variable_names(int num_vars) ///< reads variable names from names.out
	{
		
		variable_names.resize(num_vars);
		
		std::ifstream fin("names.out");
		
		for (int ii=0; ii<num_vars; ++ii){
			fin >> this->variable_names[ii];
		}
		fin.close();
		
	}
	
	
};



class witness_set : public patch_holder, public linear_holder, public point_holder, public name_holder
{
	
public:
	
	//begin data members
	

	int dim;
	int comp_num;
	int incidence_number;
	
	int num_variables;
	int num_synth_vars;
	


	

	
	boost::filesystem::path input_filename;
	function input_file;
	
	
	// end data members
	
	
	
	// overloaded operators
	
	
	// default constructor
	witness_set(){
		init();
	}
	
	
	~witness_set(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	witness_set& operator=( const witness_set& other)
	{
		reset();
		copy(other);
    return *this;
  }
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	witness_set(const witness_set & other)
	{
		init();
		copy(other);
	}
	
	void init()
	{
		this->input_filename = "unset_filename";
		
		this->num_variables = 0;
		this->num_synth_vars = 0;
		
		
		this->incidence_number = -1;
		this->comp_num = this->dim = -1;
		
		
		
	}
	
	void copy(const witness_set & other)
	{
		
		copy_skeleton(other);

		cp_names(other);
		copy_points(other);
		copy_patches(other);
		copy_linears(other);
	}
	
    void copy_skeleton(const witness_set & other)
	{
		this->input_filename = other.input_filename;
		
		this->dim = other.dim;
		this->comp_num = other.comp_num;
		this->incidence_number = other.incidence_number;
		
		this->num_variables = other.num_variables;
		this->num_synth_vars = other.num_synth_vars;
	}

    

    

	
    int num_natural_vars() const
    {
        return num_variables - num_synth_vars;
    }
    
    
	void only_natural_vars();
	void only_first_vars(int num_vars);
	void sort_for_real(tracker_config_t * T);
	void sort_for_unique(tracker_config_t * T);
	void sort_for_inside_sphere(comp_mp radius, vec_mp center);
	
	int  witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars);
	
	void reset()
	{
		clear();
	}
	
	
	
	void clear()
	{
		
		reset_names();
		
		reset_points();
		
		reset_linears();
		
		reset_patches();

		init();
	}
	
	

	
	

	

	
	
	void merge(const witness_set & W_in);///< merges W_in into this
	
	

	
	
	void print_to_screen() const; ///< prints some information about the witness set to the screen
	
	void print_to_file(boost::filesystem::path filename) const;
	
	/**
	 writes the linears in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void write_linears(boost::filesystem::path filename) const;
	
	/**
	 writes the patches in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void print_patches(boost::filesystem::path filename) const;
	void read_patches_from_file(boost::filesystem::path filename);
	
	
	void write_homogeneous_coordinates(boost::filesystem::path filename) const;
	void write_dehomogenized_coordinates(boost::filesystem::path filename) const;
	
    void send(parallelism_config & mpi_config, int target);
    void receive(int source, parallelism_config & mpi_config);
};


















// CURVE CELL DECOMP DATA TYPES
class vertex
{
public:
	
	point_mp pt_mp;
	
	
	vec_mp  projection_values;
	
	int type;  //See enum above.
	int removed;
	int input_filename_index;
	
	vertex()
	{
		init();
	}
	
	~vertex()
	{
		clear();
		
	}
	
	vertex & operator=(const vertex & other)
	{
		copy(other);
		return *this;
	}
	
	vertex(const vertex& other)
	{
		init();
		copy(other);
	}
	
	void print() const
	{
		print_point_to_screen_matlab(pt_mp,"point");
		print_point_to_screen_matlab(projection_values,"projection_values");
		std::cout << "type: " << type << std::endl;
	}
	
	void set_point(const vec_mp new_point);
	
	
	void send(int target, parallelism_config & mpi_config);
	
	void receive(int source, parallelism_config & mpi_config);
	
private:
	
	void clear()
	{
		clear_vec_mp(this->pt_mp);
		clear_vec_mp(this->projection_values);
	}
	
	void copy(const vertex & other)
	{
		set_point(other.pt_mp);
		
		vec_cp_mp(this->projection_values, other.projection_values);
		this->type = other.type;
		
		this->input_filename_index = other.input_filename_index;
		this->removed = other.removed;
	}
	
	
	void init()
	{
		init_point_mp2(this->projection_values,0,1024);
		init_point_mp2(this->pt_mp,1,64);
		this->pt_mp->size = 1;
		set_zero_mp(&pt_mp->coord[0]);
		this->projection_values->size = 0;
		this->type = UNSET;
		
		this->removed = 0;
		this->input_filename_index = -1;
	}
};



/**
 the main structure for storing vertices.  
 there are methods in place to add vertices, and perform lookups.
 */
class vertex_set
{
public:
	
	vec_mp *projections;
	int num_projections;
	int curr_projection;
	
	int curr_input_index;
	std::vector< boost::filesystem::path > filenames;
	
	std::vector<vertex> vertices;  //Isolated real points.
	int num_vertices;
	
	int num_natural_variables;  ///< the number of natural variables appearing in the problem to solve.
	
	
	
	mpf_t abs;
	mpf_t zerothresh;
	comp_mp diff;
	vec_mp checker_1;
	vec_mp checker_2;
	
	
	
	void print_to_screen(); ///< operator for displaying information to screen
	
	int add_vertex(const vertex new_vertex);
	int setup_vertices(boost::filesystem::path INfile);
	
	
	
	
	
	vertex_set(){
		init();
	}
	
	vertex_set(int num_vars){
		init();
		
		set_num_vars(num_vars);
	}
	

	

	
	
	vertex_set & operator=( const vertex_set& other) {
		init();
		copy(other);
		return *this;
	}
	
	vertex_set(const vertex_set &other)
	{
		init();
		copy(other);
	}
	
	~vertex_set()
	{
		vertex_set::clear();
	}
	
	
	
	void set_num_vars(int num_vars)
	{
		this->num_natural_variables = num_vars;
		
		change_size_vec_mp(checker_1, num_vars);
		change_size_vec_mp(checker_2, num_vars);
		checker_1->size = checker_2->size = num_vars;
	}
	
	
	
	
	void print(boost::filesystem::path outputfile) const;
	
	
	
	
	int set_curr_input(boost::filesystem::path el_nom){
		
		int nom_index = -1;
        
		if ( el_nom.string().compare("unset_filename")==0 )
		{
			std::cout << "trying to set curr_input from unset_filename" << std::endl;
			
			br_exit(-59251);
		}
		
		
		int counter = 0;
		for (std::vector<boost::filesystem::path>::iterator ii = filenames.begin(); ii!= filenames.end(); ++ii) {
			
			if (*ii == el_nom) {
				nom_index = counter;
				break;
			}
			counter++;
		}
		
		if (nom_index==-1) {
			filenames.push_back(el_nom);
			nom_index = counter;
		}
		
		curr_input_index = nom_index;
		return nom_index;
	}
	
	
    

    
    int search_for_point(vec_mp testpoint);
    int search_for_active_point(vec_mp testpoint);
    int search_for_removed_point(vec_mp testpoint);
    
	/**
     \param W witness set containing points of which we wish to retrieve projections values.
     \return  crit_downstairs the projection values of the input witness set, sorted for uniqueness and increasingness.
     \return midpoints_downstairs the bisection of each interval in crit_downstairs.
     \return index_tracker the indices of the points in W.
     \param pi the projection we are retrieving projection values with respect to.
     
     */
    int compute_downstairs_crit_midpts(const witness_set & W,
                                       vec_mp crit_downstairs,
                                       vec_mp midpoints_downstairs,
                                       std::vector< int > & index_tracker,
                                       vec_mp pi);
    
    
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value);
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index);
	
	int set_curr_projection(vec_mp new_proj){
        
        int proj_index = get_proj_index(new_proj);
        
        if (proj_index==-1) {
            int init_size = new_proj->size;
            new_proj->size = this->num_natural_variables;
            
			proj_index = add_projection(new_proj);
            
            new_proj->size = init_size;
		}
		
        curr_projection = proj_index;
        
//        std::cout << "curr_projection is now: " << curr_projection << std::endl;
		return proj_index;
	}
	
    
    int get_proj_index(vec_mp proj) const
    {
        int init_size = proj->size;
        proj->size = this->num_natural_variables;
        
        
        int proj_index = -1;
        
        for (int ii=0; ii<num_projections; ii++) {
			if (isSamePoint_inhomogeneous_input(projections[ii],proj)) {
				proj_index = ii;
				break;
			}
		}
        
        proj->size = init_size;
        
        return proj_index;
    }
    
    
	
	int add_projection(vec_mp proj){
		
		if (this->num_projections==0) {
			projections = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->projections = (vec_mp *)br_realloc(this->projections, (this->num_projections+1) * sizeof(vec_mp));
		}
		
		
		init_vec_mp2(this->projections[num_projections],num_natural_variables,DEFAULT_MAX_PREC);
		this->projections[num_projections]->size = num_natural_variables;
		
		
		if (proj->size != num_natural_variables) {
			vec_mp tempvec;
			init_vec_mp2(tempvec,num_natural_variables,DEFAULT_MAX_PREC); tempvec->size = num_natural_variables;
			for (int kk=0; kk<num_natural_variables; kk++) {
				set_mp(&tempvec->coord[kk], &proj->coord[kk]);
			}
			vec_cp_mp(projections[num_projections], tempvec);
			clear_vec_mp(tempvec);
		}
		else
		{
			vec_cp_mp(projections[num_projections], proj);
		}
		

		num_projections++;
		
		return num_projections;
	}
	
	
	
	
	void send(int target, parallelism_config & mpi_config);
	
	
	
	void receive(int source, parallelism_config & mpi_config);
	
	void reset()
	{
		for (int ii=0; ii<num_projections; ii++) {
			clear_vec_mp(projections[ii]);
		}
		num_projections = 0;
		
		filenames.resize(0);
		
		vertices.resize(0);
		num_vertices = 0;
		
		clear();
		init();
	}
protected:
	
	void init()
	{
		
		num_projections = 0;
		projections = NULL;
		curr_projection = -1;
		
		curr_input_index = -2;
		
		this->num_vertices = 0;
		this->num_natural_variables = 0;
		
		init_vec_mp(checker_1,0);
		init_vec_mp(checker_2,0);
		
		
		
		init_mp(this->diff);

		mpf_init(abs);
		mpf_init(zerothresh);
		mpf_set_d(zerothresh, 1e-8);
	}
	
	
	void copy(const vertex_set &other)
	{
		this->curr_projection = other.curr_projection;
		if (this->num_projections==0 && other.num_projections>0) {
			projections = (vec_mp *) br_malloc(other.num_projections*sizeof(vec_mp));
		}
		else if(this->num_projections>0 && other.num_projections>0) {
			projections = (vec_mp *) br_realloc(projections,other.num_projections*sizeof(vec_mp));
		}
		else if (this->num_projections>0 && other.num_projections==0){
			for (int ii=0; ii<this->num_projections; ii++) {
				clear_vec_mp(projections[ii]);
			}
			free(projections);
		}
		
		for (int ii=0; ii<other.num_projections; ii++) {
			if (ii>=this->num_projections){
				init_vec_mp2(projections[ii],1,1024); projections[ii]->size = 1;
			}
			vec_cp_mp(projections[ii],other.projections[ii]);
		}
		
		this->num_projections = other.num_projections;
		
		
		curr_input_index = other.curr_input_index;
		filenames = other.filenames;
		
		
		this->num_vertices = other.num_vertices;
		this->num_natural_variables = other.num_natural_variables;
		
		this->vertices = other.vertices;
		
		vec_cp_mp(this->checker_1,other.checker_1);
		vec_cp_mp(this->checker_2,other.checker_2);
		
	}
	
	void clear()
	{
		clear_vec_mp(checker_1);
		clear_vec_mp(checker_2);
		
		clear_mp(diff);
		
		mpf_clear(abs);
		mpf_clear(zerothresh);
		
		for (int ii=0; ii<num_projections; ii++) {
			clear_vec_mp(projections[ii]);
		}
		free(projections);
	}

};




class cell
{
	
private:
//	int n;
//	function homotopy;
	
public:
	int midpt; ///< index into vertex set
	
	
	friend std::istream & operator>>(std::istream &os, cell & c)
	{
		os >> c.midpt;
		return os;
	}
	
	void copy(const cell & other){
		this->midpt = other.midpt;
	}
	
	void send(int target, parallelism_config & mpi_config)
	{
		int buffer = midpt;
		MPI_Send(&buffer, 1, MPI_INT, target, CELL, MPI_COMM_WORLD);
	}
	
	void receive(int source, parallelism_config & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int buffer;
		MPI_Recv(&buffer, 1, MPI_INT, source, CELL, MPI_COMM_WORLD, &statty_mc_gatty);
		midpt = buffer;
	}
	
	virtual void read_from_stream( std::istream &is ) = 0;
	
};



/**
 the edge data type.  has three indices: left, right, midpt.
 */
class edge : public cell
{
public:
	int left;  ///< index into vertices
	int right; ///< index into vertices
	int midpt; ///<  index into vertices
	
	std::vector< int > removed_points;
	
	
	
	edge() : cell()
	{
		left = right = midpt = -1;
	}
	
	edge(int left_, int midpt_, int right_)
	{
		this->left = left_;
		this->right = right_;
		this->midpt = midpt_;
	}
	
	// other defaults are correct for this type.
	
	inline bool is_degenerate()
	{
		if ((left == right) || (left==midpt) || (right==midpt))
			return true;
		else
			return false;
	}
	
	void send(int target, parallelism_config & mpi_config)
	{
		int * buffer = new int[4];
		
		buffer[0] = left;
		buffer[1] = midpt;
		buffer[2] = right;
		buffer[3] = removed_points.size();
		
		MPI_Send(buffer, 4, MPI_INT, target, EDGE, MPI_COMM_WORLD);
		
		delete [] buffer;
		
		buffer = new int[removed_points.size()];
		for (unsigned int ii=0; ii!=removed_points.size(); ii++) {
			buffer[ii] = removed_points[ii];
		}
		MPI_Send(buffer, removed_points.size(), MPI_INT, target, EDGE, MPI_COMM_WORLD);
		delete [] buffer;
		

	}
	
	void receive(int source, parallelism_config & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int * buffer = new int[4];
		MPI_Recv(buffer, 4, MPI_INT, source, EDGE, MPI_COMM_WORLD, &statty_mc_gatty);
		
		left  = buffer[0];
		midpt = buffer[1];
		right = buffer[2];
		
		int temp_num_removed = buffer[3];
		
		
		delete [] buffer;
		
		buffer = new int[temp_num_removed];
		MPI_Recv(buffer, temp_num_removed, MPI_INT, source, EDGE, MPI_COMM_WORLD, &statty_mc_gatty);
		for (int ii=0; ii<temp_num_removed; ii++) {
			removed_points.push_back(buffer[ii]);
		}
		
		delete [] buffer;
		
	}
	
	
	
	virtual void read_from_stream( std::istream &os )
	{
		
	}
	
	
	
};












class decomposition : public patch_holder
{

public:
	std::map< int , int > counters;
	std::map< int , std::vector< int > > indices;
	
	witness_set W;
	
	int num_variables;
	int dimension;
	int component_num;
	
	int num_curr_projections;
	vec_mp	*pi; // the projections
	
	system_randomizer * randomizer;

	vec_mp sphere_center;
	comp_mp sphere_radius;
	bool have_sphere_radius;
	
	boost::filesystem::path input_filename;
//	function input_file;
	
	
	
    

	void add_projection(vec_mp proj){
		if (this->num_curr_projections==0) {
			pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->pi = (vec_mp *)br_realloc(this->pi, (this->num_curr_projections+1) * sizeof(vec_mp));
		}
		
		init_vec_mp2(this->pi[num_curr_projections],proj->size,DEFAULT_MAX_PREC);
		this->pi[num_curr_projections]->size = proj->size;
		
		vec_cp_mp(pi[num_curr_projections], proj);
		num_curr_projections++;
		
	}
	
	
	
    int add_witness_set(const witness_set & W, int add_type, vertex_set & V);
    
	int index_in_vertices(vertex_set &V,
                          vec_mp testpoint);
	
	int index_in_vertices_with_add(vertex_set &V,
                                   vertex vert);
	
	int setup(boost::filesystem::path INfile);
	
	virtual void print(boost::filesystem::path outputfile);
	
	
	
	
	int read_sphere(const boost::filesystem::path & bounding_sphere_filename);
	
	void compute_sphere_bounds(const witness_set & W_crit);
	
	void copy_sphere_bounds(const decomposition & other)
	{
		set_mp(this->sphere_radius, other.sphere_radius);
		vec_cp_mp(this->sphere_center, other.sphere_center);
		this->have_sphere_radius = true;
	}
	
	
	
	void output_main(const boost::filesystem::path base);
	
	
	
	
	
	
	
	void reset()
	{
		clear();
		init();
	}
	
	
	
	
	decomposition(){
		init();
	}
	
	virtual ~decomposition()
	{
		this->clear();
	}
	
	decomposition & operator=(const decomposition& other){
		
		this->init();
		
		this->copy(other);
		
		return *this;
	}
	
	decomposition(const decomposition & other){
		
		this->init();
		
		this->copy(other);
	}
	
	
	
	void send(int target, parallelism_config & mpi_config);
	void receive(int source, parallelism_config & mpi_config);
	
protected:
	
	int add_vertex(vertex_set &V, vertex source_vertex);
	
	
	
	void init(){
		
		randomizer = new system_randomizer;

		input_filename = "unset";
		pi = NULL;
		
		
		num_curr_projections = 0;
		num_variables = 0;
		dimension = -1;
		component_num = -1;
		
		init_mp2(sphere_radius,DEFAULT_MAX_PREC);
		init_vec_mp2(sphere_center,0,1024);
		sphere_center->size = 0;
		have_sphere_radius = false;
		
		set_one_mp(sphere_radius);
		neg_mp(sphere_radius,sphere_radius);
		

		


		
		
		
	}
	
	
	void copy(const decomposition & other)
	{
		
		
		patch_holder::copy(other);
		
		
		*this->randomizer = *other.randomizer;
		
		this->W = other.W;
		
		this->input_filename = other.input_filename;
		
		this->counters = other.counters;
		this->indices = other.indices;
		this->num_variables = other.num_variables;
		this->dimension = other.dimension;
		this->component_num = other.component_num;
		
		this->num_curr_projections = other.num_curr_projections;
		this->pi = (vec_mp *) br_malloc(other.num_curr_projections * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_curr_projections; ii++) {
			init_vec_mp2(this->pi[ii],other.pi[ii]->size,DEFAULT_MAX_PREC);
			this->pi[ii]->size = other.pi[ii]->size;
			vec_cp_mp(this->pi[ii], other.pi[ii])
		}
		
		
		
		
		
		copy_sphere_bounds(other);
		
		return;
	}
	
	void clear()
	{
		
		delete randomizer;
		
		if (num_curr_projections>0){
			for (int ii=0; ii<num_curr_projections; ii++)
				clear_vec_mp(pi[ii]);
			free(pi);
		}
		num_curr_projections = 0;
		
		
		counters.clear();
		indices.clear();
		
		
		clear_mp(sphere_radius);
		clear_vec_mp(sphere_center);
		

		
	}
	
	
	
	

	
}; // end decomposition














namespace color {
	
	int color_to_int(const std::string c);
	
	
	std::string bold(std::string new_color);
	
	std::string dark(std::string new_color);
	
	
	std::string underline(std::string new_color);
	
	
	std::string background(std::string new_color);
	
	
	std::string strike(std::string new_color);
			 
	std::string console_default();
	
	std::string black();
	
	std::string red();
	
	std::string green();

	
	std::string brown();
	
	std::string blue();
	
	std::string magenta();
	
	std::string cyan();
	
	std::string gray();
	
	
	//black - 30
	//red - 31
	//green - 32
	//brown - 33
	//blue - 34
	//magenta - 35
	//cyan - 36
	//lightgray - 37
	
}


#endif

