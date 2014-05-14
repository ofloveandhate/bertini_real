#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H


/** \file data_type.hpp */

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

#include "bertini_extensions.hpp"
#include "fileops.hpp"
#include "checkSelfConjugate.hpp"




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







std::string enum_lookup(int flag);








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

















/** 
 \brief a woefully incomplete class to contain systems which bertini will parse.
 
 This class is intended to hold what Bertini would need to produce a straight-line program for evaluation.
 */
class function
{
	std::string func;  //symbolic representation of function (straight from input file).
										 // this class is woefully incomplete.
};



/**
 \brief Comprehensive system randomization, based on deg.out.
 
 The system_randomizer is created based on the desired size to randomize down to, and the degrees of the functions, which are contained in the 'deg.out' file, which must pre-exist.
 
 This class does not keep track of the desired mp mode, and always populates all three randomizer matrices.
 
 This class is capable of randomizing for systems with a single hom_variable_group (set hom_var = 1), and a single variable_group.
 
 */
class system_randomizer
{
private:
	bool square_indicator; ///< a boolean indicating whether the system is square.
	bool setup_indicator; ///< a boolean indicating whether the randomizer is ready to use.
	
	mat_mp randomizer_matrix_full_prec; ///< holds the full precision randomizer matrix, from which we downsample when changing precision in MP mode.
	mat_mp randomizer_matrix_mp; ///< holds the randomizer matrix to current precision.
	mat_d randomizer_matrix_d; ///< holds the randomization matrix in double precision.
	
	int num_randomized_funcs; ///< the number of functions to which we randomize.
	int num_original_funcs; ///< the number of function from which we randomize.
	
	int max_base_degree; ///< the highest degree of any function occurring in the system.
	int max_degree_deficiency; ///< the greatest occuring deficiency in function degree relative to the max.
	std::vector<int> randomized_degrees; ///< a vector of integers keeping track of the degrees of the output functions.
	std::vector<int> original_degrees; ///< a vector of integers keeping track of the degrees of the input functions.
	
	std::vector<std::vector<int>> structure_matrix; ///< matrix of integers indicating the degree deficiency of each input function relative to the degree of the output functions.
	
	
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
	
	
	/**
	 output to a stream.  only really usable with std::cout, as it calls print_matrix_to_screen_matlab
	 
	 \param os the stream to put this text on.
	 \param s the system randomizer to write.
	 */
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
	
	/** 
	 \brief returns the number of original functions
	 \return int num_original_funcs
	 */
	int num_base_funcs() const
	{
		return num_original_funcs;
	}
	
	
	/**
	 \brief returns the number of randomized functions
	 \return int num_randomized_funcs
	 */
	int num_rand_funcs() const
	{
		return num_randomized_funcs;
	}
	
	
	/**
	 \brief returns the highest degree in the system
	 \return int max_base_degree
	 */
	int max_degree() const
	{
		return max_base_degree;
	}
	
	/**
	 \brief returns the greatest occurring deficiency
	 \return int max_degree_deficiency
	 */
	int max_deficiency() const
	{
		return max_degree_deficiency;
	}
	
	
	/**
	 \brief indicates whether the randomizer is square
	 \return bool square_indicator
	 */
	bool is_square() const
	{
		return square_indicator;
	}
	
	
	/**
	 \brief indicates whether the system_randomizer is ready to go.
	 \return bool setup_indicator
	 */
	bool is_ready() const
	{
		return setup_indicator;
	}
	
	/**
	 \brief gets the degree of function with index (loc).
	 \return int degree of original function with input index.
	 \param loc the index of the base function to get degree of.
	 */
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
	
	/**
	 \brief gets a pointer to the full precision randomizer matrix.
	 \return a pointer to the full precision randomizer matrix.
	 */
	mat_mp * get_mat_full_prec()
	{
		return &randomizer_matrix_full_prec;
	}
	
	/**
	 \brief gets a pointer to the double randomizer matrix.
	 \return a pointer to the double randomizer matrix.
	 */
	mat_d * get_mat_d()
	{
		return &randomizer_matrix_d;
	}
	
	
	/**
	 \brief gets a pointer to the current precision MP randomizer matrix.
	 \return a pointer to the current precision MP randomizer matrix.
	 */
	mat_mp * get_mat_mp()
	{
		return &randomizer_matrix_mp;
	}
	
	
	/**
	 \brief changes the precision of the system_randomizer
	 \param new_prec new precision.
	 */
	void change_prec(int new_prec);
	
	
	
	/**
	 \brief randomizes, by taking the input function values and jacobian, and multiplying the randomizer matrix.
	 \param randomized_func_vals output argument, set to R*f.
	 \param randomized_jacobian output jacobian, set to R*J.
	 \param func_vals input function values, probably produced by evaluating an SLP.
	 \param jacobian_vals input jacobian, probably produced by evaluating an SLP.
	 \param hom_var the coordinate of the single homogenizing variable.
	 */
	void randomize(vec_d randomized_func_vals, mat_d randomized_jacobian,
				   vec_d func_vals, mat_d jacobian_vals,
				   comp_d hom_var);
	
	
	/**
	 \brief randomizes, by taking the input function values and jacobian, and multiplying the randomizer matrix.
	 \param randomized_func_vals output argument, set to R*f.
	 \param randomized_jacobian output jacobian, set to R*J.
	 \param func_vals input function values, probably produced by evaluating an SLP.
	 \param jacobian_vals input jacobian, probably produced by evaluating an SLP.
	 \param hom_var the coordinate of the single homogenizing variable.
	 */
	void randomize(vec_mp randomized_func_vals, mat_mp randomized_jacobian,
				   vec_mp func_vals, mat_mp jacobian_vals,
				   comp_mp hom_var);
	
	
	/**
	 \brief sets up randomizer using 'deg.out'.
	 
	 This setup function parses "deg.out" for the degrees of the functions, and sets up the internals of this class object, for randomizing a system.
	 
	 \param num_desired_rows The number of output functions.  Must be bigger than the input num_funcs.
	 \param num_funcs The original number of functions in the system to be randomized.
	 */
	void setup(int num_desired_rows, int num_funcs);
	
	/** 
	 \brief sets up the temporaries.
	 set up the internal temporary variables.
	 */
	void setup_temps();
	
	
	/**
	 \brief single-target send
	 
	 send system_randomizer to a single target.
	 
	 \param target the id of the target.  
	 \param mpi_config container holding the mpi_config for the caller.
	 */
	void send(int target, parallelism_config & mpi_config);
	
	
	
	/**
	 \brief single-source receive
	 
	 receive system_randomizer from a single source.
	 
	 \param source the id of the source.
	 \param mpi_config container holding the mpi_config for the caller.
	 */
	void receive(int source, parallelism_config & mpi_config);
	
	
	/**
	 \brief collective MPI_COMM_WORLD broadcast send
	 
	 \see vertex_set::bcast_receive
	 
	 Send the system_randomizer to everyone in MPI_COMM_WORLD.
	 
	 \param mpi_config The current state of MPI, as represented in Bertini_real
	 */
	void bcast_send(parallelism_config & mpi_config);
	
	
	/**
	 \brief Collective MPI_COMM_WORLD broadcast receive
	 
	 Receive the system_randomizer from someone in MPI_COMM_WORLD.
	 
	 \param mpi_config The current state of MPI, as represented in Bertini_real
	 */
	void bcast_receive(parallelism_config & mpi_config);
	
	
protected:
	
	void init()
	{
		num_randomized_funcs = -1223;
		num_original_funcs = -976;
		
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

		
		
		if (other.is_ready()) {
			square_indicator = other.square_indicator;
			
			num_original_funcs = other.num_original_funcs;
			num_randomized_funcs = other.num_randomized_funcs;
			
			max_base_degree = other.max_base_degree;
			max_degree_deficiency = other.max_degree_deficiency;
			
			mat_cp_d(randomizer_matrix_d, other.randomizer_matrix_d);
			
			if (other.randomizer_matrix_mp->rows>0 && other. randomizer_matrix_mp->cols>0) {
				change_prec_mat_mp(randomizer_matrix_mp, mpf_get_prec(other.randomizer_matrix_mp->entry[0][0].r));
			}
			mat_cp_mp(randomizer_matrix_mp, other.randomizer_matrix_mp);
			
			mat_cp_mp(randomizer_matrix_full_prec, other.randomizer_matrix_full_prec);
			
			// these vectors are guaranteed to be nonempty by virtue of other being already setup.
			
			original_degrees.resize(0);
			for (auto iter=other.original_degrees.begin(); iter!=other.original_degrees.end(); ++iter) {
				original_degrees.push_back(*iter);
			}
			
			randomized_degrees.resize(0);
			for (auto iter=other.randomized_degrees.begin(); iter!=other.randomized_degrees.end(); ++iter) {
				randomized_degrees.push_back(*iter);
			}
			

			
			structure_matrix.resize(other.structure_matrix.size());
			for (unsigned int ii=0; ii<other.structure_matrix.size(); ++ii) {
				structure_matrix[ii].resize(0);
				for (auto jter = other.structure_matrix[ii].begin(); jter!= other.structure_matrix[ii].end(); ++jter) {
					structure_matrix[ii].push_back(*jter);
				}
			}
			
			setup_temps();
		}
		else
		{
			setup_indicator = false;
			
			randomized_degrees.resize(0);
			original_degrees.resize(0);
			structure_matrix.resize(0);
			
			num_original_funcs = num_randomized_funcs = max_base_degree = max_degree_deficiency = -1;
			
			change_size_mat_mp(randomizer_matrix_full_prec,0,0);
			change_size_mat_mp(randomizer_matrix_mp,0,0);
			change_size_mat_d(randomizer_matrix_d,0,0);
			randomizer_matrix_d->rows = randomizer_matrix_d->cols = randomizer_matrix_mp->rows = randomizer_matrix_mp->cols = randomizer_matrix_full_prec->rows = randomizer_matrix_full_prec->cols = 0;
			
			
			
		}
		
	}
	
	

};//re: randomizer class





/**
 \brief base class for holding a set of vec_mp's.
 
 This class gives a way to commit vec_mp's into a class object.  The major data members are pts_mp and num_pts.  The most important member function is add_point(p)
 */
class point_holder
{
	
public:
	
	vec_mp *pts_mp; ///< an array of vec_mp, which are structs and require manual initialization and clearing.
	int num_points; ///< the number of stored points.
	
	
	/**
	 copy all the points from one point_holder to this one.
	 
	 \param other an input point_holder from which to copy all the points.
	 */
	void copy_points(const point_holder & other) {
		
        for (int ii=0; ii<other.num_points; ii++)
			add_point(other.pts_mp[ii]);
    }
	
	
	/**
	 indicate whether the point_holder has at zero points in it.
	 
	 \return a boolean, true if num_points==0, false otherwise
	 */
	inline bool has_no_points() const
	{
		if (num_points==0) {
			return true;
		}
		else{
			return false;
		}
	}
	
	
	/**
	 indicate whether the point_holder has at least one point in it.
	 
	 \return a boolean, true if num_points>0, false if num_points==0
	 */
	inline bool has_points() const
	{
		if (num_points==0) {
			return false;
		}
		else{
			return true;
		}
	}
	
	
	/**
	 \brief reset to empty container.
	 
	 reset to empty container
	 */
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
	
	
	/**
	 \brief add a point to the point_holder
	 
	 \return the index of the new point.
	 \param new_point the point to add.
	 */
	int add_point(vec_mp new_point);
	
	
	
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




/**
 \brief a way to hold vec_mp's as patches.
 
 This class offers a way to hold a set of vec_mp's as patch_mp members in a class with automated initting and clearing.
 */
class patch_holder
{

public:
	
	vec_mp *patch_mp;   ///< an array of patch_mp's
	int num_patches; ///< the number of patches stored in this object.
	
	
	
	/**
	 \brief copy all the patches from another patch_holder
	 
	 Copy all the stored mp patches from another patch_holder object, without testing for uniqueness.
	 
	 \param other the patch_holder from which to copy.
	 */
	void copy_patches(const patch_holder & other) {
		
        for (int ii=0; ii<other.num_patches; ii++)
			add_patch(other.patch_mp[ii]);
    }
	
	
	
	/**
	 \brief resets the patch holder to empty.
	 
	 Reset the patch_holder to 0 patches.
	 */
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
	
	
	/**
	 \brief add a patch to this object.
	 
	 \return the index of the patch just added.
	 \param new_patch the new patch to add.
	 */
	int add_patch(vec_mp new_patch);
	
	
	
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


/**
 \brief class for holding vec_mp's as linears in an automated object.
 
 This class offers automated collection of vec_mp's as L_mp.  This is necessary because the vec_mp type must be initted and cleared manually.
 */
class linear_holder
{
	
public:
	
	vec_mp *L_mp;  ///< a pointer array of vec_mp's as linears.
	int num_linears; ///< the number of linears in this collection.
	
	
	/**
	 \brief copies all linears from another linear_holder.
	 
	 \param other the linear_holder from which to copy all the linears.
	 */
	void copy_linears(const linear_holder & other) {
        for (int ii=0; ii<other.num_linears; ii++)
			add_linear(other.L_mp[ii]);
    }
	
	
	/**
	 reset this object to an empty state.
	 */
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
	
	
	/**
	 \brief add a linear to this collection.
	 
	 \return the index of the newly added linear.
	 \param new_linear the vec_mp linear to add.
	 */
	int add_linear(vec_mp new_linear);
	
	
	
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


/** 
 \brief holds the names of variables
 */
class name_holder
{
public:
	
	std::vector< std::string > variable_names; ///< the names.
	
	
	void reset_names()
	{
		variable_names.resize(0);
	}
	
	void cp_names(const name_holder & nomnom)
	{
		this->variable_names = nomnom.variable_names;
	}
	
	
	/**
	 \brief read variable names from names.out
	 
	 Reads variable names from names.out, which must exist before calling this function.
	 
	 \param num_vars the number of variable names to read.
	 */
	void get_variable_names(int num_vars)
	{
		
		variable_names.resize(num_vars);
		
		std::ifstream fin("names.out");
		
		for (int ii=0; ii<num_vars; ++ii){
			fin >> this->variable_names[ii];
		}
		fin.close();
		
	}
	
	
};


/**
 \brief witness set holds points, patches, and linears, with the names of the variables.
 
 The witness set class collects points, patches, and linears into one object.  It offers methods for sorting for real points only, for sorting to contain only unique points.
 
 A witness set gets two numbers of variables, one is the total number appearing in it, and the other [more importantly] is the numebr of natural variables contained therein.  The witness set is assumed to be in correspondence to a variable group with a single leading homogenizing variable, and automatically dehomogenizes points for uniqueness and reality testing.
 */
class witness_set : public patch_holder, public linear_holder, public point_holder, public name_holder
{
	
public:
	
	//begin data members
	

	int dim;
	int comp_num;
	int incidence_number;
	
	int num_variables;
	int num_natural_vars;
	


	

	
	boost::filesystem::path input_filename;
	function input_file;
	
	
	// end data members
	
	
	
	// overloaded operators
	
	
	// default constructor
	
	witness_set(int nvar)
	{
		init();
		num_variables = num_natural_vars = nvar;
	};
	
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
		this->num_natural_vars = 0;
		
		
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
		this->num_natural_vars = other.num_natural_vars;
	}

    

    

	
    int num_synth_vars() const
    {
        return num_variables - num_natural_vars;
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
	
	

	
	

	

	
	/**
	 \brief merges another witness set into this, not checking for uniqueness at all.
	 
	 Should you want to straight-up merge the contents of two witness sets with the same number of natural variables, you may, using this function.  All the points, linears, and patches will be copied from the input into the existing one on which you call this function.
	 
	 \param W_in The witness set containing data you want to copy.
	 */
	void merge(const witness_set & W_in);
	
	

	
	/**
	 \brief prints some information about the witness set to the screen
	 
	 Print variable information, linears, and patches to screen  
	 This is potentially a very large amount of data depending on the set, and should be done sparingly
	 */
	void print_to_screen() const;
	
	/**
	 \brief print the witness set into a file, which can be read back in to the same format.
	 
	 Print the witness set to a file of filename
	 
	 the format is:
	 
	 [
	 num_points dim comp_num
	 
	 point 1 - in standard bertini format
	 
	 point 2 - in standard bertini format
	 
	 ...
	 
	 last point - in standard bertini format
	 
	 num_linears num_variables
	 
	 linear 1 - in standard bertini point/vec format
	 
	 ...
	 
	 linear last - in standard bertini point/vec format.
	 
	 num_patches num_variables \todo this is incorrect.
	 
	 patch 1 - in standard bertini point/vec format
	 
	 ...
	 
	 patch last - in standard bertini point/vec format.
	 ]
	 
	 
	 \param filename the name of the file to which to write.
	 */
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
	
	/**
	 read patches from a suitably set up file.
	 
	 \param filename the name of the file to read.
	 */
	void read_patches_from_file(boost::filesystem::path filename);
	
	/**
	 write all the points to a file, without dehomogenizing.
	 
	 \param filename the name of the file to which to write.
	 */
	void write_homogeneous_coordinates(boost::filesystem::path filename) const;
	
	/**
	 write all the points to a file, first dehomogenizing.
	 
	 \param filename the name of the file to which to write.
	 */
	void write_dehomogenized_coordinates(boost::filesystem::path filename) const;
	
	
	/**
	 individual send, relative to MPI_COMM_WORLD
	 
	 \param target the ID target of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void send(parallelism_config & mpi_config, int target);
	
	/**
	 individual receive, relative to MPI_COMM_WORLD
	 
	 \param source the ID source of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void receive(int source, parallelism_config & mpi_config);
};






















// CURVE CELL DECOMP DATA TYPES


/**
 a bertini_real vertex, a 0-cell.  contains a point, its projection values, its type, whether its been removed, and an index into a set of filenames contained in the vertex_set.
 
 \todo remove the metadata from this, and instead track it in the vertex set, much like the solver_output
 */
class vertex
{
public:
	
	point_mp pt_mp; ///< the main data for this class -- the point.
	
	
	vec_mp  projection_values; ///< a vector containing the projection values.
	
	int type;  ///< See enum.
	int removed; ///< boolean integer whether the vertex has been 'removed' by a merge process.
	int input_filename_index; ///< index into the vertex_set's vector of filenames.
	
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
	
	/**
	 /brief prints the vertex to the screen
	 Prints the vertex to the screen
	 */
	void print() const
	{
		print_point_to_screen_matlab(pt_mp,"point");
		print_point_to_screen_matlab(projection_values,"projection_values");
		std::cout << "type: " << type << std::endl;
	}
	
	
	/**
	 \brief sets the point.
	 
	 Set the vertex's point to the input.
	 \param new_point the input point to be set.
	 */
	void set_point(const vec_mp new_point);
	
	
	/**
	 \brief single target mpi send.
	 
	 Send the vertex to a single target.
	 
	 \see vertex::receive
	 
	 \param target The MPI ID of the target for this send
	 \param mpi_config current mpi settings
	 */
	void send(int target, parallelism_config & mpi_config);
	
	
	/**
	 \brief single source receive.
	 
	 Receive a vertex from a single source.
	 
	 \see vertex::send
	 
	 \param source The MPI ID of hte source.
	 \param mpi_config current mpi settings
	 */
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
 \brief the main structure for storing vertices.  
 
 The vertex_set is bertini_real's main method for storing data.  We essentially construct a graph of vertices, consisting of edges and faces.
 
 there are methods in place to add vertices, and perform lookups.
 */
class vertex_set
{
public:
	
	vec_mp *projections; ///< a pointer array of projection vectors.
	int num_projections; ///< the number of projections.  this should match the dimension of the object being decomposed.
	int curr_projection; ///< the projection currently being used.
	
	int curr_input_index; ///< the index of the current input file.
	std::vector< boost::filesystem::path > filenames; ///< the set of filenames from which vertices arise.
	
	std::vector<vertex> vertices;  ///< the main storage of points in the decomposition.
	int num_vertices; ///< the number of vertices found so far.
	
	int num_natural_variables;  ///< the number of natural variables appearing in the problem to solve.
	
	
	
	mpf_t abs;
	mpf_t zerothresh;
	comp_mp diff;
	vec_mp checker_1;
	vec_mp checker_2;
	
	
	
	void print_to_screen(); ///< operator for displaying information to screen
	
	
	/**
	 \brief add a new vertex to the set.
	 
	 \param new_vertex
	 \return the index of the added vertex
	 */
	int add_vertex(const vertex new_vertex);
	
	
	/**
	 \brief create a vertex_set from a file.
	 
	 Read in a vertex_set from a file.
	 
	 \param INfile the file to parse and store in a vertex_set
	 \return the number of vertices read in.
	 */
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
	
	
	
	/**
	 sets the number of variables for the vertex_set.
	 
	 \param num_vars the number of {\em natural} variables
	 */
	void set_num_vars(int num_vars)
	{
		this->num_natural_variables = num_vars;
		
		change_size_vec_mp(checker_1, num_vars);
		change_size_vec_mp(checker_2, num_vars);
		checker_1->size = checker_2->size = num_vars;
	}
	
	
	
	/**
	 \brief write vertex_set to a file, readable by bertini_real again.
	 
	 
	 write the vertex_set to a file.
	 \see setup_vertices
	 
	 Output vertex structure as follows:
	 
	 [
	 num_vertices num_projections num_natural_variables filenames.size()
	 
	 the projections, as bertini points
	 
	 the names of the files as
	 length_of_name name   pairs
	 
	 then the points as
	 
	 
	 num_variables in point \\
	 coordinates
	 
	 num_projection_coordinates \\
	 projection values 
	 
	 filename_index 
	 
	 type
	]
	 
	 
	 \param outputfile the name of the file to write the vertex_set to.
	 */
	void print(boost::filesystem::path outputfile) const;
	
	
	
	/**
	 set the name of the current input file.  while it is set to this, all added vertices will inherit the index of this name.
	 
	 \return the index of the set filename.
	 \param el_nom the name of the file
	 */
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
	
	
    

    /**
	 \brief find the index of a point.
	 
	 find the index of a point.
	 first, we check the active points, then the inactive, then give up.  the method dehomogenizes the points, and checks only the natural variables
	 
	 \return the index of the testpoint, or -1 if it is not found.
	 \param testpoint the mp point to find.
	 */
    int search_for_point(vec_mp testpoint);
	
	
	/**
	 \brief find the index of a point among active (non-removed) points only.
	 
	 find the index of a point.
	 
	 
	 \return the index of the testpoint, or -1 if it is not found.
	 \param testpoint the mp point to find.
	 */
    int search_for_active_point(vec_mp testpoint);
	
	
	
	/**
	 \brief find the index of a point among removed points only.
	 
	 find the index of a point.
	 
	 
	 \return the index of the testpoint, or -1 if it is not found.
	 \param testpoint the mp point to find.
	 */
    int search_for_removed_point(vec_mp testpoint);
    
	
	/**
	 \brief Compute projections values, and midpoints, of a set of points with respect to a projection.
	 
	 This function computes \f$\pi(x)\f$ for each of the points \f$x\f$ in witness_set W.  Then it sorts them, and computes averages.
	 
	 The output is stored in crit_downstairs and midpoints_downstairs, both pre-initialized vec_mp's.  This function also produces a std::vector<int> named index_tracker which contains the sorting of W according to \f$\pi(W)\f$.
	 
     \param W witness				set containing points of which we wish to retrieve projections values.
     \param crit_downstairs			the projection values of the input witness set, sorted for uniqueness and increasingness.
     \param midpoints_downstairs	the bisection of each interval in crit_downstairs.
     \param index_tracker			the indices of the points in W.
     \param pi						the projection we are retrieving projection values with respect to.
     \return the integer SUCCESSFUL.
     */
    int compute_downstairs_crit_midpts(const witness_set & W,
                                       vec_mp crit_downstairs,
                                       vec_mp midpoints_downstairs,
                                       std::vector< int > & index_tracker,
                                       vec_mp pi);
    
    
	/**
	 
	 sets the value of the [current] projection for each vertex which has index in the set of relevant indices.
	 
	 \param relevant_indices set of vertex indices
	 \param new_value the new value you want to set the projection value to
	 \return a vector of indices for which this operation failed, because the new and old values were too far away from each other.
	 */
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value);
	
	/**
	 
	 sets the value of the [proj_index] projection for each vertex which has index in the set of relevant indices.
	 
	 \param relevant_indices set of vertex indices
	 \param new_value the new value you want to set the projection value to
	 \param proj_index the index of the projection you want to assert.
	 \return a vector of indices for which this operation failed, because the new and old values were too far away from each other.
	 */
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index);
	
	
	/**
	 
	 
	 set the current projection.
	 
	 \return the index of the projection
	 \param new_proj the projection to set as current
	 */
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
	
    
	/**
	 \brief query the index of a projection.
	 
	 Want to find out the index of a projection in this vertex_set? pass it into this function.
	 
	 \param proj the projection to query
	 \return the index of the projection, or -1 if it doesn't exist.
	 */
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
    
    
	/**
	 \brief add a projection to the set, and get its index.
	 
	 Add a projection to the set, and get its index.  Does not test for uniqueness of the projection, assumes it is not in there yet.
	 
	 \return the index of the added projection
	 \param proj the projection to add.
	 
	 */
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
	
	
	
	/**
	 \brief single target mpi send
	 
	 Send a vertex_set to a single target process.
	 
	 \param target who to send to
	 \param mpi_config the current mpi configuration
	 */
	void send(int target, parallelism_config & mpi_config);
	
	
	
	/**
	 \brief single source mpi receive
	 
	 Receive a vertex_set from a single source.
	 
	 \param source the source of the receive
	 \param mpi_config the current mpi configuration
	 */
	void receive(int source, parallelism_config & mpi_config);
	
	
	/**
	 reset the set to empty.
	 */
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






typedef std::pair<int,int> witness_set_index;


/**
 \brief metadata for witness points, for the witness_data class.
 
 
 */
class witness_point_metadata
{
public:
	int dimension;
	
	int corank, typeflag, multiplicity, component_number, deflations_needed;
	double condition_number, smallest_nonsing_value, largest_nonsing_value;
	

	/**
	 read the metadata from the witness_data file.  call only at the appropriate point.
	 
	 \param IN a pointer to an open file from which to read.  must be set to the correct point
	 */
	void set_from_file(FILE *IN)
	{
		fscanf(IN,"%lf %d %lf %lf %d %d %d %d",
			   &condition_number,
			   &corank,
			   &smallest_nonsing_value,
			   &largest_nonsing_value,
			   &typeflag,
			   &multiplicity,
			   &component_number,
			   &deflations_needed);
	}
	
	
	/**
	 \brief constructor, setting the dimension in the process.
	 constructor, setting the dimension in the process.
	 
	 \param dim the dimension to set to.
	 */
	witness_point_metadata(int dim){dimension = dim;}
	
	
	witness_point_metadata(){};
	
	witness_point_metadata(const witness_point_metadata& other){copy(other);}
	
	witness_point_metadata& operator=( const witness_point_metadata& other)
	{
		copy(other);
		return *this;
	}
	
	void copy(const witness_point_metadata& other)
	{
		dimension = other.dimension;
		corank = other.corank;
		typeflag = other.typeflag;
		multiplicity = other.multiplicity;
		component_number = other.component_number;
		deflations_needed = other.deflations_needed;
		condition_number = other.condition_number;
		smallest_nonsing_value = other.smallest_nonsing_value;
		largest_nonsing_value = other.largest_nonsing_value;
	}
	
	
};


/**
 \brief metadata for linears read in from the witness_data file.
 
 metadata for linears read in from the witness_data file.
 */
class witness_linear_metadata
{
	
public:
	
	int dimension;
	
	
	witness_linear_metadata(int dim){dimension = dim;}
	
	witness_linear_metadata(){};
	witness_linear_metadata(const witness_linear_metadata& other)
	{
		dimension = other.dimension;
	}
	
	witness_linear_metadata& operator=(const witness_linear_metadata &other)
	{
		dimension = other.dimension;
		return *this;
	}
};

/**
 \brief metadata for patches read in from witness_data
 
 metadata for patches read in from witness_data
 */
class witness_patch_metadata
{
	
public:
	
	witness_patch_metadata(){};
	
	
	/**
	 constructor, setting the dimension in the process
	 \param dim the dimension to set.
	 */
	witness_patch_metadata(int dim){dimension = dim;}

	
	witness_patch_metadata(const witness_patch_metadata& other)
	{
		dimension = other.dimension;
	}
	
	witness_patch_metadata& operator=(const witness_patch_metadata &other)
	{
		dimension = other.dimension;
		return *this;
	}
	
	
	int dimension;
};




/**
 \brief a nearly complete class for storing bertini's witness_data file.
 
 This class reads in witness_data, and produces witness_sets based on user's choice.
 */
class witness_data : public patch_holder, public linear_holder, public point_holder
{
	
private:
	
	std::vector< witness_point_metadata > point_metadata;
	std::vector< witness_linear_metadata > linear_metadata;
	std::vector< witness_patch_metadata > patch_metadata;
	
	std::vector<int> nonempty_dimensions;
	std::map<int,std::map<int,int> > dimension_component_counter;
	std::map<int,std::map<int,std::vector<int>>> index_tracker;
	
public:
	
	int num_variables;
	
	
	
	void reset()
	{
		reset_points();
		reset_linears();
		reset_patches();
		
		point_metadata.resize(0);
		linear_metadata.resize(0);
		patch_metadata.resize(0);
		
	}
	
	
	
	void populate(); ///< fills this object with the sets in witness_data.
	
	
	
	/**
	 \brief outermost method for choosing a witness set to construct.
	 
	 This function uses information stored in BR_configuration to construct a witness set.
	 
	 \param options the current program state.
	 \return the chosen witness set.  may be empty.
	 */
	witness_set choose(BR_configuration & options);
	witness_set best_possible_automatic_set(BR_configuration & options);
	witness_set choose_set_interactive(BR_configuration & options); // lets the user choose a set, and returns a copy of it.
	
	witness_set form_specific_witness_set(int dim, int comp);
	
	
	
//	friend std::ostream & operator<<(std::ostream &os, witness_data & c)
//	{
//		for (auto iter=dimension_component_counter.begin(); iter!=dimension_component_counter.end(); ++iter) {
//			os << "dimension " << iter->first << std::endl;
//			for (auto jter=iter->second.begin(); jter!=iter->second.end(); ++jter) {
//				os << "\tcomponent " << jter->first << " degree " << jter->second << std::endl;
//			}
//		}
//		
//		
//		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
//			os << "dimension " << iter->first << " ";
//			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
//				os << "component " << jter->first << std::endl;
//				for (auto kter = jter->second.begin(); kter!=jter->second.end(); kter++) {
//					os << *kter << " ";
//				}
//				os << std::endl;
//			}
//			os << std::endl;
//		}
//		
//		
//		
//		
//		
//		return os;
//	}
	
	/**
	 print the witness_data to screen
	 */
	void print()
	{
		
		std::cout << "the nonempty dimensions:" << std::endl;
		for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
			std::cout << *iter << " ";
		}
		std::cout << std::endl;
		
		
		
		for (auto iter=dimension_component_counter.begin(); iter!=dimension_component_counter.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter=iter->second.begin(); jter!=iter->second.end(); ++jter) {
				std::cout << "\tcomponent " << jter->first << " degree " << jter->second << std::endl;
			}
		}
		
		
		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
				std::cout << "\tcomponent " << jter->first << ": ";
				for (auto kter = jter->second.begin(); kter!=jter->second.end(); kter++) {
					std::cout << *kter << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		
		for (int ii=0; ii<num_points; ii++) {
			print_point_to_screen_matlab(pts_mp[ii],"p");
		}
	}
	
	
	
//	witness_data()
//	{
//		init();
//	}
	
	
//	//copy operator.
//	witness_data(const witness_data & other)
//	{
//		init();
////		copy(other);
//	}
//	
//	
//	// assignment
//	witness_data& operator=( const witness_data& other)
//	{
//		reset();
////		copy(other);
//		return *this;
//	}
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
private:
	
	void add_linear_w_meta(vec_mp lin, const witness_linear_metadata & meta)
	{
		add_linear(lin);
		linear_metadata.push_back(meta);
	}
	
	void add_patch_w_meta(vec_mp pat, const witness_patch_metadata & meta)
	{
		add_patch(pat);
		patch_metadata.push_back(meta);
	}
	
	
	int add_solution(vec_mp pt, const witness_point_metadata & meta)
	{
		
		int ind = add_point(pt);
		point_metadata.push_back(meta);
		return ind;
	}
	
	
	
	
//	void init(){};
//	
//	void copy(const witness_data & other)
//	{
//		
//	}
	
	
	
	
};





class cell
{
	
private:

	
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
 \brief 1-cell.
 
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











/**
 \brief base decomposition class.  curves and surfaces inherit from this.
 
 The decomposition class holds the basic information for any dimensional decomposition -- a witness set which generated it, the number of variables, the dimension, which component it represents, the projections, the randomizer, the sphere, the input file name.
 
 */
class decomposition : public patch_holder
{

public:
	std::map< int , int > counters;
	std::map< int , std::vector< int > > indices;
	
	witness_set W;
	
	int num_variables; ///< the number of variables in the decomposition
	int dimension; ///< the dimension of the decomposition
	int component_num; ///< the component number.
	
	int num_curr_projections; ///< the number of projections stored in the decomposition.  should match the dimension when complete.
	vec_mp	*pi; ///< the projections used to decompose.  first ones are used to decompose nested objects.
	
	system_randomizer * randomizer; ///< the randomizer for the decomposition.

	vec_mp sphere_center; ///< the center of the sphere.
	comp_mp sphere_radius; ///< the radius of the sphere.
	bool have_sphere_radius; ///< indicates whether the decomposition has the radius set, or needs one still.
	
	boost::filesystem::path input_filename; ///< the name of the text file in which the system resides.
//	function input_file;
	
	
	
    
/**
 \brief add a projection vector to the decomposition
 
 Decompositions in Bertini_real are computed with respect to linear projections, which are stored as vectors. They are repeated in each decomposition.  
 
 \param proj the projection to add.
 */
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
	
	
	/**
	 \brief commit a set of points to the vertex set, associating them with the input file for this decomposition.
	 
	 \todo rename this function to something more accurately descriptive
	 
	 \return the number 0.  this seems pointless.
	 \param W the witness set containing the points to add.
	 \param add_type the type the points will inherit.  {\em e.g.} CRITICAL
	 \param V the vertex set to add the points to.
	 */
    int add_witness_set(const witness_set & W, int add_type, vertex_set & V);
    
	
	/**
	 \brief Find the index of a testpoint.
	 
	 Search the vertex set passed in for testpoint.  This happens relative to this decomposition for no reason whatsoever.
	 \todo change this to be a property of the vertex set
	 
	 \return the index of the point, or -1 if not found.
	 \param V the vertex set in which to search
	 \param testpoint the point for which to search
	 */
	int index_in_vertices(vertex_set &V,
                          vec_mp testpoint);
	
	
	/**
	 \brief search for a testvertex, and add it to vertex set, and this decomposition, if not found.
	 
	 Search the passed in vertex set for the testvertex -- and add it to the vertex set, and its index to the decomposition if not found.
	 
	 \return the index of the point, or -1 if not found.
	 \param V the vertex set in which to search
	 \param vert a vertex with point for which to search.
	 */
	int index_in_vertices_with_add(vertex_set &V,
                                   vertex vert);
	
	/** 
	 set up the base decomposition class from a file
	 
	 \return the number 0. seems useless.
	 \param INfile the name of the file to read from
	 */
	int setup(boost::filesystem::path INfile);
	
	
	/**
	 base method for printing decomposition to file.
	 
	 \param outputfile the name of the file to which to write.
	 */
	virtual void print(boost::filesystem::path outputfile);
	
	
	
	/**
	 read the bounding sphere from a file containing the radius and center coordinates.
	 
	 format is:
	 
	 [radius
	 
	 x_0
	 x_1
	 ...
	 x_n]
	 
	 
	 \return SUCCESSFUL value.
	 \param bounding_sphere_filename the name of the file.
	 */
	int read_sphere(const boost::filesystem::path & bounding_sphere_filename);
	
	
	/**
	 \brief set the sphere for this decomposition according to the input witness set.
	 
	 Because we need to capture parts of the component which go to inifinity, we intersect the component with a sphere containing all the critical points, which are the input to this method.  
	 
	 This method computes the centroid of the set -- which becomes the center of the sphere -- and the distance from the centroid to the outermost critical point -- 3 times which becomes the radius of the sphere.
	 
	 \param W_crit the witness set containing the critical points to capture inside the sphere
	 */
	void compute_sphere_bounds(const witness_set & W_crit);
	
	
	/**
	 \brief copy the bounds from another decomposition
	 
	 Sub-decompositions will often want to inherit the bounding sphere of another.  This method lets you copy from one to another.
	 
	 
	 \param other the decomposition which already holds sphere bounds.
	 */
	void copy_sphere_bounds(const decomposition & other)
	{
		if (!other.have_sphere_radius) {
			std::cout << "trying to copy sphere bounds from a decomposition which does not have them set!" << std::endl;
			br_exit(72471);
		}
		
		set_mp(this->sphere_radius, other.sphere_radius);
		vec_cp_mp(this->sphere_center, other.sphere_center);
		this->have_sphere_radius = true;
	}
	
	
	/**
	 \brief the main way to print a decomposition to a file.
	 
	 This method backs up the existing folder to one siffixed with "_bak", and creates a new folder with the correct name, to which it prints the decomposition in text file format.
	 
	 \param base the base folder name to print the decomposition.
	 */
	void output_main(const boost::filesystem::path base);
	
	
	
	
	
	
	/**
	 reset decomposition to empty.
	 */
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
	
	
	/**
	 \brief single target MPI send.
	 
	 Send base decomposition to another process.
	 
	 \param target the ID of the worker to which to send.
	 \param mpi_config the current configuration of MPI
	 */
	void send(int target, parallelism_config & mpi_config);
	
	/**
	 \brief single source MPI receive.
	 
	 Receive base decomposition from another process.
	 
	 \param source the ID of the process from which to receive
	 \param mpi_config the current configuration of MPI
	 */
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













/**
 \brief terminal-control of colors.
 
 namespace of color functions, using terminal controls.
 
 e.g. \033[0;30m for black.
 */
namespace color {
	
	/** get the integer associated with the name of a color 
	 
	 k - black - 30
	 r - red - 31
	 g - green - 32
	 y - brown - 33
	 b - blue - 34
	 m - magenta - 35
	 c - cyan - 36
	 l - lightgray - 37
	 
	 \return integer corresponding to color
	 \param c the single-character string of the color.
	 */
	int color_to_int(const std::string c);
	
	/**
	 set the text to bold
	 
	 \033[1;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string bold(std::string new_color);
	
	/**
	 set the text to darker color
	 
	 \033[2;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string dark(std::string new_color);
	
	/**
	 set the text to underline
	 
	 \033[4;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string underline(std::string new_color);
	
	/**
	 set the background to a color indicated by a name
	 
	 \033[7;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string background(std::string new_color);
	
	/**
	 set the text to strikethrough
	 
	 \033[9;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string strike(std::string new_color);
			
	
	
	
	/**
	 set the text to whatever the console believes is the default
	 
	 \033[0m
	 
	 \return string to print to cout to control color
	 */
	std::string console_default();
	
	/**
	 set the text to black
	 
	 \033[0;30m
	 
	 \return string to print to cout to control color
	 */
	std::string black();
	
	/**
	 set the text to red
	 
	 \033[0;31m
	 
	 \return string to print to cout to control color
	 */
	std::string red();
	
	/**
	 set the text to green
	 
	 \033[0;32m
	 
	 \return string to print to cout to control color
	 */
	std::string green();

	
	/**
	 set the text to brown
	 
	 \033[0;33m
	 
	 \return string to print to cout to control color
	 */
	std::string brown();
	
	/**
	 set the text to blue
	 
	 \033[0;34m
	 
	 \return string to print to cout to control color
	 */
	std::string blue();
	
	/**
	 set the text to magenta
	 
	 \033[0;35m
	 
	 \return string to print to cout to control color
	 */
	std::string magenta();
	
	/**
	 set the text to cyan
	 
	 \033[0;36m
	 
	 \return string to print to cout to control color
	 */
	std::string cyan();
	
	/**
	 set the text to gray
	 
	 \033[0;37m
	 
	 \return string to print to cout to control color
	 */
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

