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

#include <memory>

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
value_type map_lookup_with_default(const  std::map <key_type,value_type> & mc_mapperson, const key_type & lookup_key, const value_type& default_value )
{
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
	
	std::vector<std::vector<int> > structure_matrix; ///< matrix of integers indicating the degree deficiency of each input function relative to the degree of the output functions.
	
	
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
	

protected:
	
	vec_mp *pts_mp_; ///< an array of vec_mp, which are structs and require manual initialization and clearing.
	size_t num_pts_; ///< the number of stored points.
	
public:
	
	vec_mp * point(unsigned int index) const
	{
		if (index < num_pts_) {
			return &pts_mp_[index];
		}
		else
		{
			br_exit(-9834718);
			return NULL;
		}
		
	}
	
	
	/**
	 get the number of points in the point holder
	 
	 \return the current number of points
	 */
	size_t num_points() const
	{
		return num_pts_;
	}
	
	
	/**
	 copy all the points from one point_holder to this one.
	 
	 \param other an input point_holder from which to copy all the points.
	 */
	void copy_points(const point_holder & other)
	{
		
        for (unsigned int ii=0; ii<other.num_pts_; ii++)
			add_point(other.pts_mp_[ii]);
    }
	
	
	/**
	 indicate whether the point_holder has at zero points in it.
	 
	 \return a boolean, true if num_points==0, false otherwise
	 */
	inline bool has_no_points() const
	{
		if (num_pts_==0) {
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
		if (num_pts_==0) {
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
		for (unsigned int ii =0; ii<num_pts_; ii++)
			clear_vec_mp(pts_mp_[ii]);
		
		if (num_pts_>0) {
			free(pts_mp_);
		}
		
		num_pts_ = 0;
		pts_mp_ = NULL;
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
		pts_mp_ = NULL; // potential data loss if used improperly, which i may.
		num_pts_ = 0;
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

protected:
	
	vec_mp *patch_mp_;   ///< an array of patch_mp's
	size_t num_patches_; ///< the number of patches stored in this object.
	
public:
	
	
	/**
	 \brief  query how many patches are being held.
	 query how many patches are being held.
	 
	 \return the number of patches.
	 */
	inline size_t num_patches() const
	{
		return num_patches_;
	}
	
	
	/**
	\brief get a pointer to the patch at index
	 
	 \return a pointer to the i^th patch.
	 \param index The index of the patch you want
	 */
	inline vec_mp * patch(unsigned int index) const
	{
		return &patch_mp_[index];
	}
	
	
	/**
	 \brief copy all the patches from another patch_holder
	 
	 Copy all the stored mp patches from another patch_holder object, without testing for uniqueness.
	 
	 \param other the patch_holder from which to copy.
	 */
	void copy_patches(const patch_holder & other) {
		
        for (unsigned int ii=0; ii<other.num_patches_; ii++)
			add_patch(other.patch_mp_[ii]);
    }
	
	
	
	/**
	 \brief resets the patch holder to empty.
	 
	 Reset the patch_holder to 0 patches.
	 */
	void reset_patches()
	{
		for (unsigned int ii =0; ii<num_patches_; ii++)
			clear_vec_mp(patch_mp_[ii]);
		
		if (num_patches_>0) {
			free(patch_mp_);
		}
		
		num_patches_ = 0;
		patch_mp_ = NULL;
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
		this->patch_mp_ = NULL;
		this->num_patches_ = 0;
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
	
protected:
	
	vec_mp *L_mp_;  ///< a pointer array of vec_mp's as linears.
	size_t num_linears_; ///< the number of linears in this collection.
	
	
	
public:
	
	/**
	 \brief  query how many patches are being held.
	 query how many patches are being held.
	 
	 \return the number of patches.
	 */
	inline size_t num_linears() const
	{
		return num_linears_;
	}
	
	
	/**
	 \brief get a pointer to the patch at index
	 
	 \return a pointer to the i^th patch.
	 \param index The index of the patch you want
	 */
	inline vec_mp * linear(unsigned int index) const
	{
		return &L_mp_[index];
	}
	
	
	/**
	 \brief copies all linears from another linear_holder.
	 
	 \param other the linear_holder from which to copy all the linears.
	 */
	void copy_linears(const linear_holder & other) {
        for (unsigned int ii=0; ii<other.num_linears_; ii++)
			add_linear(other.L_mp_[ii]);
    }
	
	
	/**
	 reset this object to an empty state.
	 */
	void reset_linears()
	{
		for (unsigned int ii =0; ii<num_linears_; ii++)
			clear_vec_mp(L_mp_[ii]);
		
		if (num_linears_>0) {
			free(L_mp_);
		}
		
		num_linears_ = 0;
		L_mp_ = NULL;
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
		this->L_mp_ = NULL;
		this->num_linears_ = 0;
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

	
private:
	
	std::vector< std::string > variable_names; ///< the names.
	
public:
	
	
	/**
	 \brief get how many names there are.
	 
	 \return the number of stored names
	 */
	inline size_t num_var_names() const
	{
		return variable_names.size();
	}
	
	
	/**
	 \brief get the ith variable name
	 */
	inline std::string name(unsigned int index) const
	{
		if (index>=variable_names.size()) {
			throw std::out_of_range("trying to get variable name, requested index is out of bounds");
		}
		else
		{
			return variable_names[index];
		}
		
	}
	
	
	
	/**
	 \brief reset the names to empty.
	 */
	void reset_names()
	{
		variable_names.resize(0);
	}
	
	void copy_names(const name_holder & nomnom)
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
	
protected:
	
	//begin data members
	

	int dim_;
	int comp_num_;
	int incid_num_;
	
	int num_vars_;
	int num_natty_vars_;
	


	boost::filesystem::path input_filename_;
	function input_file_;
	
	
	// end data members
	
	
public:
	
	/**
	 \brief get the name of the bertini input file for this witness set
	 
	 \return the path of the file
	 */
	inline boost::filesystem::path input_filename() const
	{
		return input_filename_;
	}
	
	/**
	 \brief set the name of the bertini input file
	 
	 \param new_input_filename The new name of the file
	 */
	void set_input_filename(boost::filesystem::path new_input_filename)
	{
		input_filename_ = new_input_filename;
	}
	
	/**
	 \brief get the dimension of the set represented by the witness set.
	 
	 \return the integer dimension of the component.
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	/**
	 \brief set the dimension of the witness set
	 
	 \param new_dim the dimension of the set
	 */
	void set_dimension(int new_dim)
	{
		dim_ = new_dim;
	}
	
	
	
	
	
	
	/**
	 \brief get the component number of the set represented by the witness set.
	 
	 \return the index the component.
	 */
	inline int component_number() const
	{
		return comp_num_;
	}
	
	
	/**
	 \brief sets the component number
	 
	 \param new_comp_num The new component number to set.
	 */
	void set_component_number(int new_comp_num)
	{
		comp_num_ = new_comp_num;
	}
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of variables in the set
	 
	 \return the number of variables
	 */
	inline int num_variables() const
	{
		return num_vars_;
	}
	
	
	
	/**
	 \brief sets the total number of variables for the set
	 
	 \param new_num_vars the new total number of variables
	 */
	void set_num_variables(int new_num_vars)
	{
		num_vars_ = new_num_vars;
	}
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of natural variables in the set (those in the first group, including the homogenizing variable if present).
	 
	 \return the number of natural variables (those in the first variable group.
	 */
	inline int num_natural_variables() const
	{
		return num_natty_vars_;
	}
	
	/**
	 \brief set the number of natural variables.  
	 
	 \param new_num_nat_vars The new number of natural variables
	 */
	void set_num_natural_variables(int new_num_nat_vars)
	{
		num_natty_vars_ = new_num_nat_vars;
	}
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the incidence number for the witness set.  by default, is negative value.
	 
	 \return the incidence number, which is used for reading membership from bertini membership testing.
	 */
	inline int incidence_number() const
	{
		return incid_num_;
	}
	
	
	
	/**
	 \brief set the incidence number, probably after having determined it somehow.
	 
	 \param new_incidence the new incidence number to assign to the witness set.
	 */
	void set_incidence_number(int new_incidence)
	{
		incid_num_ = new_incidence;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// overloaded operators
	
	
	// default constructor
	
	witness_set(int nvar)
	{
		init();
		num_vars_ = num_natty_vars_ = nvar;
	};
	
	witness_set(){
		init();
	}
	
	
	~witness_set(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	/**
	 \brief custom assignment call.
	 
	 \return the new witness set
	 \param other the other witness set to copy from.
	 */
	witness_set& operator=( const witness_set& other)
	{
		reset();
		copy(other);
    return *this;
  }
	
	
	/**
	 \brief copy operator.  
	 
	 must be explicitly declared because the underlying c structures use pointers.
	 
	 \param other the other witness set to copy from.
	 */
	witness_set(const witness_set & other)
	{
		init();
		copy(other);
	}
	
	
	/**
	 \brief get her ready for action.
	 */
	void init()
	{
		input_filename_ = "unset_filename";
		
		num_vars_ = 0;
		num_natty_vars_ = 0;
		
		
		incid_num_ = -1;
		comp_num_ = dim_ = -1;
		
		
		
	}
	
	
	/** 
	 \brief perform a total deep copy of the witness set
	 
	 \param other the other witness set to copy from.
	 */
	void copy(const witness_set & other)
	{
		
		copy_skeleton(other);

		copy_names(other);
		copy_points(other);
		copy_patches(other);
		copy_linears(other);
	}
	
	
	
	/**
	 \brief copy only the witness_set data members, but not any inherited members.
	 
	 \param other the other witness set to copy from.
	 */
    void copy_skeleton(const witness_set & other)
	{
		this->input_filename_ = other.input_filename_;
		
		this->dim_ = other.dim_;
		this->comp_num_ = other.comp_num_;
		this->incid_num_ = other.incid_num_;
		
		this->num_vars_ = other.num_vars_;
		this->num_natty_vars_ = other.num_natty_vars_;
	}

    

    

	
    int num_synth_vars() const
    {
        return num_vars_ - num_natty_vars_;
    }
    
	
	/**
	 \brief get rid of any variables which are synthetic
	 */
	void only_natural_vars();
	
	/**
	 \brief trim off all but a number of variables.  trims linears, patches, and points.
	 
	 \param num_vars the number of variables to keep
	 */
	void only_first_vars(int num_vars);
	
	/**
	 \brief keep only the real points (upon dehomogenization)
	 
	 \param T the current state of the tracker, for the real tolerance.
	 */
	void sort_for_real(tracker_config_t * T);
	
	/**
	 \brief keep only the uniqie points (upon dehomogenization)
	 
	 \param T the current state of the tracker, for the unique tolerance.
	 */
	void sort_for_unique(tracker_config_t * T);
	
	
	
	/**
	 \brief keep only the points which are inside the sphere (upon dehomogenization)
	 
	 \param radius The radius of the sphere
	 \param center The center of the sphere.
	 */
	void sort_for_inside_sphere(comp_mp radius, vec_mp center);
	
	
	/**
	 \brief read in the witness set from a file, which MUST be formatted correctly.  for details on the format, see this code, or \see print_to_file()
	 
	 \param witness_set_file the path of the file to parse into this object
	 \param num_vars the number of variables in the witness set.  sadly, you have to set this manually at the time, as the header does not contain the information.  the is due to Bertini reasons.
	 */
	int  witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars);
	
	
	/**
	 \brief empty the set.  calls clear()
	 */
	void reset()
	{
		clear();
	}
	
	
	
	
	/**
	 \brief clear the entire contents of the witness set
	 */
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
	 write all the points to a file, first dehomogenizing.
	 
	 \param filename The name of the file to which to write.
	 \param indices Set of indices of points to write to file.
	 */
	void write_dehomogenized_coordinates(boost::filesystem::path filename,std::set<unsigned int> indices) const;

	
	
	
	/**
	 individual send, relative to MPI_COMM_WORLD
	 
	 \param target the ID target of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void send(parallelism_config & mpi_config, int target) const;
	
	/**
	 individual receive, relative to MPI_COMM_WORLD
	 
	 \param source the ID source of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void receive(int source, parallelism_config & mpi_config);
};






















// CURVE CELL DECOMP DATA TYPES


/**
 \brief 0-cell 
 
 a bertini_real vertex, a 0-cell.  contains a point, its projection values, its type, whether its been removed, and an index into a set of filenames contained in the vertex_set.
 
 \todo remove the metadata from this, and instead track it in the vertex set, much like the solver_output
 */
class vertex
{

private:
	point_mp pt_mp_; ///< the main data for this class -- the point.
	
	
	vec_mp  projection_values_; ///< a vector containing the projection values.
	
	int type_;  ///< See enum.
	bool removed_; ///< boolean integer whether the vertex has been 'removed' by a merge process.
	int input_filename_index_; ///< index into the vertex_set's vector of filenames.
	
public:
	
	
	/**
	 \brief get the index of the originating file name
	 
	 \return the integer index
	 */
	inline int input_filename_index() const
	{
		return input_filename_index_;
	}
	
	
	/**
	 \brief set the index
	 
	 \param new_index the new index to set in the vertex
	 */
	void set_input_filename_index(int new_index)
	{
		input_filename_index_ = new_index;
	}
	
	
	/**
	 \brief get the projection values.
	 \return the projection values in vec_mp form
	 */
	inline vec_mp* projection_values()
	{
		return &projection_values_;
	}
	
	
	const vec_mp* get_projection_values() const
	{
		return &projection_values_;
	}
	
	
	
	/**
	 \brief set the type of the vertex
	 
	 \param new_type the new type for the vertex
	 */
	void set_type(int new_type)
	{
		type_ = new_type;
	}
	
	
	/**
	 \brief get the type of the vertex
	 
	 \return the type, in integer form
	 */
	int type() const
	{
		return type_;
	}
	
	
	
	/**
	 \brief set the vertex to be 'removed'
	 
	 \param new_val the new value for the flag
	 */
	void set_removed(bool new_val)
	{
		removed_ = new_val;
	}
	
	
	/**
	 \brief query whether the vertex has been set to 'removed'
	 
	 \return whether it has been removed.
	 */
	bool is_removed() const
	{
		return removed_;
	}
	
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
		print_point_to_screen_matlab(pt_mp_,"point");
		print_point_to_screen_matlab(projection_values_,"projection_values");
		std::cout << "type: " << type_ << std::endl;
	}
	
	
	/**
	 \brief sets the point.
	 
	 Set the vertex's point to the input.
	 \param new_point the input point to be set.
	 */
	void set_point(const vec_mp new_point);
	
	
	const vec_mp* get_point() const
	{
		return &pt_mp_;
	}
	
	
	/**
	 \brief get the point in the vertex
	 
	 \return the point in vec_mp form
	 */
	vec_mp* point()
	{
		return &pt_mp_;
	}
	
	
	
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
		clear_vec_mp(this->pt_mp_);
		clear_vec_mp(this->projection_values_);
	}
	
	void copy(const vertex & other)
	{
		set_point(other.pt_mp_);
		
		vec_cp_mp(this->projection_values_, other.projection_values_);
		this->type_ = other.type_;
		
		this->input_filename_index_ = other.input_filename_index_;
		this->removed_ = other.removed_;
	}
	
	
	void init()
	{
		init_point_mp2(this->projection_values_,0,1024);
		init_point_mp2(this->pt_mp_,1,64);
		this->pt_mp_->size = 1;
		set_zero_mp(&pt_mp_->coord[0]);
		this->projection_values_->size = 0;
		this->type_ = UNSET;
		
		this->removed_ = false;
		this->input_filename_index_ = -1;
	}
};



/**
 \brief the main structure for storing vertices.  
 
 The vertex_set is bertini_real's main method for storing data.  We essentially construct a graph of vertices, consisting of edges and faces.
 
 there are methods in place to add vertices, and perform lookups.
 */
class vertex_set
{

protected:
	
	vec_mp *projections_; ///< a pointer array of projection vectors.
	int num_projections_; ///< the number of projections.  this should match the dimension of the object being decomposed.
	int curr_projection_; ///< the projection currently being used.
	
	int curr_input_index_; ///< the index of the current input file.
	std::vector< boost::filesystem::path > filenames_; ///< the set of filenames from which vertices arise.
	
	std::vector<vertex> vertices_;  ///< the main storage of points in the decomposition.
	size_t num_vertices_; ///< the number of vertices found so far.
	
	int num_natural_variables_;  ///< the number of natural variables appearing in the problem to solve.
	
	
	
	mpf_t abs_;
	mpf_t zerothresh_;
	comp_mp diff_;
	vec_mp checker_1_;
	vec_mp checker_2_;
	
	
public:
	
	
	
	
	/**
	 \brief get the currently active projection
	 
	 \return the index of the current projection
	 */
	inline int curr_projection() const
	{
		return curr_projection_;
	}
	
	
	
	
	/**
	 \brief get the number of natural variables
	 
	 \return the number of natural variables
	 */
	inline int num_natural_variables() const
	{
		return num_natural_variables_;
	}
	
	
	/**
	 \brief get the number of vertices.
	 
	 \return the number of vertices stored in this vertex set
	 */
	inline unsigned int num_vertices() const
	{
		return num_vertices_;
	}
	
	/**
	 \return the ith vertex
	 */
	const vertex& operator [](unsigned int index) const
	{
		return vertices_[index];
	}
	
	/**
	 \return the ith vertex
	 */
	vertex & operator [](unsigned int index)
	{
		return vertices_[index];
	}
	
	
	boost::filesystem::path filename(unsigned int index) const
	{
		if (index >= filenames_.size()) {
			throw std::out_of_range("trying to access filename out of range");
		}
		
		
		return filenames_[index];
	}
	
	
	void print_to_screen(); ///< operator for displaying information to screen
	
	
	/**
	 \brief add a new vertex to the set.
	 
	 \param new_vertex vertex to add to the set.
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
		this->num_natural_variables_ = num_vars;
		
		change_size_vec_mp(checker_1_, num_vars);
		change_size_vec_mp(checker_2_, num_vars);
		checker_1_->size = checker_2_->size = num_vars;
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
		for (std::vector<boost::filesystem::path>::iterator ii = filenames_.begin(); ii!= filenames_.end(); ++ii) {
			
			if (*ii == el_nom) {
				nom_index = counter;
				break;
			}
			counter++;
		}
		
		if (nom_index==-1) {
			filenames_.push_back(el_nom);
			nom_index = counter;
		}
		
		curr_input_index_ = nom_index;
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
            new_proj->size = this->num_natural_variables_;
            
			proj_index = add_projection(new_proj);
            
            new_proj->size = init_size;
		}
		
        curr_projection_ = proj_index;
        
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
        proj->size = this->num_natural_variables_;
        
        
        int proj_index = -1;
        
        for (int ii=0; ii<num_projections_; ii++) {
			if (isSamePoint_inhomogeneous_input(projections_[ii],proj)) {
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
		
		if (this->num_projections_==0) {
			projections_ = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->projections_ = (vec_mp *)br_realloc(this->projections_, (this->num_projections_+1) * sizeof(vec_mp));
		}
		
		
		init_vec_mp2(this->projections_[num_projections_],num_natural_variables_,DEFAULT_MAX_PREC);
		this->projections_[num_projections_]->size = num_natural_variables_;
		
		
		if (proj->size != num_natural_variables_) {
			vec_mp tempvec;
			init_vec_mp2(tempvec,num_natural_variables_,DEFAULT_MAX_PREC); tempvec->size = num_natural_variables_;
			for (int kk=0; kk<num_natural_variables_; kk++) {
				set_mp(&tempvec->coord[kk], &proj->coord[kk]);
			}
			vec_cp_mp(projections_[num_projections_], tempvec);
			clear_vec_mp(tempvec);
		}
		else
		{
			vec_cp_mp(projections_[num_projections_], proj);
		}
		

		num_projections_++;
		
		return num_projections_;
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
		for (int ii=0; ii<num_projections_; ii++) {
			clear_vec_mp(projections_[ii]);
		}
		num_projections_ = 0;
		
		filenames_.resize(0);
		
		vertices_.resize(0);
		num_vertices_ = 0;
		
		clear();
		init();
	}
protected:
	
	void init()
	{
		
		num_projections_ = 0;
		projections_ = NULL;
		curr_projection_ = -1;
		
		curr_input_index_ = -2;
		
		this->num_vertices_ = 0;
		this->num_natural_variables_ = 0;
		
		init_vec_mp(checker_1_,0);
		init_vec_mp(checker_2_,0);
		
		
		
		init_mp(this->diff_);

		mpf_init(abs_);
		mpf_init(zerothresh_);
		mpf_set_d(zerothresh_, 1e-8);
	}
	
	
	void copy(const vertex_set &other)
	{
		this->curr_projection_ = other.curr_projection_;
		if (this->num_projections_==0 && other.num_projections_>0) {
			projections_ = (vec_mp *) br_malloc(other.num_projections_*sizeof(vec_mp));
		}
		else if(this->num_projections_>0 && other.num_projections_>0) {
			projections_ = (vec_mp *) br_realloc(projections_,other.num_projections_*sizeof(vec_mp));
		}
		else if (this->num_projections_>0 && other.num_projections_==0){
			for (int ii=0; ii<this->num_projections_; ii++) {
				clear_vec_mp(projections_[ii]);
			}
			free(projections_);
		}
		
		for (int ii=0; ii<other.num_projections_; ii++) {
			if (ii>=this->num_projections_){
				init_vec_mp2(projections_[ii],1,1024); projections_[ii]->size = 1;
			}
			vec_cp_mp(projections_[ii],other.projections_[ii]);
		}
		
		this->num_projections_ = other.num_projections_;
		
		
		curr_input_index_ = other.curr_input_index_;
		filenames_ = other.filenames_;
		
		
		this->num_vertices_ = other.num_vertices_;
		this->num_natural_variables_ = other.num_natural_variables_;
		
		this->vertices_ = other.vertices_;
		
		vec_cp_mp(this->checker_1_,other.checker_1_);
		vec_cp_mp(this->checker_2_,other.checker_2_);
		
	}
	
	void clear()
	{
		clear_vec_mp(checker_1_);
		clear_vec_mp(checker_2_);
		
		clear_mp(diff_);
		
		mpf_clear(abs_);
		mpf_clear(zerothresh_);
		
		for (int ii=0; ii<num_projections_; ii++) {
			clear_vec_mp(projections_[ii]);
		}
		free(projections_);
	}

};






typedef std::pair<int,int> witness_set_index;


/**
 \brief metadata for witness points, for the witness_data class.
 
 
 */
class witness_point_metadata
{

	int dimension_;
	
	int corank_, typeflag_, multiplicity_, component_number_, deflations_needed_;
	double condition_number_, smallest_nonzero_sing_value_, largest_zero_sing_value_;
	
	
	
	
public:
	
	
	/**
	 \brief get the largest zero singular value
	 \return largest zero singular value
	 */
	inline double largest_zero_sing_value() const
	{
		return largest_zero_sing_value_;
	}
	
	
	
	/**
	 \brief get the smallest nonzero singular value
	 \return smallest nonsingular value
	 */
	inline double smallest_nonzero_sing_value() const
	{
		return smallest_nonzero_sing_value_;
	}
	
	
	
	/**
	 \brief get the condition number
	 \return condition number
	 */
	inline double condition_number() const
	{
		return condition_number_;
	}
	
	

	
	/**
	 \brief get the number of deflations needed
	 \return num_deflations_needed
	 */
	inline int num_deflations_needed() const
	{
		return deflations_needed_;
	}
	
	
	
	/**
	 \brief get the multiplicity
	 \return multiplicity	 */
	inline int multiplicity() const
	{
		return multiplicity_;
	}
	
	
	
	
	/**
	 \brief get the type
	 \return type
	 */
	inline int type() const
	{
		return typeflag_;
	}
	
	/**
	 \brief get the corank
	 \return corank
	 */
	inline int corank() const
	{
		return corank_;
	}
	
	
	/**
	 \brief get the dimension
	 \return dimension
	 */
	inline int dimension() const
	{
		return dimension_;
	}
	
	
	
	/**
	 \brief get the component number
	 \return component number
	 */
	inline int component_number() const
	{
		return component_number_;
	}
	
	
	
	
	/**
	 read the metadata from the witness_data file.  call only at the appropriate point.
	 
	 \param IN a pointer to an open file from which to read.  must be set to the correct point
	 */
	void set_from_file(FILE *IN)
	{
		fscanf(IN,"%lf %d %lf %lf %d %d %d %d",
			   &condition_number_,
			   &corank_,
			   &smallest_nonzero_sing_value_,
			   &largest_zero_sing_value_,
			   &typeflag_, // 10 is nonsingular, 15 is singular
			   &multiplicity_,
			   &component_number_,
			   &deflations_needed_);
	}
	
	
	
	
	/**
	 output to a stream.  only really usable with std::cout.
	 \param os the stream to put this text on.
	 \param s the system randomizer to write.
	 */
	friend std::ostream & operator<<(std::ostream &os, witness_point_metadata & s)
	{
		os << "condition_number " << s.condition_number_ << "\n";
		os << "corank " << s.corank_ << "\n";
		os << "smallest_nonzero_sing_value " << s.smallest_nonzero_sing_value_ << "\n";
		os << "largest_zero_sing_value " << s.largest_zero_sing_value_ << "\n";
		os << "typeflag " << s.typeflag_ << "\n";
		os << "multiplicity " << s.multiplicity_ << "\n";
		os << "component_number " << s.component_number_ << "\n";
		os << "deflations_needed " << s.deflations_needed_ << std::endl;
		
		return os;
	}
	
	
	
	
	
	/**
	 \brief constructor, setting the dimension in the process.
	 constructor, setting the dimension in the process.
	 
	 \param new_dim the dimension to set to.
	 */
	witness_point_metadata(int new_dim){dimension_ = new_dim;}
	
	
	witness_point_metadata(){};
	
	witness_point_metadata(const witness_point_metadata& other){copy(other);}
	
	witness_point_metadata& operator=( const witness_point_metadata& other)
	{
		copy(other);
		return *this;
	}
	
	void copy(const witness_point_metadata& other)
	{
		dimension_ = other.dimension_;
		corank_ = other.corank_;
		typeflag_ = other.typeflag_;
		multiplicity_ = other.multiplicity_;
		component_number_ = other.component_number_;
		deflations_needed_ = other.deflations_needed_;
		condition_number_ = other.condition_number_;
		smallest_nonzero_sing_value_ = other.smallest_nonzero_sing_value_;
		largest_zero_sing_value_ = other.largest_zero_sing_value_;
	}
	
	
};


/**
 \brief metadata for linears read in from the witness_data file.
 
 metadata for linears read in from the witness_data file.
 */
class witness_linear_metadata
{
	int dim_;
	
public:
	
	
	/**
	 get the dimension
	 \return the dimension
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	
	witness_linear_metadata(int new_dim){dim_ = new_dim;}
	
	witness_linear_metadata(){};
	witness_linear_metadata(const witness_linear_metadata& other)
	{
		dim_ = other.dim_;
	}
	
	witness_linear_metadata& operator=(const witness_linear_metadata &other)
	{
		dim_ = other.dim_;
		return *this;
	}
};









/**
 \brief metadata for patches read in from witness_data
 
 metadata for patches read in from witness_data
 */
class witness_patch_metadata
{
	int dim_;
	
public:
	
	/**
	 get the dimension
	 \return the dimension
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	
	witness_patch_metadata(){};
	
	
	/**
	 constructor, setting the dimension in the process
	 \param new_dim the dimension to set.
	 */
	witness_patch_metadata(int new_dim){dim_ = new_dim;}

	
	witness_patch_metadata(const witness_patch_metadata& other)
	{
		dim_ = other.dim_;
	}
	
	witness_patch_metadata& operator=(const witness_patch_metadata &other)
	{
		dim_ = other.dim_;
		return *this;
	}
	
	
	
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
	std::map<int,std::map<int,std::vector<int> > > index_tracker;
	
	std::vector< std::vector< int> > homogenization_matrix_;
	int num_variables_;
	
public:
	
	
	/**
	 \brief get the number of variables
	 
	 \return the number of variables
	 */
	inline int num_variables() const
	{
		return num_variables_;
	}
	
	
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
	
	witness_set form_specific_witness_set(int dim, int comp)	;
	
	
	
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
	void print() const
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
		
		for (unsigned int ii=0; ii<num_pts_; ii++) {
			print_point_to_screen_matlab(pts_mp_[ii],"p");
		}
	}
	
	
	
	
	
	
	
	
	
	
	
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
	
	
	
	
	
};




/**
 \brief The base cell class.  All n-cells should inherit from this.
 
 Contains a midpoint, methods to send and receive, etc.
 */
class cell
{
	
protected:

	int midpt_; ///< index into vertex set
	
	
public:
	
	/**
	 \brief get the midpoint
	 
	 \return the index of the midpoint
	 */
	inline int midpt() const
	{
		return midpt_;
	}
	
	/**
	 \brief set the midpoint
	 
	 \param new_mid the new index of the midpoint
	 \return the new index of the midpoint
	 */
	int midpt(int new_mid)
	{
		return midpt_ = new_mid;
	}
	
	
	
	
	/**
	 \brief get a cell from an input stream
	 */
	friend std::istream & operator>>(std::istream &is,  cell & c)
	{
		is >> c.midpt_;
		return is;
	}
	
	/**
	 \brief copy to another cell
	 \param other the other cell to copy to
	 */
	inline void copy(const cell & other){
		midpt(other.midpt());
	}
	
	/**
	 \brief send cell to target in communicator
	 \param target the integer id of the target of the send
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, parallelism_config & mpi_config)
	{
		int buffer = midpt();
		MPI_Send(&buffer, 1, MPI_INT, target, CELL, mpi_config.comm());
	}
	
	
	/**
	 \brief receive cell from source in communicator
	 \param source the integer id of the source of the send
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, parallelism_config & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int buffer;
		MPI_Recv(&buffer, 1, MPI_INT, source, CELL, mpi_config.comm(), &statty_mc_gatty);
		midpt(buffer);
	}
	
	/**
	 \brief virtual read_from_stream function which enforces all cells to have this functtion.
	 \param is input stream from whom to read
	 */
	virtual void read_from_stream( std::istream &is ) = 0;
	
};



/**
 \brief 1-cell.
 
 the edge data type.  has three indices: left, right, midpt.
 */
class edge : public cell
{
	int left_;  ///< index into vertices
	int right_; ///< index into vertices
	
	std::vector< int > removed_points_;
	
	
public:
	
	
	typedef std::vector< int >::iterator removed_iterator;
	typedef std::vector< int >::const_iterator removed_const_iterator;
	
	
	removed_iterator removed_begin(){return removed_points_.begin();}
	
	removed_const_iterator removed_begin() const {return removed_points_.begin();}
	
	removed_iterator removed_end(){return removed_points_.end();}
	
	removed_const_iterator removed_end() const {return removed_points_.end();}
	
	
	
	
	/**
	 \brief adds a point as a removed point.  tacks on to the end of the vector
	 
	 \param new_removed_point the index of the point to add
	 \return the index of the point
	 */
	int add_removed_point(int new_removed_point)
	{
		removed_points_.push_back(new_removed_point);
		return new_removed_point;
	}
	
	
	/**
	 \brief get the right point
	 
	 \return the index of the right point
	 */
	inline int right() const
	{
		return right_;
	}
	
	
	/**
	 \brief set the left point
	 \param new_right the new index of the left point
	 \return the index of the left point
	 */
	int right(int new_right)
	{
		return right_ = new_right;
	}
	
	
	
	
	
	/**
	 \brief get the left point
	 
	 \return the index of the left point
	 */
	inline int left() const
	{
		return left_;
	}
	
	
	/**
	 \brief set the left point
	 \param new_left the new index of the left point
	 \return the index of the left point
	 */
	int left(int new_left)
	{
		return left_ = new_left;
	}
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	edge() : cell()
	{
		left_ = right_ = -1;
	}
	
	
	/**
	 \brief construct edge from left mid and right indices
	 \param new_left the new left index for constructed edge
	 \param new_midpt the new mid index for constructed edge
	 \param new_right the new right index for constructed edge
	 */
	edge(int new_left, int new_midpt, int new_right)
	{
		left(new_left);
		right(new_right);
		midpt(new_midpt);
	}
	
	
	
	
	// other defaults are correct for this type.
	
	/**
	 check whether the edge is degenerate
	 \return true if the edge is degenerate, false if not.
	 */
	inline bool is_degenerate()
	{
		if ((left() == right()) || (left()==midpt()) || (right()==midpt()))
			return true;
		else
			return false;
	}
	
	
	
	/**
	 \brief send to a single target
	 \param target to whom to send this edge
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, parallelism_config & mpi_config)
	{
		int * buffer = new int[4];
		
		buffer[0] = left();
		buffer[1] = midpt();
		buffer[2] = right();
		buffer[3] = removed_points_.size();
		
		MPI_Send(buffer, 4, MPI_INT, target, EDGE, mpi_config.comm());
		
		delete [] buffer;
		
		buffer = new int[removed_points_.size()];
		for (unsigned int ii=0; ii!=removed_points_.size(); ii++) {
			buffer[ii] = removed_points_[ii];
		}
		MPI_Send(buffer, removed_points_.size(), MPI_INT, target, EDGE, mpi_config.comm());
		delete [] buffer;
		

	}
	
	
	/**
	 \brief receive from a single source
	 \param source from whom to receive this edge
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, parallelism_config & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int * buffer = new int[4];
		MPI_Recv(buffer, 4, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		
		left(buffer[0]);
		midpt(buffer[1]);
		right(buffer[2]);
		int temp_num_removed = buffer[3];
		
		
		delete [] buffer;
		
		buffer = new int[temp_num_removed];
		MPI_Recv(buffer, temp_num_removed, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		for (int ii=0; ii<temp_num_removed; ii++) {
			removed_points_.push_back(buffer[ii]);
		}
		
		delete [] buffer;
		
	}
	
	
	/**
	 \brief get edge from input stream.  this function is defunct, and needs implementation apparently.
	 \param is the stream from whom to read
	 */
	virtual void read_from_stream( std::istream &is )
	{
		
		is >> left_ >> midpt_ >> right_;
	}
	
	
	
};











/**
 \brief base decomposition class.  curves and surfaces inherit from this.
 
 The decomposition class holds the basic information for any dimensional decomposition -- a witness set which generated it, the number of variables, the dimension, which component it represents, the projections, the randomizer, the sphere, the input file name.
 
 */
class decomposition : public patch_holder
{


	
	
	
public:
	
	
    
/**
 \brief add a projection vector to the decomposition
 
 Decompositions in Bertini_real are computed with respect to linear projections, which are stored as vectors. They are repeated in each decomposition.  
 
 \param proj the projection to add.
 */
	void add_projection(vec_mp proj){
		if (this->num_curr_projections_==0) {
			pi_ = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->pi_ = (vec_mp *)br_realloc(this->pi_, (this->num_curr_projections_+1) * sizeof(vec_mp));
		}
		
		init_vec_mp2(this->pi_[num_curr_projections_],proj->size,DEFAULT_MAX_PREC);
		this->pi_[num_curr_projections_]->size = proj->size;
		
		vec_cp_mp(pi_[num_curr_projections_], proj);
		num_curr_projections_++;
		
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
	 
	 \throws invalid_argument if the input decomposition has no sphere yet
	 
	 \param other the decomposition which already holds sphere bounds.
	 */
	void copy_sphere_bounds(const decomposition & other)
	{
		if (!other.have_sphere_) {
			throw std::invalid_argument("trying to copy sphere bounds from a decomposition which does not have them set!");
		}
		
		set_mp(this->sphere_radius_, other.sphere_radius_);
		vec_cp_mp(this->sphere_center_, other.sphere_center_);
		this->have_sphere_ = true;
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
	
	
	
	/**
	 default constructor.
	 */
	decomposition(){
		init();
	}
	
	
	/**
	 default destructor.
	 */
	virtual ~decomposition()
	{
		this->clear();
	}
	
	
	/**
	 assignment
	 */
	decomposition & operator=(const decomposition& other){
		
		this->init();
		
		this->copy(other);
		
		return *this;
	}
	
	
	/**
	 copy
	 */
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
	
	
	
	
	
	
	
	
	
	
	

	/**
	 \brief get a shared pointer to the randomizer within
	 
	 \return the shared pointer to the system randomizer
	 */
	std::shared_ptr<system_randomizer> randomizer() const
	{
		return randomizer_;
	}
	
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of variables in the decomposition
	 
	 \return the number of variables
	 */
	inline int num_variables() const
	{
		return num_variables_;
	}
	
	/**
	 \brief set the number of variables
	 
	 \param new_num_variables the new number of variables
	 \return the number of variables
	 */
	int set_num_variables(int new_num_variables)
	{
		return num_variables_ = new_num_variables;
	}
	
	
	
	/**
	 \brief get the name of the bertini input file for this witness set
	 
	 \return the path of the file
	 */
	inline boost::filesystem::path input_filename() const
	{
		return input_filename_;
	}
	
	/**
	 \brief set the name of the bertini input file
	 
	 \param new_input_filename The new name of the file
	 */
	void set_input_filename(boost::filesystem::path new_input_filename)
	{
		input_filename_ = new_input_filename;
	}
	
	
	
	
	
	/**
	 \brief get the dimension of the set represented by the witness set.
	 
	 \return the integer dimension of the component.
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	/**
	 \brief set the dimension of the witness set
	 
	 \param new_dim the dimension of the set
	 */
	void set_dimension(int new_dim)
	{
		dim_ = new_dim;
	}
	
	
	
	
	
	
	/**
	 \brief get the component number of the set represented by the witness set.
	 
	 \return the index the component.
	 */
	inline int component_number() const
	{
		return comp_num_;
	}
	
	
	/**
	 \brief sets the component number
	 
	 \param new_comp_num The new component number to set.
	 */
	void set_component_number(int new_comp_num)
	{
		comp_num_ = new_comp_num;
	}
	
	
	
	
	
	
	
	
	
	
	
	


	

	
	
	/**
	 \brief get the witness set associated with the decomposition
	 
	 \return the witness set which generated the decomposition
	 */
	witness_set get_W() const
	{
		return W_;
	}
	
	
	/**
	 \brief get the number of current projections
	 \return the number of currently held projections
	 */
	inline int num_curr_projections() const
	{
		return num_curr_projections_;
	}
	
	
	/**
	 \brief get a pointer to the beginning of the array of projections.  
	 \throws out of range if set of projections is empty when this is requested
	 \return pointer to the 0th projection.  will be NULL if have no projections.
	 */
	inline vec_mp* pi() const
	{
		if (num_curr_projections_>0) {
			return &pi_[0];
		}
		else
		{
			throw std::out_of_range("trying get pointer for projections, but have no projections");
			return NULL;
		}
	}
	
	
	/**
	 \brief get a pointer to the ith projections.
	 \throws out of range if trying to get out of range projection
	 \return pointer to the ith projection.  will throw if out of range
	 */
	inline vec_mp* pi(int index) const
	{
		if (index >= num_curr_projections_) {
			throw std::out_of_range("trying to access an out of range projection in decomposition");
		}
		else
		{
			return &pi_[index];
		}
	}
	
	
	
	/**
	 \brief test for whether the sphere is set
	 */
	inline bool have_sphere() const
	{
		return have_sphere_;
	}
	
	/**
	 get pointer to the sphere radius
	 
	 \return pointer to the sphere radius
	 */
	comp_mp* sphere_radius()
	{
		return &sphere_radius_;
	}
	
	/**
	 \brief get pointer to the sphere's center
	 
	 \return pointer to the sphere's center
	 */
	vec_mp* sphere_center()
	{
		return &sphere_center_;
	}
	
	
	/**
	 \brief set the radius of the sphere
	 \param new_radius the radius of the sphere
	 */
	void set_sphere_radius(comp_mp new_radius)
	{
		set_mp(sphere_radius_, new_radius);
	}
	
	
	/**
	 \brief set the center of the sphere
	 \param new_center the center of the sphere
	 */
	void set_sphere_center(vec_mp new_center)
	{
		vec_cp_mp(sphere_center_,new_center);
	}
	
	
	
protected:
	
	vec_mp sphere_center_; ///< the center of the sphere.
	comp_mp sphere_radius_; ///< the radius of the sphere.
	bool have_sphere_; ///< indicates whether the decomposition has the radius set, or needs one still.
	
	
	
	int num_curr_projections_; ///< the number of projections stored in the decomposition.  should match the dimension when complete.
	vec_mp	*pi_; ///< the projections used to decompose.  first ones are used to decompose nested objects.
	
	
	

	
	
	
	/**
	 \brief set the witness set.
	 
	 \param new_w the new witness set to set.
	 */
	void set_W(const witness_set & new_w)
	{
		W_ = new_w;
	}
	
	
	
	witness_set W_; ///< generating witness set
	
	boost::filesystem::path input_filename_; ///< the name of the text file in which the system resides.
	
	
	int dim_; ///< the dimension of the decomposition
	int comp_num_; ///< the component number.
	
	
	//	function input_file;
	int num_variables_; ///< the number of variables in the decomposition
	

	
	std::shared_ptr<system_randomizer> randomizer_; ///< the randomizer for the decomposition.
	

	
	
	int add_vertex(vertex_set &V, vertex source_vertex);
	
	
	
	void init(){
		
		randomizer_ = std::make_shared<system_randomizer> (*(new system_randomizer()));

		input_filename_ = "unset";
		pi_ = NULL;
		
		
		num_curr_projections_ = 0;
		num_variables_ = 0;
		dim_ = -1;
		comp_num_ = -1;
		
		init_mp2(sphere_radius_,DEFAULT_MAX_PREC);
		init_vec_mp2(sphere_center_,0,1024);
		sphere_center_->size = 0;
		have_sphere_ = false;
		
		set_one_mp(sphere_radius_);
		neg_mp(sphere_radius_,sphere_radius_);
		

		


		
		
		
	}
	
	
	void copy(const decomposition & other)
	{
		
		
		patch_holder::copy(other);
		
		
		this->randomizer_ = other.randomizer_;//make_shared<system_randomizer>( *other.randomizer_ );
		
		this->W_ = other.W_;
		
		this->input_filename_ = other.input_filename_;
		
		this->num_variables_ = other.num_variables_;
		this->dim_ = other.dim_;
		this->comp_num_ = other.comp_num_;
		
		
		if (this->num_curr_projections_==0) {
			this->pi_ = (vec_mp *) br_malloc(other.num_curr_projections_ * sizeof(vec_mp));
		}
		else{
			for (int ii=0; ii<num_curr_projections_; ii++) {
				clear_vec_mp(pi_[ii]);
			}
			this->pi_ = (vec_mp *) br_realloc(this->pi_,other.num_curr_projections_ * sizeof(vec_mp));
		}
		
		this->num_curr_projections_ = other.num_curr_projections_;
		for (int ii = 0; ii<other.num_curr_projections_; ii++) {
			init_vec_mp2(this->pi_[ii],other.pi_[ii]->size,DEFAULT_MAX_PREC);
			this->pi_[ii]->size = other.pi_[ii]->size;
			vec_cp_mp(this->pi_[ii], other.pi_[ii])
		}
		
		
		
		
		
		copy_sphere_bounds(other);
		
		return;
	}
	
	void clear()
	{
		
		
		if (num_curr_projections_>0){
			for (int ii=0; ii<num_curr_projections_; ii++)
				clear_vec_mp(pi_[ii]);
			free(pi_);
		}
		num_curr_projections_ = 0;
		
	
		
		
		clear_mp(sphere_radius_);
		clear_vec_mp(sphere_center_);
		

		
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

