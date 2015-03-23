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





class ParallelismConfig; // a forward declaration
class BertiniRealConfig;   // a forward declaration

/*** low-level data types. ***/














/**
 a templated function for looking up a value in a map by key, without accidentally creating a key-value pair when it didn't previously exist, and when it doesn't exist, gives back a default value.
 
 \return the value of the key, if it exists.  if nexists, returns default_value (which was input)
 \param mc_mapperson the map to look into.
 \param lookup_key the key to search for in the map.
 \param default_value If the map doesn't hold an entry for the lookup_key, returns this value.
 */
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
 \brief An odometer of odometers, so to speak.
 
 */
class DoubleOdometer
{
	
private:
	
	/**
	 roll one mile.
	 \return 1 if rolled over, 0 if not.  if 1, then need to increment the functions.
	 */
	int increment_registers(){
		
		int carry = 1; // seed carry so that it causes addition of at least the last entry of the odometer
		for (int ii=num_active_registers-1; ii>=0; ii--) { // count down from the end of the indexes
			
			if (carry==1)
				registers[ii]++;
			
			if ( registers[ii]>=(bases[active_registers[ii]]) ) {
				registers[ii] = 0;
				carry = 1;
			}
			else{
				carry = 0;
				break;
			}
		}
		return carry;  // if return 1, then need to increment the functions.
		
	};
	
	
	int increment_active_registers(){
		int carry = 1; // seed 1
		
		for (int ii=num_active_registers-1; ii>=0; ii--) {
			
			if (carry==1){
				
				active_registers[ii]++;
				for (int jj=ii+1; jj<num_active_registers; jj++) {
					active_registers[jj] = active_registers[jj-1]+1;
				}
			}
			
			if (active_registers[num_active_registers-1]>=num_total_registers) {
				carry = 1;
			}
			else{
				carry = 0;
				break; // all done!
			}
			
		}
		
		int local_counter = 0;
		for (int ii=0; (ii<num_total_registers) && (local_counter<num_inactive_registers) ; ii++) {
			if (std::find(active_registers.begin(), active_registers.end(), ii)==active_registers.end())
			{
				inactive_registers[local_counter] = ii;
				local_counter++;
			}
		}
		
		return carry;
	};
	
public:
	
	int num_total_registers;
	int num_active_registers;
	int num_inactive_registers;
	// create and seed the function indices -- keep track of which functions we are working on
	
	std::vector< int > inactive_registers; // of length total - active
	std::vector< int > active_registers; // of length num_active_registers
	
	std::vector< int > bases; // of length num_total_registers
	std::vector< int > registers;
	
	
	/**
	 constructor
	 */
	DoubleOdometer()
	{
		num_total_registers = num_active_registers = num_inactive_registers = 0;
	}
	
	
	/**
	 constructor
	 
	 \param num_total_ the number of total registers there will be.
	 \param num_active_ the number of active registers there will be total.  this must be ≤ num_total_.  If == num_total_, this is a traditional car odometer.
	 \param uniform_base the base of all registers.  if 10, this is a traditional car odometer
	 */
	DoubleOdometer(int num_total_, int num_active_, int uniform_base)
	{
		num_total_registers = num_total_;
		num_active_registers = num_active_;
		num_inactive_registers = num_total_registers - num_active_registers;
		
		for (int ii=0; ii<num_active_registers; ii++)
			active_registers.push_back(ii);
		
		for (int ii=num_active_registers; ii<num_total_registers; ii++)
			inactive_registers.push_back(ii);
		
		for (int ii=0; ii<num_active_registers; ii++)
			registers.push_back(0);
		
		for (int ii=0; ii<num_total_registers; ii++)
			bases.push_back(uniform_base);
		
	}
	
	/**
	 constructor
	 
	 \param num_total_ the number of total registers there will be.
	 \param num_active_ the number of active registers there will be total.  this must be ≤ num_total_.  If == num_total_, this is a traditional car odometer.
	 \param new_bases the base of all registers.  if all 10, this is a traditional car odometer.  The number in here must match num_total_.
	 */
	DoubleOdometer(int num_total_, int num_active_, const std::vector<int> & new_bases)
	{
		num_total_registers = num_total_;
		num_active_registers = num_active_;
		num_inactive_registers = num_total_registers - num_active_registers;
		
		if ( num_total_!= int(new_bases.size()) ) {
			throw std::logic_error("size mismatch in creation of double odometer.  num_total must equal size of base");
		}
		
		for (int ii=0; ii<num_active_registers; ii++)
			active_registers.push_back(ii);
		
		for (int ii=num_active_registers; ii<num_total_registers; ii++)
			inactive_registers.push_back(ii);
		
		for (int ii=0; ii<num_active_registers; ii++)
			registers.push_back(0);
		
		for (int ii=0; ii<num_total_registers; ii++){
			bases.push_back(new_bases[ii]);
			std::cout << bases[ii] << std::endl;
		}
		
	}
	
	
	
	/**
	 get the register value at position
	 \param reggie the index of the register to get.
	 */
	int reg_val(int reggie){
		return registers[reggie];
	}
	
	/**
	 get the active register value at position
	 \param reggie the index of the ACTIVE register to get.
	 */
	int act_reg(int reggie){
		return active_registers[reggie];
	}
	
	
	/**
	 get the inactive register value at position
	 \param reggie the index of the INACTIVE register to get.
	 */
	int inact_reg(int reggie){
		return inactive_registers[reggie];
	}
	
	
	/**
	 a wrapper for incrementing.
	 
	 \return -1 if completely done. 0 if didn't roll over.  1 if rolled over, but not done, just moving to next active register.
	 
	 
	 */
	int increment(){
		
		if (DoubleOdometer::increment_registers()!=0) {
			if (DoubleOdometer::increment_active_registers()!=0)
				return -1;
			else
				return 1;
		}
		else
			return 0;
	};
	
	
	/**
	 print to cout
	 */
	void print(){
		std::cout << "active: ";
		for (int ii=0; ii<num_active_registers; ii++)
			std::cout << active_registers[ii] << " ";
		std::cout << "\t|\t";
		
		
		std::cout << "inactive: ";
		for (int ii=0; ii<num_inactive_registers; ii++)
			std::cout << inactive_registers[ii] << " ";
		std::cout << "\t|\t";
		
		std::cout << "register values: ";
		for (int ii=0; ii<num_active_registers; ii++)
			std::cout << registers[ii] << " ";
		std::cout << "\n";
	}
};







/** 
 \brief a woefully incomplete class to contain systems which bertini will parse.
 
 This class is intended to hold what Bertini would need to produce a straight-line program for evaluation.
 */
class Function
{
	std::string func;  ///< symbolic representation of function (straight from input file).
										 // this class is woefully incomplete.
};



/**
 \brief Comprehensive system randomization, based on deg.out.
 
 The SystemRandomizer is created based on the desired size to randomize down to, and the degrees of the functions, which are contained in the 'deg.out' file, which must pre-exist.
 
 This class does not keep track of the desired mp mode, and always populates all three randomizer matrices.
 
 This class is capable of randomizing for systems with a single hom_variable_group (set hom_var = 1), and a single variable_group.
 
 */
class SystemRandomizer
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
	 not all functions are of the same degree, and when randomizing, if there is a deficiency of a function with respect to another, you must homogenize.  this helps give that info.
	 
	 \param randomized_func The index of the output function.
	 \param base_index The index of the input function.
	 \return the level of deficiency of randomized_func with respect to base_index.
	 */
	int deficiency(unsigned int randomized_func, unsigned int base_index)
	{
		return structure_matrix[randomized_func][base_index];
	}
	
	
	/**
	 get a reference to the degrees of the output functions for this randomizer.
	 
	 \return the degrees in vector form.
	 */
	std::vector<int> & rand_degrees()
	{
		return randomized_degrees;
	}
	
	
	/**
	 get the degree of a particular randomized function.
	 
	 \throws out_of_range if there isn't any such function.
	 \param index the index of the output function.
	 \return the degree of the output function
	 */
	int randomized_degree(unsigned int index)
	{
		if (index>= randomized_degrees.size()) {
			throw std::out_of_range("trying to access a randomized degree out of range");
		}
		
		return randomized_degrees[index];
	}
	
	

	
	/**
	 output to a stream.  only really usable with std::cout, as it calls print_matrix_to_screen_matlab
	 
	 \param os the stream to put this text on.
	 \param s the system randomizer to write.
	 */
	friend std::ostream & operator<<(std::ostream &os, SystemRandomizer & s)
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
	
	
	/**
	 constructor
	 */
	SystemRandomizer()
	{
		init();
	}
	
	
	/**
	 assignment operator
	 \param other another SystemRandomizer to copy from.
	 */
	SystemRandomizer & operator=( const SystemRandomizer & other)
	{
		copy(other);
		return *this;
	}
	
	/**
	 copy constructor
	 */
	SystemRandomizer(const SystemRandomizer & other)
	{
		init();
		copy(other);
	} // re: copy
	
	
	/**
	 destructor
	 
	 because the SystemRandomizer uses Bertini data types, there are lots of pointers, and this must be done explicitly.
	 */
	~SystemRandomizer()
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
	 \brief indicates whether the SystemRandomizer is ready to go.
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
			throw std::out_of_range("trying to access an original degree out of range");
		}
		
		return original_degrees[loc];
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
	 \brief changes the precision of the SystemRandomizer
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
	 
	 send SystemRandomizer to a single target.
	 
	 \param target the id of the target.  
	 \param mpi_config container holding the mpi_config for the caller.
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	
	
	/**
	 \brief single-source receive
	 
	 receive SystemRandomizer from a single source.
	 
	 \param source the id of the source.
	 \param mpi_config container holding the mpi_config for the caller.
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
	/**
	 \brief collective MPI_COMM_WORLD broadcast send
	 
	 \see VertexSet::bcast_receive
	 
	 Send the SystemRandomizer to everyone in MPI_COMM_WORLD.
	 
	 \param mpi_config The current state of MPI, as represented in Bertini_real
	 */
	void bcast_send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief Collective MPI_COMM_WORLD broadcast receive
	 
	 Receive the SystemRandomizer from someone in MPI_COMM_WORLD.
	 
	 \param mpi_config The current state of MPI, as represented in Bertini_real
	 */
	void bcast_receive(ParallelismConfig & mpi_config);
	
	
protected:
	
	
	/**
	 because the SystemRandomizer uses Bertini data types, there are lots of pointers, and this must be done explicitly.
	 */
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
	
	
	
	/**
	 because the SystemRandomizer uses Bertini data types, there are lots of pointers, and this must be done explicitly.
	 */
	void copy(const SystemRandomizer & other)
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
class PointHolder
{
	

protected:
	
	vec_mp *pts_mp_; ///< an array of vec_mp, which are structs and require manual initialization and clearing.
	size_t num_pts_; ///< the number of stored points.
	
public:
	
	
	/**
	 get the ith point
	 
	 \return a reference to the point at the requested index.
	 \param index the index of the point to retrieve.
	 \throws out_of_range if the requested point does not exist.
	 */
	vec_mp & point(unsigned int index) const
	{
		if (index < num_pts_) {
			return pts_mp_[index];
		}
		else
		{
			throw std::out_of_range("trying to get a point in PointHolder which is out of range");
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
	 copy all the points from one PointHolder to this one.
	 
	 \param other an input PointHolder from which to copy all the points.
	 */
	void copy_points(const PointHolder & other)
	{
		
        for (unsigned int ii=0; ii<other.num_pts_; ii++)
			add_point(other.pts_mp_[ii]);
    }
	
	
	/**
	 indicate whether the PointHolder has at zero points in it.
	 
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
	 indicate whether the PointHolder has at least one point in it.
	 
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
	 \brief add a point to the PointHolder
	 
	 \return the index of the new point.
	 \param new_point the point to add.
	 */
	int add_point(vec_mp new_point);
	
	
	/**
	 constructor
	 */
	PointHolder(){
		init();
	}
	
	
	/**
	 destructor
	 */
	~PointHolder(){ // the
		
		clear();
		
	}
	
	/**
	 assignment
	 
	 \param other the other to copy from.
	 */
	PointHolder& operator=( const PointHolder & other) {
		
		reset_points();
		
		copy(other);
		return *this;
	}
	
	/**
	 copy operator.  must be explicitly declared because the underlying c structures use pointers.
	 
	 \param other the other to copy from
	 */
	PointHolder(const PointHolder & other){
		init();
		copy(other);
	}
	
	
	/**
	 copy from one to another.
	 
	 \param other the other to copy from
	 */
	void copy(const PointHolder & other)
	{
		copy_points(other);
	}
	
	
	
private:
	
	
	/**
	 set up the pointers to good value.
	 */
	void init()
	{
		pts_mp_ = NULL; // potential data loss if used improperly.
		num_pts_ = 0;
	}
	
	
	/**
	 purge the held data.
	 */
	void clear(){
		reset_points();
	}
};




/**
 \brief a way to hold vec_mp's as patches.
 
 This class offers a way to hold a set of vec_mp's as patch_mp members in a class with automated initting and clearing.
 */
class PatchHolder
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
	inline vec_mp & patch(unsigned int index) const
	{
		if (index >= num_patches_) {
			throw std::out_of_range("trying to access an out-of-range patch");
		}
		else
			return patch_mp_[index];
	}
	
	
	/**
	 \brief copy all the patches from another PatchHolder
	 
	 Copy all the stored mp patches from another PatchHolder object, without testing for uniqueness.
	 
	 \param other the PatchHolder from which to copy.
	 */
	void copy_patches(const PatchHolder & other) {
		
        for (unsigned int ii=0; ii<other.num_patches_; ii++)
			add_patch(other.patch_mp_[ii]);
    }
	
	
	
	/**
	 \brief resets the patch holder to empty.
	 
	 Reset the PatchHolder to 0 patches.
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
	
	
	/**
	 constructor
	 */
	PatchHolder(){
		init();
	}
	
	
	
	/**
	 destructor
	 */
	~PatchHolder(){ // the destructor
		
		clear();
		
	}
	
	
	/**
	 assignment
	 */
	PatchHolder& operator=( const PatchHolder& other) {
		
		reset_patches();
		
		copy(other);
		return *this;
	}
	
	/**
	 copy operator.  must be explicitly declared because the underlying c structures use pointers.
	*/
	 PatchHolder(const PatchHolder & other){
		init();
		copy(other);
	}
	
	
	/**
	 copy from another
	 
	 \param other the other to copy from.
	 */
	void copy(const PatchHolder & other)
	{
		copy_patches(other);
	}
	
	
	
private:
	
	
	/**
	 initialize
	 */
	void init()
	{
		this->patch_mp_ = NULL;
		this->num_patches_ = 0;
	}
	
	
	/**
	 purge old patches
	 */
	void clear(){
		reset_patches();
	}
};


/**
 \brief class for holding vec_mp's as linears in an automated object.
 
 This class offers automated collection of vec_mp's as L_mp.  This is necessary because the vec_mp type must be initted and cleared manually.
 */
class LinearHolder
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
	inline vec_mp & linear(unsigned int index) const
	{
		if (index>= num_linears_) {
			throw std::out_of_range("trying to access out-of-range linear");
		}
		else
			return L_mp_[index];
	}
	
	
	/**
	 \brief copies all linears from another LinearHolder.
	 
	 \param other the LinearHolder from which to copy all the linears.
	 */
	void copy_linears(const LinearHolder & other) {
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
	
	
	
	LinearHolder(){
		init();
	}
	
	
	~LinearHolder(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	LinearHolder& operator=( const LinearHolder& other) {
		reset_linears();
		copy(other);
		return *this;
	}
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	LinearHolder(const LinearHolder & other){
		init();
		copy(other);
	}
	
	
	
	void copy(const LinearHolder & other)
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
class NameHolder
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
	 \param index The index of the name desired.
	 \return the ith name.
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
	
	void copy_names(const NameHolder & nomnom)
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
class WitnessSet : public PatchHolder, public LinearHolder, public PointHolder, public NameHolder
{
	
protected:
	
	//begin data members
	

	int dim_;
	int comp_num_;
	int incid_num_;
	
	int num_vars_;
	int num_natty_vars_;
	


	boost::filesystem::path input_filename_;
	Function input_file_;
	
	
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
	
	WitnessSet(int nvar)
	{
		init();
		num_vars_ = num_natty_vars_ = nvar;
	};
	
	WitnessSet(){
		init();
	}
	
	
	~WitnessSet(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	/**
	 \brief custom assignment call.
	 
	 \return the new witness set
	 \param other the other witness set to copy from.
	 */
	WitnessSet& operator=( const WitnessSet& other)
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
	WitnessSet(const WitnessSet & other)
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
	void copy(const WitnessSet & other)
	{
		
		copy_skeleton(other);

		copy_names(other);
		copy_points(other);
		copy_patches(other);
		copy_linears(other);
	}
	
	
	
	/**
	 \brief copy only the WitnessSet data members, but not any inherited members.
	 
	 \param other the other witness set to copy from.
	 */
    void copy_skeleton(const WitnessSet & other)
	{
		this->input_filename_ = other.input_filename_;
		
		this->dim_ = other.dim_;
		this->comp_num_ = other.comp_num_;
		this->incid_num_ = other.incid_num_;
		
		this->num_vars_ = other.num_vars_;
		this->num_natty_vars_ = other.num_natty_vars_;
	}

    

    

	/**
	 get the number of synthetic variables
	 
	 \return the number of synthetic variables held.
	 */
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
	 
	 \return 0, for no good reason.
	 \param witness_set_file the path of the file to parse into this object
	 \param num_vars the number of variables in the witness set.  sadly, you have to set this manually at the time, as the header does not contain the information.  the is due to Bertini reasons.
	 */
	int  Parse(const boost::filesystem::path witness_set_file, const int num_vars);
	
	
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
	 \param T Bertini's tracker_config_t object, with settings for telling whether two points are the same.
	 */
	void merge(const WitnessSet & W_in, tracker_config_t * T);
	
	

	
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
    void send(ParallelismConfig & mpi_config, int target) const;
	
	/**
	 individual receive, relative to MPI_COMM_WORLD
	 
	 \param source the ID source of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void receive(int source, ParallelismConfig & mpi_config);
};






















// CURVE CELL DECOMP DATA TYPES


/**
 \brief 0-cell 
 
 a bertini_real Vertex, a 0-cell.  contains a point, its projection values, its type, whether its been removed, and an index into a set of filenames contained in the VertexSet.
 
 \todo remove the metadata from this, and instead track it in the Vertex set, much like the solver_output
 */
class Vertex
{

private:
	point_mp pt_mp_; ///< the main data for this class -- the point.
	
	
	vec_mp  projection_values_; ///< a vector containing the projection values.
	
	int type_;  ///< See enum.
	bool removed_; ///< boolean integer whether the Vertex has been 'removed' by a merge process.
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
	 
	 \param new_index the new index to set in the Vertex
	 */
	void set_input_filename_index(int new_index)
	{
		input_filename_index_ = new_index;
	}
	
	
	/**
	 \brief get the projection values.
	 \return the projection values in vec_mp form
	 */
	inline vec_mp& projection_values()
	{
		return projection_values_;
	}
	
	
	const vec_mp& get_projection_values() const
	{
		return projection_values_;
	}
	
	
	
	/**
	 \brief set the type of the Vertex
	 
	 \param new_type the new type for the Vertex
	 */
	void set_type(int new_type)
	{
		type_ = new_type;
	}
	
	
	/**
	 \brief get the type of the Vertex
	 
	 \return the type, in integer form
	 */
	int type() const
	{
		return type_;
	}
	
	
	
	/**
	 \brief set the Vertex to be 'removed'
	 
	 \param new_val the new value for the flag
	 */
	void set_removed(bool new_val)
	{
		removed_ = new_val;
	}
	
	
	/**
	 \brief query whether the Vertex has been set to 'removed'
	 
	 \return whether it has been removed.
	 */
	bool is_removed() const
	{
		return removed_;
	}
	
	Vertex()
	{
		init();
	}
	
	~Vertex()
	{
		clear();
		
	}
	
	Vertex & operator=(const Vertex & other)
	{
		copy(other);
		return *this;
	}
	
	Vertex(const Vertex& other)
	{
		init();
		copy(other);
	}
	
	/**
	 /brief prints the Vertex to the screen
	 Prints the Vertex to the screen
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
	
	
	const vec_mp & get_point() const
	{
		return pt_mp_;
	}
	
	
	/**
	 \brief get the point in the Vertex
	 
	 \return the point in vec_mp form
	 */
	vec_mp& point()
	{
		return pt_mp_;
	}
	
	
	
	/**
	 \brief single target mpi send.
	 
	 Send the Vertex to a single target.
	 
	 \see Vertex::receive
	 
	 \param target The MPI ID of the target for this send
	 \param mpi_config current mpi settings
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	
	/**
	 \brief single source receive.
	 
	 Receive a Vertex from a single source.
	 
	 \see Vertex::send
	 
	 \param source The MPI ID of hte source.
	 \param mpi_config current mpi settings
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
private:
	
	void clear()
	{
		clear_vec_mp(this->pt_mp_);
		clear_vec_mp(this->projection_values_);
	}
	
	void copy(const Vertex & other)
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
 
 The VertexSet is bertini_real's main method for storing data.  We essentially construct a graph of vertices, consisting of edges and faces.
 
 there are methods in place to add vertices, and perform lookups.
 */
class VertexSet
{

protected:
	
	vec_mp *projections_; ///< a pointer array of projection vectors.
	int num_projections_; ///< the number of projections.  this should match the dimension of the object being decomposed.
	int curr_projection_; ///< the projection currently being used.
	
	int curr_input_index_; ///< the index of the current input file.
	std::vector< boost::filesystem::path > filenames_; ///< the set of filenames from which vertices arise.
	
	std::vector<Vertex> vertices_;  ///< the main storage of points in the real numerical cellular decomposition.
	size_t num_vertices_; ///< the number of vertices found so far.
	
	int num_natural_variables_;  ///< the number of natural variables appearing in the problem to solve.
	
	
	double same_point_tolerance_;
	mpf_t abs_;
	mpf_t zerothresh_;
	comp_mp diff_;
	vec_mp checker_1_;
	vec_mp checker_2_;
	
	tracker_config_t * T_;
public:
	
	/**
	 get the tracker config struct
	 
	 \return a pointer to the tracker config.
	 */
	tracker_config_t * T() const
	{
		return T_;
	}
	
	/**
	 set the tracker config pointer
	 
	 \param new_T a pointer to the tracker config.
	 */
	void set_tracker_config(tracker_config_t * new_T)
	{
		T_ = new_T;
	}
	
	/**
	 get the tolerance for two points being the same
	 \return the L2 tolerance for whether two points are the same.
	 */
	double same_point_tolerance()
	{
		return same_point_tolerance_;
	}
	
	/**
	 \brief set the tolerance for points being the same
	 \param new_tolerance The new tolerance.  This is for telling whether two points are the same.
	 */
	void set_same_point_tolerance(double new_tolerance)
	{
		same_point_tolerance_ = new_tolerance;
	}
	
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
	 
	 \return the number of vertices stored in this Vertex set
	 */
	inline unsigned int num_vertices() const
	{
		return num_vertices_;
	}
	
	/**
	 \return the ith Vertex, or a reference to it.
	 \param index The index of the vertex to get.
	 */
	const Vertex& operator [](unsigned int index) const
	{
		if (index >= vertices_.size()) {
			throw std::out_of_range("trying to access Vertex out of range in VertexSet.");
		}
		return vertices_[index];
	}
	
	/**
	 \return the ith Vertex, or a reference to it.
	 \param index The index of the vertex to get.
	 */
	Vertex & operator [](unsigned int index)
	{
		if (index >= vertices_.size()) {
			throw std::out_of_range("trying to access Vertex out of range in VertexSet.");
		}
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
	 \brief add a new Vertex to the set.
	 
	 \param new_vertex Vertex to add to the set.
	 \return the index of the added Vertex
	 */
	int add_vertex(const Vertex & new_vertex);
	
	
	/**
	 \brief create a VertexSet from a file.
	 
	 Read in a VertexSet from a file.
	 
	 \param INfile the file to parse and store in a VertexSet
	 \return the number of vertices read in.
	 */
	int setup_vertices(boost::filesystem::path INfile);
	
	
	
	
	
	VertexSet(){
		init();
	}
	
	VertexSet(int num_vars){
		init();
		
		set_num_vars(num_vars);
	}
	

	

	
	
	VertexSet & operator=( const VertexSet& other) {
		init();
		copy(other);
		return *this;
	}
	
	VertexSet(const VertexSet &other)
	{
		init();
		copy(other);
	}
	
	~VertexSet()
	{
		VertexSet::clear();
	}
	
	
	
	/**
	 sets the number of variables for the VertexSet.
	 
	 \param num_vars the number of {\em natural} variables
	 */
	void set_num_vars(int num_vars)
	{
		this->num_natural_variables_ = num_vars;
		
		change_size_vec_mp(checker_1_, num_vars-1);
		change_size_vec_mp(checker_2_, num_vars-1);
		checker_1_->size = checker_2_->size = num_vars-1;
	}
	
	
	
	/**
	 \brief write VertexSet to a file, readable by bertini_real again.
	 
	 
	 write the VertexSet to a file.
	 \see setup_vertices
	 
	 Output Vertex structure as follows:
	 
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
	 
	 
	 \param outputfile the name of the file to write the VertexSet to.
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
			throw std::logic_error("trying to set curr_input from unset_filename");
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
	 
	 This function computes \f$\pi(x)\f$ for each of the points \f$x\f$ in WitnessSet W.  Then it sorts them, and computes averages.
	 
	 The output is stored in crit_downstairs and midpoints_downstairs, both pre-initialized vec_mp's.  This function also produces a std::vector<int> named index_tracker which contains the sorting of W according to \f$\pi(W)\f$.
	 
     \param W witness				set containing points of which we wish to retrieve projections values.
     \param crit_downstairs			the projection values of the input witness set, sorted for uniqueness and increasingness.
     \param midpoints_downstairs	the bisection of each interval in crit_downstairs.
     \param index_tracker			the indices of the points in W.
     \param pi						the projection we are retrieving projection values with respect to.
	 \param T pointer to a Bertini tracker_config_t object, holding the necessary settings for this method.
     \return the integer SUCCESSFUL.
     */
    int compute_downstairs_crit_midpts(const WitnessSet & W,
                                       vec_mp crit_downstairs,
                                       vec_mp midpoints_downstairs,
                                       std::vector< int > & index_tracker,
									   vec_mp pi,
									   tracker_config_t * T);
	
	
	/**
	 
	 sets the value of the [current] projection for each Vertex which has index in the set of relevant indices.
	 
	 \param relevant_indices set of Vertex indices
	 \param new_value the new value you want to set the projection value to
	 \return a vector of indices for which this operation failed, because the new and old values were too far away from each other.
	 */
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value);
	
	/**
	 
	 sets the value of the [proj_index] projection for each Vertex which has index in the set of relevant indices.
	 
	 \param relevant_indices set of Vertex indices
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
	 
	 Want to find out the index of a projection in this VertexSet? pass it into this function.
	 
	 \param proj the projection to query
	 \return the index of the projection, or -1 if it doesn't exist.
	 */
    int get_proj_index(vec_mp proj) const
    {
        int init_size = proj->size;
        proj->size = this->num_natural_variables_;
        
        
        int proj_index = -1;
        
        for (int ii=0; ii<num_projections_; ii++) {
			if (isSamePoint_inhomogeneous_input(projections_[ii],proj,same_point_tolerance_)) {
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
		
		
		init_vec_mp2(this->projections_[num_projections_],num_natural_variables_,T_->AMP_max_prec);
		this->projections_[num_projections_]->size = num_natural_variables_;
		
		
		if (proj->size != num_natural_variables_) {
			vec_mp tempvec;
			init_vec_mp2(tempvec,num_natural_variables_,T_->AMP_max_prec); tempvec->size = num_natural_variables_;
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
	 
	 Send a VertexSet to a single target process.
	 
	 \param target who to send to
	 \param mpi_config the current mpi configuration
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	
	
	/**
	 \brief single source mpi receive
	 
	 Receive a VertexSet from a single source.
	 
	 \param source the source of the receive
	 \param mpi_config the current mpi configuration
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
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
		T_ = NULL;
		
		same_point_tolerance_ = 1e-5;
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
	
	
	void copy(const VertexSet &other)
	{
		set_tracker_config(other.T());
		
		
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






/**
 \brief metadata for witness points, for the NumericalIrreducibleDecomposition class.
 
 
 */
class WitnessPointMetadata
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
	friend std::ostream & operator<<(std::ostream &os, WitnessPointMetadata & s)
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
	WitnessPointMetadata(int new_dim){dimension_ = new_dim;}
	
	
	WitnessPointMetadata(){};
	
	WitnessPointMetadata(const WitnessPointMetadata& other){copy(other);}
	
	WitnessPointMetadata& operator=( const WitnessPointMetadata& other)
	{
		copy(other);
		return *this;
	}
	
	void copy(const WitnessPointMetadata& other)
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
class WitnessLinearMetadata
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
	
	
	WitnessLinearMetadata(int new_dim){dim_ = new_dim;}
	
	WitnessLinearMetadata(){};
	WitnessLinearMetadata(const WitnessLinearMetadata& other)
	{
		dim_ = other.dim_;
	}
	
	WitnessLinearMetadata& operator=(const WitnessLinearMetadata &other)
	{
		dim_ = other.dim_;
		return *this;
	}
};









/**
 \brief metadata for patches read in from witness_data
 
 metadata for patches read in from witness_data
 */
class WitnessPatchMetadata
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
	
	
	WitnessPatchMetadata(){};
	
	
	/**
	 constructor, setting the dimension in the process
	 \param new_dim the dimension to set.
	 */
	WitnessPatchMetadata(int new_dim){dim_ = new_dim;}

	
	WitnessPatchMetadata(const WitnessPatchMetadata& other)
	{
		dim_ = other.dim_;
	}
	
	WitnessPatchMetadata& operator=(const WitnessPatchMetadata &other)
	{
		dim_ = other.dim_;
		return *this;
	}
	
	
	
};




/**
 \brief a nearly complete class for storing bertini's witness_data file.
 
 This class reads in witness_data, and produces witness_sets based on user's choice.
 */
class NumericalIrreducibleDecomposition : public PatchHolder, public LinearHolder, public PointHolder
{
	
private:
	
	std::vector< WitnessPointMetadata > point_metadata;
	std::vector< WitnessLinearMetadata > linear_metadata;
	std::vector< WitnessPatchMetadata > patch_metadata;
	
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
	
	
	/** fills this object with the sets in witness_data.
	 \param T the current state of the tracker configuration
	 */
	void populate(tracker_config_t * T);
	
	
	
	/**
	 \brief outermost method for choosing a witness set to construct.
	 
	 This function uses information stored in BertiniRealConfig to construct a witness set.
	 
	 \param options The current state of the program.  If the user passed in a particular component or dimension, this is how it gets into this method.
	 \return the chosen witness set.  may be empty.
	 */
	WitnessSet choose(BertiniRealConfig & options);
	
	
	/**
	 form a witness set automagically, based on the user's call time options.
	 
	 \param options The current state of the program.  If the user passed in a particular component or dimension, this is how it gets into this method.
	 \return the best possible witness set, based on the user's choices at call time to Bertini_real.
	 */
	WitnessSet best_possible_automatic_set(BertiniRealConfig & options);
	
	/**
	 If there are multiple dimensions and components, then the user needs to choose which he wishes to decompose.  This method is that choice.
	 
	 \param options The current state of the program.  If the user passed in a particular component or dimension, this is how it gets into this method.
	 \return A formed witness set, from the user's choices.
	 */
	WitnessSet choose_set_interactive(BertiniRealConfig & options); // lets the user choose a set, and returns a copy of it.
	
	
	/**
	 retreive from the NID the specific witness set of particular dimension and component number.  If either doesn't exist, then the returned witness set is empty.
	 
	 \return The specific witness set desired.  May be empty if the dimension or component number correspond to something don't exist.
	 \param dim The desired dimension
	 \param comp The desired component number.  Starts at 0.
	 */
	WitnessSet form_specific_witness_set(int dim, int comp)	;
	
	
	
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
	 print the witness_data to std::cout
	 
	 \todo replace this with a friend operator << ()
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
		
	}
	
	
	
	
	
	
	
	
	
	
	
private:
	
	void add_linear_w_meta(vec_mp lin, const WitnessLinearMetadata & meta)
	{
		add_linear(lin);
		linear_metadata.push_back(meta);
	}
	
	void add_patch_w_meta(vec_mp pat, const WitnessPatchMetadata & meta)
	{
		add_patch(pat);
		patch_metadata.push_back(meta);
	}
	
	
	int add_solution(vec_mp pt, const WitnessPointMetadata & meta)
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
class Cell
{
	
protected:

	int midpt_; ///< index into Vertex set
	
	
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
	 \return reference to the input stream, so you can chain inputs together.
	 \param is the input stream to get from.
	 \param c Cell to get from stream.
	 */
	friend std::istream & operator>>(std::istream &is,  Cell & c)
	{
		is >> c.midpt_;
		return is;
	}
	
	/**
	 \brief copy to another cell
	 \param other the other cell to copy to
	 */
	inline void copy(const Cell & other){
		midpt(other.midpt());
	}
	
	/**
	 \brief send cell to target in communicator
	 \param target the integer id of the target of the send
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, ParallelismConfig & mpi_config)
	{
		int buffer = midpt();
		MPI_Send(&buffer, 1, MPI_INT, target, CELL, mpi_config.comm());
	}
	
	
	/**
	 \brief receive cell from source in communicator
	 \param source the integer id of the source of the send
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, ParallelismConfig & mpi_config)
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
 \brief base Decomposition class.  curves and surfaces inherit from this.
 
 The Decomposition class holds the basic information for any dimensional decomposition -- a witness set which generated it, the number of variables, the dimension, which component it represents, the projections, the randomizer, the sphere, the input file name.
 
 */
class Decomposition : public PatchHolder
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
		
		init_vec_mp2(this->pi_[num_curr_projections_],proj->size,proj->curr_prec);
		this->pi_[num_curr_projections_]->size = proj->size;
		
		vec_cp_mp(pi_[num_curr_projections_], proj);
		num_curr_projections_++;
		
	}
	
	
	/**
	 \brief commit a set of points to the Vertex set, associating them with the input file for this decomposition.
	 
	 \todo rename this function to something more accurately descriptive
	 
	 \return the number 0.  this seems pointless.
	 \param W the witness set containing the points to add.
	 \param add_type the type the points will inherit.  {\em e.g.} CRITICAL
	 \param V the Vertex set to add the points to.
	 */
    int add_witness_set(const WitnessSet & W, int add_type, VertexSet & V);
    
	
	/**
	 \brief Find the index of a testpoint.
	 
	 Search the Vertex set passed in for testpoint.  This happens relative to this decomposition for no reason whatsoever.
	 \todo change this to be a property of the Vertex set
	 
	 \return the index of the point, or -1 if not found.
	 \param V the Vertex set in which to search
	 \param testpoint the point for which to search
	 */
	int index_in_vertices(VertexSet &V,
                          vec_mp testpoint);
	
	
	/**
	 \brief search for a testvertex, and add it to Vertex set, and this decomposition, if not found.
	 
	 Search the passed in Vertex set for the testvertex -- and add it to the Vertex set, and its index to the decomposition if not found.
	 
	 \return the index of the point, or -1 if not found.
	 \param V the Vertex set in which to search
	 \param vert a Vertex with point for which to search.
	 */
	int index_in_vertices_with_add(VertexSet &V,
                                   Vertex vert);
	
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
	void compute_sphere_bounds(const WitnessSet & W_crit);
	
	
	/**
	 \brief copy the bounds from another Decomposition
	 
	 Sub-decompositions will often want to inherit the bounding sphere of another.  This method lets you copy from one to another.
	 
	 \throws invalid_argument if the input Decomposition has no sphere yet
	 
	 \param other the Decomposition which already holds sphere bounds.
	 */
	void copy_sphere_bounds(const Decomposition & other)
	{
		if (!other.have_sphere_) {
			throw std::invalid_argument("trying to copy sphere bounds from a Decomposition which does not have them set!");
		}
		
		set_mp(this->sphere_radius_, other.sphere_radius_);
		vec_cp_mp(this->sphere_center_, other.sphere_center_);
		this->have_sphere_ = true;
	}
	
	
	/**
	 \brief the main way to print a Decomposition to a file.
	 
	 This method backs up the existing folder to one siffixed with "_bak", and creates a new folder with the correct name, to which it prints the Decomposition in text file format.
	 
	 \param base the base folder name to print the Decomposition.
	 */
	void output_main(const boost::filesystem::path base);
	
	
	/**
	 \brief copy the component number, filename, number of variables, and witness set itself into the Decomposition.
	 
	 \param W the witness set with the data.
	 */
	void copy_data_from_witness_set(const WitnessSet & W);
	
	
	/**
	 reset Decomposition to empty.
	 */
	void reset()
	{
		clear();
		init();
	}
	
	
	
	/**
	 default constructor.
	 */
	Decomposition(){
		init();
	}
	
	
	/**
	 default destructor.
	 */
	virtual ~Decomposition()
	{
		this->clear();
	}
	
	
	/**
	 assignment
	 \return The assigned decomposition.
	 \param other The input decomposition, from which to assign.
	 */
	Decomposition & operator=(const Decomposition& other){
		
		this->init();
		
		this->copy(other);
		
		return *this;
	}
	
	
	/**
	 copy-constructor
	 \param other Another decomposition from which to copy-construct.
	 */
	Decomposition(const Decomposition & other){
		
		this->init();
		
		this->copy(other);
	}
	
	
	/**
	 \brief single target MPI send.
	 
	 Send base Decomposition to another process.
	 
	 \param target the ID of the worker to which to send.
	 \param mpi_config the current configuration of MPI
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	/**
	 \brief single source MPI receive.
	 
	 Receive base Decomposition from another process.
	 
	 \param source the ID of the process from which to receive
	 \param mpi_config the current configuration of MPI
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
	
	
	
	
	
	
	
	
	

	/**
	 \brief get a shared pointer to the randomizer within
	 
	 \return the shared pointer to the system randomizer
	 */
	std::shared_ptr<SystemRandomizer> randomizer() const
	{
		return randomizer_;
	}
	
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of variables in the Decomposition
	 
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
	 \brief get the witness set associated with the Decomposition
	 
	 \return the witness set which generated the Decomposition
	 */
	const WitnessSet& get_W() const
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
	inline vec_mp& pi() const
	{
		if (num_curr_projections_>0) {
			return pi_[0];
		}
		else
		{
			throw std::out_of_range("trying get pointer for projections, but have no projections");
		}
	}
	
	
	/**
	 \brief get a pointer to the ith projections.
	 \throws out of range if trying to get out of range projection
	 \param index The index of the projection you want to get.
	 \return pointer to the ith projection.  will throw if out of range
	 */
	inline vec_mp& pi(int index) const
	{
		if (index >= num_curr_projections_) {
			throw std::out_of_range("trying to access an out of range projection in Decomposition");
		}
		else
		{
			return pi_[index];
		}
	}
	
	
	
	/**
	 \brief test for whether the sphere is set
	 \return Whether have the sphere radius and center set, whether from file or computed from a set of points.
	 */
	inline bool have_sphere() const
	{
		return have_sphere_;
	}
	
	/**
	 get pointer to the sphere radius
	 
	 \return pointer to the sphere radius
	 */
	comp_mp& sphere_radius()
	{
		return sphere_radius_;
	}
	
	/**
	 \brief get pointer to the sphere's center
	 
	 \return pointer to the sphere's center
	 */
	vec_mp& sphere_center()
	{
		return sphere_center_;
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
	
	
	
	/**
	 make a deep copy of another Decomposition
	 
	 \param other The other Decomposition, from which to clone into this.
	 */
	void clone(const Decomposition & other)
	{
		PatchHolder::copy(other);
		
		
		this->randomizer_ = std::make_shared<SystemRandomizer>( *other.randomizer_ );
		
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
			init_vec_mp2(this->pi_[ii],other.pi_[ii]->size,other.pi_[ii]->curr_prec);
			this->pi_[ii]->size = other.pi_[ii]->size;
			vec_cp_mp(this->pi_[ii], other.pi_[ii])
		}
		
		
		
		
		
		copy_sphere_bounds(other);

	}
	
	
	
protected:
	
	vec_mp sphere_center_; ///< the center of the sphere.
	comp_mp sphere_radius_; ///< the radius of the sphere.
	bool have_sphere_; ///< indicates whether the Decomposition has the radius set, or needs one still.
	
	
	
	int num_curr_projections_; ///< the number of projections stored in the Decomposition.  should match the dimension when complete.
	vec_mp	*pi_; ///< the projections used to decompose.  first ones are used to decompose nested objects.
	
	
	

	
	
	
	/**
	 \brief set the witness set.
	 
	 \param new_w the new witness set to set.
	 */
	void set_W(const WitnessSet & new_w)
	{
		W_ = new_w;
	}
	
	
	
	WitnessSet W_; ///< generating witness set
	
	boost::filesystem::path input_filename_; ///< the name of the text file in which the system resides.
	
	
	int dim_; ///< the dimension of the Decomposition
	int comp_num_; ///< the component number.
	
	
	//	function input_file;
	int num_variables_; ///< the number of variables in the Decomposition
	
	
	std::shared_ptr<SystemRandomizer> randomizer_; ///< the randomizer for the Decomposition.
	

	
	
	int add_vertex(VertexSet &V, Vertex source_vertex);
	
	
	
	void init(){
		
		randomizer_ = std::make_shared<SystemRandomizer> (*(new SystemRandomizer()));

		input_filename_ = "unset";
		pi_ = NULL;
		
		
		num_curr_projections_ = 0;
		num_variables_ = 0;
		dim_ = -1;
		comp_num_ = -1;
		
		init_mp2(sphere_radius_,1024);
		init_vec_mp2(sphere_center_,0,1024);
		sphere_center_->size = 0;
		have_sphere_ = false;
		
		set_one_mp(sphere_radius_);
		neg_mp(sphere_radius_,sphere_radius_);
	}
	
	
	
	
	void copy(const Decomposition & other)
	{
		PatchHolder::copy(other);
		
		
		this->randomizer_ = other.randomizer_;
		
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
			init_vec_mp2(this->pi_[ii],other.pi_[ii]->size,1024);
			this->pi_[ii]->size = other.pi_[ii]->size;
			vec_cp_mp(this->pi_[ii], other.pi_[ii])
		}
		
		
		if (other.have_sphere_) {
			copy_sphere_bounds(other);
		}
		else{
			this.have_sphere_ = false;
		}
		
		
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
	
	
	
	

	
}; // end Decomposition













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

