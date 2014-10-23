#ifndef NULLSPACE_H_
#define NULLSPACE_H_

/** 
 \file nullspace_left.hpp 
 
 \brief Methods for computing the left nullspace of a jacobian matrix, coupled with the system which generated it, and several linear projections.s
 */

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


#include "fileops.hpp"
#include "data_type.hpp"
#include "solver_multilintolin.hpp"
#include "solver_nullspace_left.hpp"




#include "derivative_systems.hpp"


/**
 \brief An odometer of odometers, so to speak.
 
 */
class double_odometer
{
	
private:
	
	
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
	
	double_odometer()
	{
		num_total_registers = num_active_registers = num_inactive_registers = 0;
	}
	
	double_odometer(int num_total_, int num_active_, int uniform_base)
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
	
	int reg_val(int reggie){
		return registers[reggie];
	}
	
	int act_reg(int reggie){
		return active_registers[reggie];
	}
	
	int inact_reg(int reggie){
		return inactive_registers[reggie];
	}
	
	int increment(){
		
		if (double_odometer::increment_registers()!=0) {
			if (double_odometer::increment_active_registers()!=0)
				return -1;
			else
				return 1;
		}
		else
			return 0;
	};
	
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
 \brief the main function for computing critical sets.
 
 \return SUCCESSFUL
 \param solve_out The output class for all solvers.
 \param W								input witness_set.
 \param randomizer	randomizes the system down to the correct number of equations.
 \param pi							the set of projections to use.
 \param ambient_dim The dimension of the object containing this critical set.
 \param target_dim The dimension of the critical object
 \param target_crit_dim The dimension of the critical object.
 \param program_options				holds the configuration for the main program.  is a pointer so that it is mutable.
 \param solve_options					holds the configuration for any solvers called.  is a pointer so that it is mutable.
 \param ns_config	nullspace_config object.  this is populated in this method.  must be empty as input.
 */
int compute_crit_nullspace(solver_output & solve_out, // the returned value
						   const witness_set & W,
						   std::shared_ptr<system_randomizer> randomizer,
						   vec_mp *pi,
						   int ambient_dim,
						   int target_dim, // this should also be the number of vectors in the *pi entry
						   int target_crit_dim,
						   BR_configuration & program_options,
						   solver_configuration & solve_options,
						   nullspace_config *ns_config);





/**
 \brief Put a few finishing details on the output from compute_crit_nullspace
 
 \param solve_out Returning object from solver.
 \param W input witness set
 \param ns_config The nullspace solver config object.
 */
void ns_concluding_modifications(solver_output & solve_out,
								 const witness_set & W,
								 nullspace_config * ns_config);






/**
 \brief performs the setup for the nullspace_config which is used in the compute_crit_nullspace method, and is passed into the solverNullspace.
 
 \param ns_config						the data structure we are setting up.
 \param pi the projections to use.  there could be more than used.
 \param ambient_dim the dimension of the containing object.
 \param target_dim					the dimension of the object we are detecting.
 \param target_crit_codim COdimension of the critical set.
 \param max_degree  computed -- the highest degree of any derivative of the system passed in.
 \param randomizer		how the system is randomized to the correct number of equations.
 \param W										the input witness_set
 \param solve_options The current solver setup.
 */
void nullspace_config_setup(nullspace_config *ns_config,
							vec_mp *pi, // an array of projections, the number of which is the target dimensions
							int ambient_dim,
							int target_dim,
							int target_crit_codim,
							int *max_degree, // a pointer to the value
							std::shared_ptr<system_randomizer> randomizer,
							const witness_set & W,
							solver_configuration & solve_options);




/**
 \brief Create a Bertini input file with the nullspace system.
 
 \param output_name The desired name of the output file.
 \param input_name The name of the file from which to construct the new file.
 \param program_options The current state of Bertini_real
 \param ns_config nullspace config object.
 */
void create_nullspace_system(boost::filesystem::path output_name,
							 boost::filesystem::path input_name,
							 BR_configuration & program_options,
							 nullspace_config *ns_config);


/**
 /brief Create a matlab file which will create left nullspace equations, and write it to a text file.
 
 
 \param output_name the desired output file's name
 \param input_name bertini input file out of which to create the new file.
 \param ns_config the nullspace configuration.
 \param numVars the number of variables
 \param vars The names of the variables
 \param lineVars the lines on which the variables appear
 \param numConstants the number of constants
 \param consts The names of the constants
 \param lineConstants the lines on which the constants appear
 \param numFuncs the number of functions
 \param funcs The names of the functions
 \param lineFuncs the lines on which the functions appear
 */
void createMatlabDerivative(boost::filesystem::path output_name,
							boost::filesystem::path input_name,
							nullspace_config *ns_config,
							int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs);


/**
 \brief Create a matlab file which will take the determinant of the jacobian matrix, and write it to a text file.
 
 
 \param output_name the desired output file's name
 \param input_name bertini input file out of which to create the new file.
 \param ns_config the nullspace configuration.
 \param numVars the number of variables
 \param vars The names of the variables
 \param lineVars the lines on which the variables appear
 \param numConstants the number of constants
 \param consts The names of the constants
 \param lineConstants the lines on which the constants appear
 \param numFuncs the number of functions
 \param funcs The names of the functions
 \param lineFuncs the lines on which the functions appear
 */
void create_matlab_determinantal_system(boost::filesystem::path output_name,
										boost::filesystem::path input_name,
										nullspace_config *ns_config,
										int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs);



#endif


