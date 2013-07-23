
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

#include <getopt.h> 

#include <boost/filesystem.hpp>
#include <iostream>


#define BERTINI_REAL_VERSION_STRING "0.0.101"



#ifndef _PROGRAM_STARTUP_H
#define _PROGRAM_STARTUP_H


extern "C" {
#include "polysolve.h"
}
extern "C" {
#include "cascade.h"
}

#include "solver.hpp"
#include "data_type.hpp"
#include "missing_bertini_headers.hpp"

// enum for crit solver choice
enum {NULLSPACE, LINPRODTODETJAC};

///////////
//
//    PROGRAM CONFIGURATION
//
//////////

#define MAX_STRLEN 200


class prog_config
{
	
private:
	
	
	
public:
	int verbose_level;
	
	boost::filesystem::path called_dir;
	boost::filesystem::path working_dir;
	boost::filesystem::path output_dir;
	
	void move_to_temp();
	
	void move_to_called();
};




class BR_configuration : public prog_config
{
public:
	int max_deflations;
	
	int stifle_membership_screen; //< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text; // std::string
	
	int user_randomization; // bool
	int user_projection; // bool
	
	int MPType;
	
	// all these char* should be boost::filesystem::path
	boost::filesystem::path projection_filename;
	boost::filesystem::path randomization_filename;
	boost::filesystem::path input_filename;
	boost::filesystem::path witness_set_filename;
	boost::filesystem::path input_deflated_filename;
	 
	
	boost::filesystem::path current_working_filename;
	
	
	std::string matlab_command;
	int verbose_level;
	
	int use_gamma_trick; // bool
	int use_bounding_box; // bool
	
	int crit_solver;
	

	
	
	/** display options to user. */
	void print_usage();
	
	
	/** get the BR_configuration from the command line. */
	int parse_commandline(int argC, char *args[]);
	
	
	/** check to make sure files are in place, etc.  */
	int startup();
	
	
	/** displays the bertini_real splash screen */
	void splash_screen();
	
	/** prints the current configuration to the screen, and pauses. */
	void display_current_options();
	
	
	BR_configuration()
	{
		
		this->max_deflations = 10;
		
		this->user_projection = 0;
		this->projection_filename = "";
		
		this->user_randomization = 0;
		this->randomization_filename = "";
		
		this->input_filename = "input";
		this->current_working_filename = this->input_filename;
		
		this->witness_set_filename = "witness_set";
		
		this->output_dir = boost::filesystem::absolute("output");
		
		
		this->crit_solver = NULLSPACE;
		
		this->stifle_membership_screen = 1;
		this->stifle_text = " > /dev/null ";
		
		this->matlab_command = "matlab -nosplash";
		this->verbose_level = 0; // default to 0
		
		this->MPType = 2;
		
		this->use_bounding_box = 0;
		this->use_gamma_trick = 0;
		return;

	};
	
	
	
}; //re: BR_configuration




class sampler_configuration
{
public:
	int stifle_membership_screen; //< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text;
	
	int verbose_level;
	
	int maximum_num_iterations;
	
	int use_gamma_trick;
	mpf_t TOL;
	
	
	/** get the sampler_configuration from the command line. */
	int  parse_commandline(int argc, char **argv);
	void splash_screen();
	void print_usage();
	int  parse_options(int argc, char **argv);
	
	
	sampler_configuration()
	{
		this->stifle_membership_screen = 1;
		this->stifle_text = (char *)bmalloc(MAX_STRLEN*sizeof(char));
		this->stifle_text = " > /dev/null ";
		
		this->verbose_level = 0; // default to 0
		
		this->maximum_num_iterations = 10;
		
		mpf_init(this->TOL);
		mpf_set_d(this->TOL, 1e-1); // this should be made adaptive to the span of the projection values or the endpoints
		
		this->use_gamma_trick = 0;
	};
	
	~sampler_configuration()
	{
		mpf_clear(this->TOL);	
	}

	
	
	
};






/**
 splits the bertini input file into several files for later use.
 */
void parse_input_file(boost::filesystem::path filename);

void parse_input_file(boost::filesystem::path filename, int * MPType);

void parse_preproc_data(boost::filesystem::path filename, preproc_data *PPD);



/**
 reads in projection from file if user specified, creates one otherwise.
 --
// currently defaults to create a random real projection with homogeneous value 0;
 */
void get_projection(vec_mp *pi,
										BR_configuration program_options,
										solver_configuration solve_options,
										int num_vars,
										int num_projections);






















#endif


