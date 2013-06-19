
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

#include <boost/filesystem/path.hpp>
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

typedef struct
{
	
	
	
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
	boost::filesystem::path output_basename; 
	
	
	int verbose_level;
	
	int use_gamma_trick; // bool
	int use_bounding_box; // bool
	
	int crit_solver;
} program_configuration;




typedef struct
{
	int stifle_membership_screen; //< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text;
	
	int verbose_level;
	
	int maximum_num_iterations;
	
	int use_gamma_trick;
	mpf_t TOL;
} sampler_configuration;






/**
 splits the bertini input file into several files for later use.
 */
void parse_input_file(boost::filesystem::path filename, int *MPType);


void parse_preproc_data(boost::filesystem::path filename, preproc_data *PPD);



/**
 reads in projection from file if user specified, creates one otherwise.
 --
// currently defaults to create a random real projection with homogeneous value 0;
 */
void get_projection(vec_mp *pi,
										program_configuration program_options,
										solver_configuration solve_options,
										int num_vars,
										int num_projections);




/** display options to user. */
void BR_print_usage();


/** get the program_configuration from the command line. */
int BR_parse_commandline(int argC, char *args[], program_configuration *options);


/** check to make sure files are in place, etc.  */
int BR_startup(program_configuration options);


/** displays the bertini_real splash screen */
void BR_splash_screen();

/** prints the current configuration to the screen, and pauses. */
void BR_display_current_options(program_configuration options);


void BR_init_config(program_configuration *options);
void BR_clear_config(program_configuration *options);



/** get the sampler_configuration from the command line. */
int  sampler_parse_commandline(int argc, char **argv, sampler_configuration *options);
void sampler_splash_screen();
void sampler_print_usage();
int  sampler_parse_options(int argc, char **argv, sampler_configuration *options);


void sampler_init_config(sampler_configuration *options);
void sampler_clear_config(sampler_configuration *options);













#endif


