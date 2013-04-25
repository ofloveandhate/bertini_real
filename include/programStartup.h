
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

#define BERTINI_REAL_VERSION_STRING "0.0.101"



#ifndef _PROGRAM_STARTUP_H
#define _PROGRAM_STARTUP_H

#include "polysolve.h"
#include "fileops.h"
#include "data_type.h"

/**
 displays the bertini_real splash screen
 */
void splash_screen();

/**
 prints the current configuration to the screen, and pauses.
 */
void display_current_options(program_configuration options);

/**
 get the program_configuration from the command line.
 */
int parse_options(int argC, char *args[], program_configuration *options);

/**
 get the sampler_configuration from the command line.
 */
int sampler_parse_options(int argc, char **argv, sampler_configuration *options);

/**
	check to make sure files are in place, etc.
	*/
int startup(program_configuration options);


/**
 splits the bertini input file into several files for later use.
 */
void parse_input_file(char filename[], int *MPType);

/**
 reads the tracker_config_t from file.
 */
void get_tracker_config(solver_configuration *solve_options,int MPType);


/**
 reads in projection from file if user specified, creates one otherwise.
--
 currently hardcoded to project onto the first coordinate for the decomposition.
 */
void get_projection(vec_mp pi_mp,
										program_configuration program_options,
										solver_configuration solve_options,
										int num_var);


void sampler_splash_screen();
void sampler_print_usage();
int sampler_parse_options(int argc, char **argv, sampler_configuration *options);
#endif


