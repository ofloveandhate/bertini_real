
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


#include "missing_bertini_headers.hpp"


#include "data_type.hpp"


// enum for worker mode choice
enum {NULLSPACE = 3000, LINPRODTODETJAC, DETJACTODETJAC, LINTOLIN, MULTILIN};

enum {TERMINATE = 2000, INITIAL_STATE};

enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION};


///////////
//
//    PROGRAM CONFIGURATION
//
//////////

#define MAX_STRLEN 200

class parallelism_config
{
	
public:
	
	
	parallelism_config(){
		

		
		
		MPI_Comm_size(MPI_COMM_WORLD, &this->numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &this->my_id);
		
		
		if (my_id==0)
			worker_level = 0;
		else
			worker_level = 1;
		
		headnode = 0;
		
		
		
		
		
		
//		MPI_Group orig_group, new_group;
//		
//		/* Extract the original group handle */
//		
//		MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
//		
//		/* Create new communicator and then perform collective communications */
//		
//		MPI_Comm_create(MPI_COMM_WORLD, orig_group, &my_communicator);
//		
//		MPI_Comm_size(my_communicator, &this->numprocs);
//		MPI_Comm_rank(my_communicator, &this->my_id);
		
		
		
		
	}
	
	bool is_head()
	{
		if (my_id == headnode)
			return true;
		else
			return false;
	}
	
	
	bool use_parallel(){
		if (numprocs>1)
			return true;
		else
			return false;
	}
	
	
	
	
	int headnode;
	int my_id, my_id_global;
	int numprocs;
	MPI_Comm   my_communicator;
	
	int worker_level; // higher worker level means more tedious work, in a sense.  worker_level 0 is uber-master.  worker_level 1 will be the next level down in management, etc.  the exact usage of this is relative to the process being run.
	
	void abort(int why){
		MPI_Abort(my_communicator,why);
	}
	
	int head(){return headnode;};
	int id(){ return my_id;};
	int level(){return worker_level;};
	int size(){ return numprocs;};
	
		
};




class prog_config : public parallelism_config
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
	
	std::string bertini_command;
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
	
	
	BR_configuration() : prog_config()
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
		
		this->bertini_command = "bertini_serial";
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


























#endif


