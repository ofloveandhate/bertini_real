#ifndef _PROGRAM_STARTUP_H
#define _PROGRAM_STARTUP_H

/** \file programConfiguration.hpp */

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
#include <queue>
#include <map>

#include "bertini_extensions.hpp"

#define BERTINI_REAL_VERSION_STRING "0.1.0"
#define SAMPLER_VERSION_STRING "0.9.9"

enum {INACTIVE = 500, ACTIVE};
enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};

#include "fileops.hpp"




///////////
//
//    PROGRAM CONFIGURATION
//
//////////


/**
 \brief Carries the configuration of MPI, including ID number, communicator, id of headnode, etc.
 
 
 */
class parallelism_config
{
	
protected:
	
	
public:
	bool force_no_parallel;
	int headnode;
	int my_id, my_id_global;
	int numprocs;
	MPI_Comm   my_communicator;
	
	

	
	int worker_level; // higher worker level means more tedious work, in a sense.  worker_level 0 is uber-master.  worker_level 1 will be the next level down in management, etc.  the exact usage of this is relative to the process being run.
	

	
	std::map< int, int> worker_status;
	
	std::queue< int > available_workers;
	
	
	parallelism_config(){
		init();
	}
	
	
	/**
	 \brief determine whether the process is the head.
	 
	 \return Indicator of whether the process thinks it is currently the head.
	 */
	inline bool is_head()
	{
		if (my_id == headnode)
			return true;
		else
			return false;
	}
	
	
	/**
	 
	 \brief Determine whether the process thinks it should be in parallel mode.
	 
	 \return Indicator of whether to use parallel mode.
	 */
	inline bool use_parallel(){
		if (numprocs>1 && force_no_parallel!=true)
			return true;
		else
			return false;
		
		
	}
	
	
	
	
	/**
	 \brief Call MPI_Abort using the internally stored communicator.
	 
	 \param why An integer to feed to MPI_Abort.
	 */
	void abort(int why){
		MPI_Abort(my_communicator,why);
	}
	
	/**
	 Get the current communicator
	 \return the currently stored communicator
	 */
	inline MPI_Comm comm(){return my_communicator;}
	
	/**
	 \brief Get the ID of the supervisor
	 
	 \return the ID of the head node, supervisor, or whatever you want to call it.
	 */
	inline int head(){return headnode;}
	
	/**
	 \brief Get the ID of this process
	 
	 \return the ID
	 */
	inline int id(){ return my_id;}
	
	/**
	 \brief  Get the worker level
	 
	 \return the worker level
	 */
	inline int level(){return worker_level;}
	
	/**
	 \brief Get how many workers there are in the current communicator
	 
	 \return The number of workers.
	 */
	inline int size(){ return numprocs;}
	
	
	/**
	 \brief Set up the vector of available workers, based on how many processors there are.
	 
	 It also sets up the workers to be listed as inactive.
	 */
	void init_active_workers()
	{

		for (int ii=1; ii<this->numprocs; ii++) {
			available_workers.push(ii);
			worker_status[ii] = INACTIVE;
		}
		
		
	}
	
	/**
	 \brief Get the next available worker, relist it as active, and return its id.
	 
	 Available workers are stored as a queue of integers, and work is assigned to the front of the vector.  This method pops the front entry of the queue, relists it as active, and returns its ID.
	 \return the ID of the next worker.
	 */
	int activate_next_worker()
	{
		int worker_id = available_workers.front();
		available_workers.pop();
		
		if (worker_status[worker_id] == ACTIVE) {
			std::cout << "master tried making worker" << worker_id << " active when it was already active" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		worker_status[worker_id] = ACTIVE;
		

		return worker_id;
	}
	
	/**
	 \brief Take an active worker, relist it as inactive, and make it available.  If the worker is not active, it calls MPI_Abort
	 
	 \param worker_id The id of the worker to deactivate.
	 */
	void deactivate(int worker_id)
	{
		if (worker_status[worker_id] == INACTIVE) {
			std::cout << "master tried decativating worker" << worker_id << " when it was already inactive" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,2);
		}
		worker_status[worker_id] = INACTIVE;
		
		
		available_workers.push(worker_id);
	}
	
	
	/**
	 \brief Individually send a number to all available workers.
	 
	 pops available workers off the queue, and deliveres the them numtosend.
	 
	 \param numtosend The integer number to send.
	 */
	void send_all_available(int numtosend)
	{
		while (available_workers.size()>0)  {
			int sendtome = available_workers.front();
			MPI_Send(&numtosend, 1, MPI_INT, sendtome, NUMPACKETS, MPI_COMM_WORLD);
			available_workers.pop();
		}
	}
	
	
	/**
	 \brief MPI_Bcast a number to everyone in the ring, calling for help.
	 
	 \param solver_type The case-index of the type of help the head wants.
	 */
	void call_for_help(int solver_type)
	{

		MPI_Bcast(&solver_type, 1, MPI_INT, head(), MPI_COMM_WORLD);
		
		init_active_workers();
		
	}
	
	/**
	 \brief check if there are available workers.
	 
	 \return a boolean indicating if there are available workers.
	 */
	bool have_available()
	{
		if (available_workers.size()==0) {
			return false;
		}
		else
		{
			return true;
		}
		
	}
	
	/**
	 \brief check if there are workers working.
	 
	 \return A boolean indicating whether there are workers with the status 'ACTIVE'.
	 */
	bool have_active()
	{
		bool yep = false;
		for (int ii=1; ii<this->numprocs; ii++) {
			if (this->worker_status[ii]==ACTIVE) {
				yep = true;
				break;
			}
		}
		return yep;
	}
	
	/**
	 \brief Get the number of workers listed as active.
	 
	 \return the number of workers listed as Active.
	 */
	int num_active()
	{
		int num = 0;
		for (int ii=1; ii<this->numprocs; ii++) {
			if (this->worker_status[ii]==ACTIVE) {
				num++;
			}
		}
		return num;
	}
	
	
private:
	
	void init()
	{
		
		force_no_parallel = false;
		numprocs = 1;
		headnode = 0;
		
		MPI_Comm_size(MPI_COMM_WORLD, &this->numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &this->my_id);
		
		
		if (is_head())
			worker_level = 0;
		else
			worker_level = 1;
		
		
		
		my_communicator = MPI_COMM_WORLD; // default communicator is MPI_COMM_WORLD
		
        
        
        
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
		return;
	}
};



/**
 \brief Base class for program configuations.
 
 Both sampler_configuration and BR_configuration inherit from this.
 */
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



/**
 \brief holds the current state of configuration for Bertini_real.
 
 */
class BR_configuration : public prog_config
{
public:
	
	int max_deflations; ///< the maximum allowable number of deflation iterations before it gives up.
	int debugwait; ///< flag for whether to wait 30 seconds before starting, and print the master process ID to screen.
	int stifle_membership_screen; ///< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text; ///< string to append to system commands to stifle screen output.
	
	int quick_run;  ///< indicator of whether to use the robust solver wherever possible
	bool user_sphere; ///< flag for whether to read the sphere from a file, rather than compute it.

	bool user_projection; // indicator for whether to read the projection from a file, rather than randomly choose it.
	
	int MPType; ///< store M.O.
	
	boost::filesystem::path bounding_sphere_filename; ///< name of file to read if user_sphere==true
	boost::filesystem::path projection_filename; ///name of file to read if user_projection==true
	boost::filesystem::path input_filename; ///< name of the input file to read -- by default it's "input"

	boost::filesystem::path input_deflated_filename; ///< the name of the file post-deflation
	boost::filesystem::path sphere_filename; ///< the name of the sphere file.  this seems like a duplicate.
	

	
	int target_dimension;  ///< the dimension to shoot for
	int target_component;  ///< the integer index of the component to decompose.  by default, it's -2, which indicates 'ask me'.
	
	
	std::string bertini_command; ///< the string of what to call for bertini.
	std::string matlab_command; ///< the string for how to call matlab.
	
	bool use_gamma_trick; ///< indicator for whether to use the gamma trick in a particular solver.
	
	bool merge_edges; ///< a mode switch, indicates whether should be merging.
	
	

	
	
	/** 
	 \brief  display options to user. */
	void print_usage();
	
	
	/** 
	 \brief get the BR_configuration from the command line. 
	 
	 \return the number 0.
	 \param argC the command count, from main()
	 \param args the input command string from main()
	 */
	int parse_commandline(int argC, char *args[]);
	
	
	/** 
	 \brief check to make sure files are in place, etc.  
	 
	 checks for write priveledges, and for the existence of bounding_sphere_filename and projection_filename.
	 
	 
	 \return the number 0.
	 */
	int startup();
	
	
	/** 
	 \brief displays the bertini_real splash screen */
	void splash_screen();
	
	/** 
	 \brief prints the current configuration to the screen, and pauses. */
	void display_current_options();
	
	
	BR_configuration() : prog_config()
	{
		
        init();
	};
	
	void init();
	
}; //re: BR_configuration




class sampler_configuration : public prog_config
{
public:
	int stifle_membership_screen; ///< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text; ///< the text to append to system() commands to stifle screen output
	
	
	int maximum_num_iterations; ///< the maximum number of passes for iterative adaptive sampling
	
	int use_gamma_trick; ///< indicator for whether to use the gamma trick.
	mpf_t TOL; ///< the distance-tolerance for spatial-adaptive sampling
	
	bool no_duplicates; ///< a flag for whether to never duplicate points in the vertex_set as it is constructed.
	
	bool use_distance_condition; ///< switch for adaptive modes, between distance or movement breaking of while loop.
	bool use_fixed_sampler; ///< mode switch between adaptive and fixed-number.
	int target_num_samples; ///< the number of samples per cell, more or less.
	
	bool use_projection_binning; ///< switch for whether to use a projection-value based binning method for stitching ribs of unequal length, or to use a distance-based method.
	
	/** 
	 \brief get the sampler_configuration from the command line. */
	int  parse_commandline(int argc, char **argv);
	
	/**
	 \brief print a splash opening message, including the version number and author list
	 */
	void splash_screen();
	
	/**
	 \brief print a message to screen about how to use the program
	 */
	void print_usage();
	
	/**
	 \brief parse the command line for options, using getopt_long_*
	 */
	int  parse_options(int argc, char **argv);
	
	
	/**
	 \brief default constructor, contains some default settings.
	 */
	sampler_configuration()
	{
		no_duplicates = true;
		use_fixed_sampler = false;
		use_distance_condition = false;
		
		target_num_samples = 10;
		
		use_projection_binning = false;
		
		stifle_membership_screen = 1;
		stifle_text = " > /dev/null ";
		
		verbose_level = 0; // default to 0
		
		maximum_num_iterations = 10;
		
		mpf_init(TOL);
		mpf_set_d(TOL, 1e-1); // this should be made adaptive to the span of the projection values or the endpoints
		
		use_gamma_trick = 0;
	};
	
	~sampler_configuration()
	{
		mpf_clear(TOL);
	}

	
	
	
};






/**
 \brief splits the bertini input file into several files for later use.
 
 calls MPI_Bcast(&PARSING, 1, MPI_INT, 0, MPI_COMM_WORLD); to be able to let workers carry through.  sadly, the parse_input() method in Bertini calls an MPI_Bcast, which will trip up the workers if it is not caught.
 
 basically, the files made are
 • num.out
 • arr.out
 • deg.out
 • names.out
 • func_input
 • preproc_data
 and maybe others.
 
 \param filename the name of the file to parse.
 */
void parse_input_file(boost::filesystem::path filename);

/**
 \brief splits the bertini input file into several files for later use.
 
 calls MPI_Bcast(&PARSING, 1, MPI_INT, 0, MPI_COMM_WORLD); to be able to let workers carry through.  sadly, the parse_input() method in Bertini calls an MPI_Bcast, which will trip up the workers if it is not caught.
 
 basically, the files made are
 • num.out
 • arr.out
 • deg.out
 • names.out
 • func_input
 • preproc_data
 and maybe others.
 
 \param filename the name of the file to parse.
 \param MPType a set-integer by pointer, this function splits the file and gets the MPType
 */

void parse_input_file(boost::filesystem::path filename, int * MPType);


/**
 \brief a wrapper around setupPreProcData(), and populates a preproc_data
 
 \param filename the name of the file to parse.
 \param PPD the preproc_data to populate
 */
void parse_preproc_data(boost::filesystem::path filename, preproc_data *PPD);


























#endif


