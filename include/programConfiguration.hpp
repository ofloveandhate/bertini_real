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


enum {BERTINIREAL=-9000,CRIT=-8999};

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
	
	bool force_no_parallel_;
	int headnode_;
	int my_id_, my_id_global_;
	int numprocs_;
	MPI_Comm   my_communicator_;
	
	int worker_level_; // higher worker level means more tedious work, in a sense.  worker_level 0 is uber-master.  worker_level 1 will be the next level down in management, etc.  the exact usage of this is relative to the process being run.
	
	
	
	std::map< int, int> worker_status_;
	
	std::queue< int > available_workers_;
	
	
	
public:
	
	/**
	 \brief get the state of whether forcing no parallel at the current time
	 \return indicator of whether parallelism currently turned off
	 */
	bool force_no_parallel()
	{
		return force_no_parallel_;
	}
	
	
	/**
	 \brief set the value of force_no_parallel, to turn on/off parallelism
	 
	 \param new_val the new value to set it to
	 */
	void force_no_parallel(bool new_val)
	{
		force_no_parallel_ = new_val;
	}
	
	
	
	
	parallelism_config(){
		init();
	}
	
	
	/**
	 \brief determine whether the process is the head.
	 
	 \return Indicator of whether the process thinks it is currently the head.
	 */
	inline bool is_head()
	{
		if (my_id_ == headnode_)
			return true;
		else
			return false;
	}
	
	
	/**
	 
	 \brief Determine whether the process thinks it should be in parallel mode.
	 
	 \return Indicator of whether to use parallel mode.
	 */
	inline bool use_parallel(){
		if (numprocs_>1 && force_no_parallel_!=true)
			return true;
		else
			return false;
		
		
	}
	
	
	
	
	/**
	 \brief Call MPI_Abort using the internally stored communicator.
	 
	 \param why An integer to feed to MPI_Abort.
	 */
	void abort(int why){
		MPI_Abort(my_communicator_,why);
	}
	
	/**
	 Get the current communicator
	 \return the currently stored communicator
	 */
	inline MPI_Comm comm(){return my_communicator_;}
	
	/**
	 \brief Get the ID of the supervisor
	 
	 \return the ID of the head node, supervisor, or whatever you want to call it.
	 */
	inline int head(){return headnode_;}
	
	/**
	 \brief Get the ID of this process
	 
	 \return the ID
	 */
	inline int id(){ return my_id_;}
	
	/**
	 \brief  Get the worker level
	 
	 \return the worker level
	 */
	inline int level(){return worker_level_;}
	
	/**
	 \brief Get how many workers there are in the current communicator
	 
	 \return The number of workers.
	 */
	inline int num_procs(){ return numprocs_;}
	
	
	
	
	/**
	 \brief Get how many workers there are in the current communicator
	 
	 \return The number of workers.
	 */
	inline int size(){ return numprocs_;}
	
	
	/**
	 \brief Set up the vector of available workers, based on how many processors there are.
	 \throws logic_error if the set of active workers is not empty.
	 
	 It also sets up the workers to be listed as inactive.
	 */
	void init_active_workers()
	{
		if (!available_workers_.empty()) {
			throw std::logic_error("set of available workers is not empty at call of init_active_workers...  it must be empty.  some previous process did not finish properly, dismissing all workers at the end.");
		}
		
		while (!available_workers_.empty())
			available_workers_.pop();
		
		
		for (int ii=1; ii<this->numprocs_; ii++) {
			available_workers_.push(ii);
			worker_status_[ii] = INACTIVE;
		}
		
		
	}
	
	/**
	 \brief Get the next available worker, relist it as active, and return its id.
	 
	 Available workers are stored as a queue of integers, and work is assigned to the front of the vector.  This method pops the front entry of the queue, relists it as active, and returns its ID.
	 \return the ID of the next worker.
	 */
	int activate_next_worker()
	{
		int worker_id = available_workers_.front();
		available_workers_.pop();
		
		if (worker_status_[worker_id] == ACTIVE) {
			std::cout << "master tried making worker" << worker_id << " active when it was already active" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		worker_status_[worker_id] = ACTIVE;
		

		return worker_id;
	}
	
	/**
	 \brief Take an active worker, relist it as inactive, and make it available.  If the worker is not active, it calls MPI_Abort
	 
	 \param worker_id The id of the worker to deactivate.
	 */
	void deactivate(int worker_id)
	{
		if (worker_status_[worker_id] == INACTIVE) {
			std::cout << "master tried decativating worker" << worker_id << " when it was already inactive" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,2);
		}
		worker_status_[worker_id] = INACTIVE;
		
		
		available_workers_.push(worker_id);
	}
	
	
	/**
	 \brief Individually send a number to all available workers.
	 
	 pops available workers off the queue, and deliveres the them numtosend.
	 
	 \param numtosend The integer number to send.
	 */
	void send_all_available(int numtosend)
	{
		while (available_workers_.size()>0)  {
			int sendtome = available_workers_.front();
			MPI_Send(&numtosend, 1, MPI_INT, sendtome, NUMPACKETS, MPI_COMM_WORLD);
			available_workers_.pop();
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
		if (available_workers_.size()==0) {
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
		for (int ii=1; ii<this->numprocs_; ii++) {
			if (this->worker_status_[ii]==ACTIVE) {
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
		for (int ii=1; ii<this->numprocs_; ii++) {
			if (this->worker_status_[ii]==ACTIVE) {
				num++;
			}
		}
		return num;
	}
	
	
private:
	
	void init()
	{
		
		force_no_parallel_ = false;
		numprocs_ = 1;
		headnode_ = 0;
		
		MPI_Comm_size(MPI_COMM_WORLD, &this->numprocs_);
		MPI_Comm_rank(MPI_COMM_WORLD, &this->my_id_);
		
		
		if (is_head())
			worker_level_ = 0;
		else
			worker_level_ = 1;
		
		
		
		my_communicator_ = MPI_COMM_WORLD; // default communicator is MPI_COMM_WORLD
		
        
        
        
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
	
	
	int verbose_level_;
	
	boost::filesystem::path called_dir_;
	boost::filesystem::path working_dir_;
	boost::filesystem::path output_dir_;
	
protected:
	
	/**
	 \brief set the name of the called directory
	 \param new_name the new name of the called directory.
	 */
	void set_called_dir(boost::filesystem::path new_name)
	{
		called_dir_ = new_name;
	}
	
public:
	
	/**
	 get the directory from which the program was called
	 
	 \return the directory in which the user was, when the prog_config was created.
	 */
	boost::filesystem::path called_dir() const
	{
		return called_dir_;
	}
	
	
	/**
	 get the working directory
	 
	 \return the working directory
	 */
	boost::filesystem::path working_dir() const
	{
		return working_dir_;
	}
	
	
	/**
	 \brief set the name of the working directory
	 
	 \param new_name the new name of the working directory
	 */
	void working_dir(boost::filesystem::path new_name)
	{
		working_dir_ = new_name;
	}
	
	
	
	/**
	 get the output directory
	 
	 \return the output directory
	 */
	boost::filesystem::path output_dir() const
	{
		return output_dir_;
	}
	
	
	/**
	 \brief set the name of the output directory
	 
	 \param new_name the new name of the output directory
	 */
	void output_dir(boost::filesystem::path new_name)
	{
		output_dir_ = new_name;
	}
	
	
	
	/**
	 \brief get the level of verbosity
	 
	 \return the level of verbosity
	 */
	inline int verbose_level() const
	{
		return verbose_level_;
	}
	
	/**
	 \brief set the level of verbosity
	 
	 \param new_level the new level of verbosity
	 */
	int verbose_level(int new_level)
	{
		return verbose_level_ = new_level;
	}
	
	
	void move_to_temp();
	
	void move_to_called();
	
};



/**
 \brief holds the current state of configuration for Bertini_real.
 
 */
class BR_configuration : public prog_config
{
	bool orthogonal_projection_;
	
	bool debugwait_; ///< flag for whether to wait 30 seconds before starting, and print the master process ID to screen.
	int max_deflations_; ///< the maximum allowable number of deflation iterations before it gives up.
	
	bool stifle_membership_screen_; ///< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text_; ///<
	
	int quick_run_;  ///< indicator of whether to use the robust solver wherever possible
	
	
	bool user_sphere_; ///< flag for whether to read the sphere from a file, rather than compute it.
	
	bool user_projection_; ///< indicator for whether to read the projection from a file, rather than randomly choose it.
	
	bool merge_edges_; ///< a mode switch, indicates whether should be merging.
public:
	
	/** 
	 \brief query whether should merge edges
	 \return whether we should merge edges or not
	 */
	bool merge_edges()
	{
		return merge_edges_;
	}
	
	
	/**
	 \brief set whether should merge edges
	 \param new_val set whether we should merge edges or not
	 */
	void merge_edges(bool new_val)
	{
		merge_edges_ = new_val;
	}
	
	int primary_mode; ///< mode of operation -- bertini_real is default, but there is also crit method for computing critical points.
	
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
	
	
	
	
	
	
	
	
	/**
	 \brief whether to read projection from file.
	 \return whether to read the projection from a user-created and specified file.
	 */
	bool user_projection()
	{
		return user_projection_;
	}
	
	/**
	 \brief set whether to read projection from file.
	 \param new_val whether to read the projection from a user-created and specified file.
	 */
	void user_projection(bool new_val)
	{
		user_projection_ = new_val;
	}
	
	
	/**
	 \brief whether to read sphere parameters from file.
	 \return whether to read the sphere parameters from a user-created and specified file.
	 */
	bool user_sphere()
	{
		return user_sphere_;
	}
	
	
	
	/**
	 \brief set whether to read sphere parameters from file.
	 \param new_val whether to read the sphere parameters from a user-created and specified file.
	 */
	void user_sphere(bool new_val)
	{
		user_sphere_ = new_val;
	}
	
	/**
	 \brief get the level of quick.  higher == faster (less robust)
	 \return the level of quickness
	 */
	int quick_run()
	{
		return quick_run_;
	}
	
	/**
	 \brief set the level of quickness
	 \param new_val the new level
	 */
	void quick_run(int new_val)
	{
		quick_run_ = new_val;
	}
	
	
	
	
	/**
	 \brief the stifling text for system commands
	 \return string to append to system commands to stifle screen output.
	 */
	std::string stifle_text() const
	{
		return stifle_text_;
	}
	
	/**
	 \brief set the stifling text
	 \param new_val the new stifling text
	 */
	void stifle_text(std::string new_val)
	{
		stifle_text_ = new_val;
	}
	
	
	
	/**
	 \brief should we wait for 30 seconds before running program?
	 \return whether we should.  true==yes
	 */
	bool debugwait() const
	{
		return debugwait_;
	}
	
	/**
	 \brief set value for question -- should we wait for 30 seconds before running program?
	 \param new_val the new value for the debugwait parameter
	 */
	void debugwait(bool new_val)
	{
		debugwait_ = new_val;
	}
	
	
	/**
	 how many times it is ok to deflate
	 \return how many times to deflate, at maximum.  default is 10.
	 */
	int max_deflations() const
	{
		return max_deflations_;
	}
	
	/**
	 \brief set the maximum number of deflations
	 */
	void max_deflations(int new_val)
	{
		max_deflations_ = new_val;
	}
	
	
	/**
	 \brief query whether we should stifle some screen output
	 \return whether we should.
	 */
	bool stifle_membership_screen() const
	{
		return stifle_membership_screen_;
	}
	
	/**
	 \brief set whether should stifle the membership testing screen output
	 */
	void stifle_membership_screen(bool new_val)
	{
		stifle_membership_screen_ = new_val;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 get whether should use orthogonal projection.  default is yes
	 */
	bool orthogonal_projection()
	{
		return orthogonal_projection_;
	}
	
	
	
	
	
	

	
	
	
	
	
	
	
	
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


