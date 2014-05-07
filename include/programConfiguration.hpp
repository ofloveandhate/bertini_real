
#ifndef _PROGRAM_STARTUP_H
#define _PROGRAM_STARTUP_H


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


#define BERTINI_REAL_VERSION_STRING "0.1.0"

enum {INACTIVE = 500, ACTIVE};
enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};

#include "fileops.hpp"




///////////
//
//    PROGRAM CONFIGURATION
//
//////////


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
	
	
	
	inline bool is_head()
	{
		if (my_id == headnode)
			return true;
		else
			return false;
	}
	
	
	inline bool use_parallel(){
		if (numprocs>1 && force_no_parallel!=true)
			return true;
		else
			return false;
		
		
	}
	
	
	
	
	
	void abort(int why){
		MPI_Abort(my_communicator,why);
	}
	
	inline MPI_Comm comm(){return my_communicator;}
	inline int head(){return headnode;}
	inline int id(){ return my_id;}
	inline int level(){return worker_level;}
	inline int size(){ return numprocs;}
	
	void init_active_workers()
	{

		for (int ii=1; ii<this->numprocs; ii++) {
			available_workers.push(ii);
			worker_status[ii] = INACTIVE;
		}
		
		
	}
	
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
	
	
	void deactivate(int worker_id)
	{
		if (worker_status[worker_id] == INACTIVE) {
			std::cout << "master tried decativating worker" << worker_id << " when it was already inactive" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,2);
		}
		worker_status[worker_id] = INACTIVE;
		
		
		available_workers.push(worker_id);
	}
	
	void send_all_available(int numtosend)
	{
		while (available_workers.size()>0)  {
			int sendtome = available_workers.front();
			MPI_Send(&numtosend, 1, MPI_INT, sendtome, NUMPACKETS, MPI_COMM_WORLD);
			available_workers.pop();
		}
	}
	
	
	void call_for_help(int solver_type)
	{

		MPI_Bcast(&solver_type, 1, MPI_INT, head(), MPI_COMM_WORLD);
		
		init_active_workers();
		
	}
	
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
	int debugwait;
	int stifle_membership_screen; ///< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text; //
	
	bool quick_run;
	bool user_sphere;

	bool user_projection; // bool
	
	int MPType;
	
	boost::filesystem::path bounding_sphere_filename;
	boost::filesystem::path projection_filename;
	boost::filesystem::path input_filename;

	boost::filesystem::path input_deflated_filename;
	boost::filesystem::path sphere_filename;
	

	
	int target_dimension;
	int target_component;
	
	
	std::string bertini_command;
	std::string matlab_command;
	
	bool use_gamma_trick; // bool
	
	bool merge_edges;
	
	

	
	
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
		
        init();
	};
	
	void init();
	
}; //re: BR_configuration




class sampler_configuration : public prog_config
{
public:
	int stifle_membership_screen; //< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text;
	
	int verbose_level;
	
	int maximum_num_iterations;
	
	int use_gamma_trick;
	mpf_t TOL;
	
	bool no_duplicates;
	
	bool use_fixed_sampler;
	int target_num_samples;
	
	/** get the sampler_configuration from the command line. */
	int  parse_commandline(int argc, char **argv);
	void splash_screen();
	void print_usage();
	int  parse_options(int argc, char **argv);
	
	
	sampler_configuration()
	{
		no_duplicates = true;
		use_fixed_sampler = false;
		target_num_samples = 10;
		
		
		this->stifle_membership_screen = 1;
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


