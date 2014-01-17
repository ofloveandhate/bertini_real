#include "bertini_real.hpp"

int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{


	boost::timer::auto_cpu_timer t;
    
	////
	//  INITIALIZATION
	////
	MPI_Init(&argC,&args);
	
	
	//instantiate options
	BR_configuration program_options;
	
	program_options.parse_commandline(argC, args); // everybody gets to parse the command line.
	
	
	
	
	
//	// set up the solver configuration

	
	
	// split the input_file.  this must be called before setting up the solver config.
	
	
	// set up the solver configuration
	solver_configuration solve_options; 
	
	
	if (program_options.debugwait) {
		std::cout << "in debug mode, waiting to start so you can attach to this process" << std::endl;
		if (solve_options.is_head()) {
			std::cout << "master PID: " << getpid() << std::endl;
			for (int ii=0; ii<30; ii++) {
				std::cout << "starting program in " << 30-ii << " seconds" << std::endl;
				sleep(1);
				std::cout << "\033[F"; // that's a line up movement.  only works with some terminals
			}
		}
		else{
			sleep(30);
		}
		
		
		if (solve_options.use_parallel()) {
			int arbitrary_int = 1;
			MPI_Bcast(&arbitrary_int,1,MPI_INT,0,MPI_COMM_WORLD);
		}
		//		int z = 1;
		//		while (z) ;
	}
	
	
	int MPType;
	
	if (solve_options.is_head()) {
		parse_input_file(program_options.input_filename, &MPType);
	}
	
	MPI_Bcast(&MPType, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	get_tracker_config(solve_options,MPType);
	
	initMP(solve_options.T.Precision); // set up some globals.

	solve_options.use_midpoint_checker = 0;
	solve_options.T.ratioTol = 0.99999999; // manually assert to be more permissive.  i don't really like this.
	solve_options.verbose_level = program_options.verbose_level;
	solve_options.use_gamma_trick = program_options.use_gamma_trick;
	
	if (program_options.verbose_level>=2)
		solve_options.show_status_summary = 1;
	else
		solve_options.show_status_summary = 0;
	
	
	
	if (solve_options.is_head()) {
		
		ubermaster_process current_process(program_options, solve_options);
		current_process.main_loop();
	}
	else{
		
		worker_process current_process(program_options, solve_options);
		current_process.main_loop();
	}
	
	
	

	MPI_Finalize();
	
	
	return 0;
}



