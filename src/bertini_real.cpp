#include "bertini_real.hpp"

int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{


	
	
	
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
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}
	
	
	int MPType;
	
	if (solve_options.is_head()) {
		parse_input_file(program_options.input_filename, &MPType);
		get_tracker_config(solve_options,MPType);
		
	}
	else
	{ // catch the bcast from parallel parsing. (which cannot be disabled)
		int arbitrary_int;
		MPI_Bcast(&arbitrary_int,1,MPI_INT,0,MPI_COMM_WORLD);
	}
	
	
	
	if (solve_options.use_parallel()) {
		MPI_Bcast(&MPType, 1, MPI_INT, 0, MPI_COMM_WORLD);
		bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	}
	
	
	initMP(solve_options.T.Precision); // set up some globals.
	
	
	solve_options.use_midpoint_checker = 0;
	solve_options.T.ratioTol = 0.9999999999999999999999999; // manually assert to be more permissive.  i don't really like this.
	solve_options.verbose_level = program_options.verbose_level;
	solve_options.use_gamma_trick = program_options.use_gamma_trick;
	

	
	
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



