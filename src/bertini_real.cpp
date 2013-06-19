#include "bertini_real.hpp"

int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	int rV,num_vars=0,self_conjugate;  //1=self-conjugate; 0=not
	
	
	witness_set W;
	int max_deflations = 10;
	int num_deflations, *deflation_sequence = NULL;
	int MPType;
	
	
	////
	//  INITIALIZATION
	////
	
	BR_splash_screen();
	
	//instantiate options
	program_configuration program_options;  BR_init_config(&program_options);
	BR_parse_commandline(argC, args, &program_options);
	
	//parse the options
	BR_startup(program_options); // tests for existence of necessary files, etc.
	
	//if desired, display the options
	if (program_options.verbose_level>=3)
		BR_display_current_options(program_options);
	
	
	srand(time(NULL));
	
	
	// split the input_file.  this must be called before setting up the solver config.
	parse_input_file(program_options.input_filename, &MPType);
	
	
	// set up the solver configuration
	solver_configuration solve_options;  solver_init_config(&solve_options);
	get_tracker_config(&solve_options,MPType);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	initMP(solve_options.T.Precision); // set up some globals.
	
	
	num_vars = get_num_vars_PPD(solve_options.PPD);
	
	
	init_witness_set(&W);
	witnessSetParse(&W,program_options.witness_set_filename,num_vars);
	W.num_var_gps = solve_options.PPD.num_var_gp;
	W.MPType = program_options.MPType = MPType;
	
	get_variable_names(&W);
	
	
	W.incidence_number = get_incidence_number(W,num_vars,
											  program_options.input_filename,
											  program_options.stifle_text);
	
	
	if (program_options.verbose_level>=2) {
		printf("performing isosingular deflation\n");
	}
	
	// perform an isosingular deflation
	write_dehomogenized_coordinates(W, "witness_points_dehomogenized"); // write the points to file
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, program_options.input_filename, "witness_points_dehomogenized", "bertini", "matlab -nosplash", max_deflations, W.dim, W.comp_num);
	

	
	program_options.input_deflated_filename = program_options.input_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W.dim << "_comp_" << W.comp_num << "_deflated";
	program_options.input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename, &MPType);
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	
	
	if (program_options.verbose_level>=2) {
		printf("checking if component is self-conjugate\n");
	}
	self_conjugate = checkSelfConjugate(W,num_vars,program_options.input_filename, program_options.stifle_text);  //later:  could be passed in from user, if we want
	
	
	
	
	//regenerate the various files, since we ran bertini since then.
	parse_input_file(program_options.input_deflated_filename, &MPType);
	
	
	
	switch (W.dim) {
		case 1:
			//curve
			curve_main(W, self_conjugate, &program_options, &solve_options);
			break;
			
		case 2:
			surface_main(W, self_conjugate, &program_options, &solve_options);
			break;
		default:
			break;
	}
	return 0;
}



