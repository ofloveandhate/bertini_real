#include "curve.hpp"


void curve_main(vertex_set & V,
								curve_decomposition & C,
								witness_set & W,
								vec_mp *pi,
								BR_configuration & program_options,
								solver_configuration & solve_options)
{
	
	int num_vars = W.num_variables;
	//HERE CHANGE THIS
	
	C.num_variables = num_vars;
	C.add_projection(pi[0]);


	// perform an isosingular deflation
	// note: you do not need witness_data to perform isosingular deflation
	if (program_options.verbose_level>=2)
		printf("performing isosingular deflation\n");
	
	W.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
	int num_deflations, *deflation_sequence = NULL;
	
	
	isosingular_deflation(&num_deflations, &deflation_sequence, program_options, program_options.current_working_filename,
																 "witness_points_dehomogenized", program_options.max_deflations, W.dim, W.comp_num);
	free(deflation_sequence);
	
	
	program_options.input_deflated_filename = program_options.current_working_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W.dim << "_comp_" << W.comp_num << "_deflated";
	program_options.input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	std::cout << "!@#$ before the parse command\n";
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename);
	std::cout << "%^&* after the parse command\n";
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	int self_conjugate = 1;
	if (W.num_synth_vars==0) {

		if (program_options.verbose_level>=2) {
			printf("checking if component is self-conjugate\n");
		}
		self_conjugate = checkSelfConjugate(W,num_vars,program_options, program_options.current_working_filename);  //later:  could be passed in from user, if we want
		
		//regenerate the various files, since we ran bertini since then.
		parse_input_file(program_options.input_deflated_filename);

		if (verify_projection_ok(W,
														 *pi,
														 solve_options)==1){
			if (program_options.verbose_level>=1) {
				printf("verified projection is ok\n");
			}
		}
		else{
			printf("the projection is invalid, in that the jacobian of the randomized system\nbecomes singular at a random point, when the projection is concatenated\n");
			
			print_point_to_screen_matlab(pi[0], "pi[0]");
			
			exit(196);
		}
		
	}
	
	

	
	

	
	if (self_conjugate==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		printf("\n\nentering not-self-conjugate case\n\n");
		computeCurveNotSelfConj(W, pi[0], C, V, num_vars,program_options.input_deflated_filename,
														program_options, solve_options);//This is Wenrui's !!!
		printf("Bertini_real found %d vertices (vertex)\n",C.counters[ISOLATED]);
	}
	else
	{
		//Call self-conjugate case code
		printf("\n\nentering self-conjugate case\n\n");
		computeCurveSelfConj(program_options.input_deflated_filename,
												 W,
												 pi,
												 C,V,
												 num_vars,
												 program_options, solve_options);  //This is Dans', at least at first !!!
	}
	

	
}






