#include "curve.hpp"


void curve_main(witness_set W,
								int self_conjugate,
								program_configuration *program_options,
								solver_configuration *solve_options)
{
	
	int num_vars = W.num_variables;
	
	
	
	//get the random projection \pi
	vec_mp *pi = (vec_mp *) br_malloc(1*sizeof(vec_mp));  //random projection
	init_vec_mp(pi[0],num_vars); pi[0]->size = num_vars; // should include the homogeneous variable
	get_projection(pi, *program_options, *solve_options, num_vars, 1);
	
	
	
	
	
	//initialize the data structure which collets the output
	vertex_set V;
  curveDecomp_d C;  // stores edges, etc.
	init_curveDecomp_d(&C);
	init_vertex_set(&V);

	C.num_variables = num_vars;
	init_vec_mp(C.pi_mp, W.num_variables); C.pi_mp->size = W.num_variables;
	vec_cp_mp(C.pi_mp, pi[0]);
	
	
	
	
	solve_options->verbose_level = program_options->verbose_level;
	solve_options->T.ratioTol = 1; // manually assert to be more permissive.
	solve_options->use_gamma_trick = program_options->use_gamma_trick;
	
	if (self_conjugate==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		printf("\n\nentering not-self-conjugate case\n\n");
		computeCurveNotSelfConj(W, pi[0], &C, &V, num_vars,program_options->input_deflated_filename,
														program_options, solve_options);//This is Wenrui's !!!
		printf("Bertini_real found %d vertices (vertex)\n",C.num_V0);
	}
	else
	{
		//Call self-conjugate case code
		printf("\n\nentering self-conjugate case\n\n");
		computeCurveSelfConj(program_options->input_deflated_filename,
												 W,
												 pi,
												 &C,&V,
												 num_vars,W.num_var_gps,
												 program_options, solve_options);  //This is Dans', at least at first !!!
	}
	
	
	
	
	
	
	if (program_options->verbose_level>=2) {
		printf("outputting data\n");
	}
	Output_Main(*program_options, W, C, V);
	
	
	if (program_options->verbose_level>=2) {
		printf("clearing witness_set\n");
	}
	clear_witness_set(W);
	
	if (program_options->verbose_level>=2) {
		printf("clearing C\n");
	}
	clear_curveDecomp_d(&C);
	
	//	printf("clearing program_options\n");
	//	clear_program_config(&program_options);
	
  //TMP END


}