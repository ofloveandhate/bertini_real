#include "surface.h"



void surface_main(witness_set W,
									int self_conjugate,
									program_configuration *program_options,
									solver_configuration *solve_options)
{
	
	int ii;  // counters
	int ambient_dim = 2;
	
	vertex_set V;
  surface_decomposition surf; 

	vec_mp *pi = (vec_mp *) br_malloc(ambient_dim*sizeof(vec_mp ));
	for (ii=0; ii<ambient_dim; ii++) {
		init_vec_mp2(pi[ii],W.num_variables, solve_options->T.AMP_max_prec);
		pi[ii]->size = W.num_variables;
	}
	get_projection(pi, *program_options, *solve_options, W.num_variables, ambient_dim);
	
	
	
	
	//create the matrix
	mat_mp randomizer_matrix;
	init_mat_mp2(randomizer_matrix,W.num_variables-1-ambient_dim,solve_options->PPD.num_funcs,solve_options->T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)bmalloc((W.num_variables-1-ambient_dim)*sizeof(int));
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, &randomized_degrees, W.num_variables-1-ambient_dim, solve_options->PPD.num_funcs);
	
	
	
	// 2(b) compute NID of critical curve, deflating if necessary.
	witness_set W_crit_real;  init_witness_set(&W_crit_real);
	
	compute_crit_nullspace(&W_crit_real, // the returned value
												 W,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 1,   // COdimension of the critical set to find.
												 program_options,
												 solve_options);
	
	write_dehomogenized_coordinates(W_crit_real, "W_crit_real"); // write the points to file
	
	//
	
	
	// 2(c) Find singular points \emph{not} on the critical curve, and keep them in the set of vertices.
	
	
	witness_set W_crit_real2;  init_witness_set(&W_crit_real2);
	
	compute_crit_nullspace(&W_crit_real2, // the returned value
												 W,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 2,   // COdimension of the critical set to find.
												 program_options,
												 solve_options);
	
	write_dehomogenized_coordinates(W_crit_real2, "W_crit_real2"); // write the points to file
	
	
	
	
	
	
	
	for (ii=0; ii<ambient_dim; ii++) 
		clear_vec_mp(pi[ii]);
	
}


